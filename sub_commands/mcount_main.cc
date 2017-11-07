/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <unistd.h>
#include <assert.h>
#include <signal.h>
#include <mpi.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <memory>
#include <chrono>

#include <jellyfish/err.hpp>
#include <jellyfish/mer_heap.hpp>
#include <jellyfish/jellyfish.hpp>
#include <jellyfish/merge_files.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <sub_commands/count_main_cmdline.hpp>

#include "mimir.h"
using namespace MIMIR_NS;

#ifdef VAL64BITS
using CountType = uint64_t;
#else
using CountType = uint32_t;
#endif

uint64_t uniq = 0, distinct = 0, total = 0, maxcount = 0;

static count_main_cmdline args; // Command line switches and arguments

namespace err = jellyfish::err;

using std::chrono::system_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

template<typename DtnType>
inline double as_seconds(DtnType dtn) { return duration_cast<duration<double>>(dtn).count(); }

using jellyfish::mer_dna;
using jellyfish::mer_dna_bloom_counter;
using jellyfish::mer_dna_bloom_filter;
typedef std::vector<const char*> file_vector;

int rank, size;

enum OPERATION { COUNT, PRIME, UPDATE };

enum FileType { NULL_TYPE, FASTA_TYPE, FASTQ_TYPE };

struct MimirParam{
    struct filter *filter;
    int  key_len, val_len;
    mer_dna m, rcm;
    unsigned int filled;
    FileType ftype;
};

MimirParam param;

// k-mer filters. Organized in a linked list, interpreted as a &&
// (logical and). I.e. all filter must return true for the result to
// be true. By default, filter returns true.
struct filter {
  filter* prev_;
  filter(filter* prev = 0) : prev_(prev) { }
  virtual ~filter() { }
  virtual bool operator()(const mer_dna& x) { return and_res(true, x); }
  bool and_res(bool r, const mer_dna& x) const {
    return r ? (prev_ ? (*prev_)(x) : true) : false;
  }
};

struct filter_bc : public filter {
  const mer_dna_bloom_counter& counter_;
  filter_bc(const mer_dna_bloom_counter& counter, filter* prev = 0) :
    filter(prev),
    counter_(counter)
  { }
  bool operator()(const mer_dna& m) {
    unsigned int c = counter_.check(m);
    return and_res(c > 1, m);
  }
};

struct filter_bf : public filter {
  mer_dna_bloom_filter& bf_;
  filter_bf(mer_dna_bloom_filter& bf, filter* prev = 0) :
    filter(prev),
    bf_(bf)
  { }
  bool operator()(const mer_dna& m) {
    unsigned int c = bf_.insert(m);
    return and_res(c > 0, m);
  }
};

void report_error(const char *msg) {
    fprintf(stderr, "%d[%d] Error: %s\n", rank, size, msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
}

// add a mer to hash array
void add_mer(mer_dna &m, uint64_t v,
             mer_hash& ary, filter& filter, OPERATION op)
{
    if(filter(m)){
        switch(op){
        case COUNT: {
            ary.add(m, v); break;
        }
        case PRIME: ary.set(m); break;
        case UPDATE: {
            mer_dna tmp;
            ary.update_add(m, v);
        }; break;
        }
    }
}

class MerDatabase : public BaseDatabase<char, CountType>
{
  public:
    typedef typename mer_array::key_type key_type;
    typedef typename mer_array::region_iterator iterator;
    typedef typename jellyfish::mer_heap::heap<typename mer_array::key_type, iterator> heap_type;
    typedef typename heap_type::const_item_t    heap_item;

    uint64_t get_record_count() {
        return 0;
    }

    int key_len() {
        return ary->key_len();
    }

    MerDatabase(OPERATION op,
                filter* filter = new struct filter) {
        this->op = op;
        this->filter = filter;
        ary = new mer_hash(args.size_arg,
                           args.mer_len_arg * 2,
                           args.counter_len_arg, 1,
                           args.reprobes_arg);
        ary->do_size_doubling(true);
        min_val = 0;
        max_val = std::numeric_limits<uint64_t>::max();
        key_len_ = ary->ary()->key_len() / 8 + (ary->ary()->key_len() % 8 != 0);
        val_len_ = args.out_counter_len_arg;
        max_val_ = (((uint64_t)1 << (8 * val_len_)) - 1);
        block_info = ary->ary()->blocks_for_records(5 * ary->ary()->max_reprobe_offset());
        heap = new heap_type(ary->ary()->max_reprobe_offset());
    }

    ~MerDatabase() {
        delete ary;
        delete heap;
    }

    int open() {
        id = 0;
        it = NULL;
        if (id * block_info.second < ary->ary()->size()) {
            it = new iterator(ary->ary(),
                              id * block_info.second,
                              (id + 1) * block_info.second,
                              key);
            heap->fill((*it));
        }
        return true;
    }

    void close() {
        if (it != NULL) {
            delete it;
            it = NULL;
        }
    }

    void analysis() {
        uint64_t uniq = 0, distinct = 0, total = 0, maxcount = 0;
        uint64_t global_uniq = 0, global_distinct = 0, global_total = 0, global_max = 0;

        for (id = 0; id * block_info.second < ary->ary()->size(); id ++) {
            it = new iterator(ary->ary(),
                              id * block_info.second,
                              (id + 1) * block_info.second,
                              key);
            heap->fill((*it));
            while (heap->is_not_empty()) {
                heap_item item = heap->head();
                if(item->val_ >= min_val && item->val_ <= max_val) {
                    uniq  += item->val_ == 1;
                    total += item->val_;
                    maxcount    = std::max(maxcount, item->val_);
                    ++distinct;
                }
                heap->pop();
                if(it->next()) heap->push(*it);
            }
            delete it;
        }

        MPI_Reduce(&uniq, &global_uniq, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&distinct, &global_distinct, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total, &global_total, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&maxcount, &global_max, 1, MPI_INT64_T, MPI_MAX, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            std::cout << "Unique:    " << global_uniq << std::endl
                << "Distinct:  " << global_distinct << std::endl
                << "Total:     " << global_total << std::endl
                << "Max_count: " << global_max << std::endl;
        }
    }

    int read(char *mr_key, CountType *mr_val) {

        while (id * block_info.second < ary->ary()->size()) {
            while (heap->is_not_empty()) {
                heap_item item = heap->head();
                if(item->val_ >= min_val && item->val_ <= max_val) {
                    uint64_t v = std::min(max_val_, item->val_);
                    heap->pop();
                    if(it->next()) heap->push(*it);
                    memcpy(mr_key, (char*)(item->key_).data(), key_len_);
                    *mr_val = v;
                    return true;
                }
                heap->pop();
                if(it->next()) heap->push(*it);
            }
            id += 1;
            if (it != NULL) delete it;
            it = new iterator(ary->ary(),
                              id * block_info.second,
                              (id + 1) * block_info.second,
                              key);
            heap->fill(*it);
        };

        return false;
    }

    int write(char *mr_key, CountType *mr_val) {

        mer_dna mer(mer_dna::k());
        memcpy(mer.data__(), mr_key, param.key_len);
        CountType count = *mr_val;

        if((*filter)(mer)){
            switch(op){
                case COUNT: {
                    ary->add(mer, count); break;
                 }
                case PRIME: ary->set(mer); break;
                case UPDATE: {
                    ary->update_add(mer, count);
                }; break;
            }
        }

        return true;
    }

    int remove() {
        fprintf(stderr, "Donot support remove funtion!\n");
        exit(1);
        return false;
    }

    int seek(DB_POS pos) {
        fprintf(stderr, "Donot support seek funtion!\n");
        exit(1);
        return false;
    }

    void set_min_val(CountType min) {
        this->min_val = min;
    }

    void set_max_val(CountType max) {
        this->max_val = max;
    }

  private:
    mer_hash* ary;
    struct filter* filter;
    OPERATION op;

    iterator *it;
    heap_type*  heap;
    typename mer_array::key_type key;
    size_t id;

    std::pair<size_t, size_t> block_info;
    uint64_t min_val, max_val;
    int key_len_, val_len_;
    uint64_t max_val_;
};

extern mer_dna_bloom_counter* load_bloom_filter(const char* path);

int fasta_padding(const char *buffer, int len) {
    int padding = std::min((unsigned int)len, mer_dna::k() - 1);
    param.filled = 0;
    for (int j = 0; j < padding; j++) {
        char ch = buffer[j];
        int code = param.m.code(ch);
        if (code >= 0) {
            param.m.shift_left(code);
            if (args.canonical_flag)
                param.rcm.shift_right(param.rcm.complement(code));
            param.filled = std::min(param.filled + 1, mer_dna::k());
        }
        else {
            param.filled = 0;
        }
    }
    return padding;
}

int repartition(uint64_t foff, const char *buffer, int len, bool islast) {
    int padding = 0;
    int pre_pos = 0;
    char pre = 0;
    int pos = 0;
    int line_num = 5;

    for (int i=0; i < line_num; i++) {
        while (pos < len && buffer[pos] != '\n') pos ++;
        if (pos == len && !islast) {
            report_error("Cannot judge the border.");
        }
        if (pos == len || pos == (len - 1))
            return (int)len;
        pos ++;

        if (pre == '>') {
            if (buffer[pos] == '@') {
                param.ftype = FASTQ_TYPE;
                return pos;
            }
            else {
                param.ftype = FASTA_TYPE;
                return pre_pos;
            }
        }

        if (pre == '@') {
            param.ftype = FASTQ_TYPE;
            if (buffer[pos] == '@') return pos;
            else return pre_pos;
        }

        if (buffer[pos] == '@' || buffer[pos] == '>') {
            pre_pos = pos;
            pre = buffer[pos];
        }
    }

    param.ftype = FASTA_TYPE;

    pos = 0;
    while (buffer[pos] != '\n') pos ++;
    pos ++;
    padding = fasta_padding(buffer + pos, (int)len - pos);
    return pos + padding;
}

void add_mers(const char *seq,
              Writable<char,CountType> *output, void *ptr)
{
    int len = (int)strlen(seq);
    for (int i = 0; i < len; i++) {
        char ch = seq[i];
        int code = param.m.code(ch);
        if (code >= 0) {
            param.m.shift_left(code);
            if (args.canonical_flag)
                param.rcm.shift_right(param.rcm.complement(code));
            param.filled = std::min(param.filled + 1, mer_dna::k());
            if (param.filled >= mer_dna::k()) {
                mer_dna mer = !args.canonical_flag 
                    || param.m < param.rcm ? param.m : param.rcm;
                if((*(param.filter))(mer)) {
                    CountType one = 1;
                    output->write((char*)mer.data(), &one);
                }
            }
        } else {
            param.filled = 0;
        }
    }

}

void parse_sequence(Readable<char*,void> *input,
                    Writable<char,CountType> *output, void *ptr)
{
    int ret = 0;
    char *seq = NULL;
    char *line = NULL;

    while (input->read(&line, NULL) == true) {
        //printf("line=%s\n", line);
        if (line[0] == '>'){
            param.filled = 0;
            param.ftype = FASTA_TYPE;
            continue;
        } else if (line[0] == '@') {
            param.ftype = FASTQ_TYPE;
            ret = input->read(&line, NULL);
            if (ret != true) return;
            seq = line;
            add_mers(seq, output, ptr);
            // '+' line
            ret = input->read(&line, NULL);
            if (ret != true) report_error("Fastq file format incorrect!!!");
            // score line
            ret = input->read(&line, NULL);
            if (ret != true) report_error("Fastq file format incorrect!!!");
            param.filled = 0;
        } else if (param.ftype == FASTA_TYPE) {
            seq = line;
            add_mers(seq, output, ptr);
        } else {
            report_error("File format incorrect!!!");
        }

    }
}

void add_qual_mers(std::string &seq, std::string &qual,
                   Writable<char,CountType> *output, void *ptr)
{
    if (seq.size() != qual.size()) {
        report_error("The quality score does not mactch sequence length!");
    }

    std::string::iterator seq_iter = seq.begin();
    std::string::iterator qual_iter = qual.begin();
    for (;seq_iter != seq.end(); seq_iter++, qual_iter++) {
        const int code = param.m.code(*seq_iter);
        const char qual = *qual_iter;
        if (code >=0 && qual >= args.min_qual_char_arg[0]) {
            param.m.shift_left(code);
            if (args.canonical_flag)
                param.rcm.shift_right(param.rcm.complement(code));
            param.filled = std::min(param.filled + 1, mer_dna::k());
            if (param.filled >= param.m.k()) {
                mer_dna mer = !args.canonical_flag || param.m < param.rcm ? param.m : param.rcm;
                if((*(param.filter))(mer)) {
                    CountType one = 1;
                    output->write((char*)mer.data(), &one);
                }
            }
        }else{
            param.filled = 0;
        }
    }

    seq.clear();
    qual.clear();
}

void parse_qual_sequence(Readable<char*, void> *input,
                         Writable<char, CountType> *output, void *ptr)
{
    char *line = NULL;
    std::string seq, qual;
    int ret = 0;

    while (input->read(&line, NULL) == true) {

        if (line[0] == '@') {
            ret = input->read(&line, NULL);
            if (ret != true) return;
            // sequence line
            seq += line;
            // '+' line
            ret = input->read(&line, NULL);
            if (ret != true) report_error("Fastq file format incorrect!!!");
            // score line
            ret = input->read(&line, NULL);
            if (ret != true) report_error("Fastq file format incorrect!!!");
            qual += line;
            param.filled = 0;
        } else {
            report_error("Fastq file format incorrect!!!");
        }

        add_qual_mers(seq, qual, output, ptr);
    }
}

void combine_fn(Combinable<char,CountType> *combiner,
                char* key, CountType* val1, CountType* val2,
                CountType* val3, void* ptr) {
    *val3 = *val1 + *val2;
}

void analysis_fn(char *key, CountType *val, void *ptr) {
    uniq ++;
    if (*val == 1) distinct ++;
    total += *val;
    if (*val > maxcount) maxcount = *val;
}

void analysis() {
    uint64_t global_uniq = 0, global_distinct = 0, global_total = 0, global_max = 0;

    MPI_Reduce(&uniq, &global_uniq, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&distinct, &global_distinct, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total, &global_total, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&maxcount, &global_max, 1, MPI_INT64_T, MPI_MAX, 0, MPI_COMM_WORLD);

     if (rank == 0) {
        std::cout << "Unique:    " << global_uniq << std::endl
            << "Distinct:  " << global_distinct << std::endl
            << "Total:     " << global_total << std::endl
            << "Max_count: " << global_max << std::endl;
    }
}

void output_txt(Readable<char, CountType> *input,
                Writable<const char*, CountType> *output, void *ptr) {
    char key[param.key_len];
    CountType count;
    std::string key_txt;
    const char *key_ptr = NULL;
    while (input->read(key, &count) == true) {
        mer_dna mer(mer_dna::k());
        memcpy(mer.data__(), key, param.key_len);
        key_txt = mer.to_str();
        key_ptr = key_txt.c_str();
        output->write(&key_ptr, &count);
    }
}

int mcount_main(int argc, char *argv[])
{
  auto start_time = system_clock::now();

  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  jellyfish::file_header header;
  header.fill_standard();
  header.set_cmdline(argc, argv);

  args.parse(argc, argv);

  if(args.min_qual_char_given && args.min_qual_char_arg.size() != 1)
    count_main_cmdline::error("[-Q, --min-qual-char] must be one character.");

  mer_dna::k(args.mer_len_arg);

  if(args.generator_given) {
    fprintf(stderr, "mcount currently does not support -g, --generator=path parameter.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if(args.threads_arg > 1){
    fprintf(stdout, "mcount currently does not support multithreading (1 therad is used)\n");
  }

  header.canonical(args.canonical_flag);
  if(args.disk_flag){
    //ary.do_size_doubling(false);
    fprintf(stderr, "mcount currently does not support --disk parameter.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  auto after_init_time = system_clock::now();

  OPERATION do_op = COUNT;
  if(args.if_given) {
    fprintf(stderr, "mcount currently does not support --if=path parameter.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // Bloom counter read from file to filter out low frequency
  // k-mers. Two pass algorithm.
  std::unique_ptr<filter> mer_filter(new filter);
  std::unique_ptr<mer_dna_bloom_counter> bc;
  if(args.bc_given) {
    bc.reset(load_bloom_filter(args.bc_arg));
    mer_filter.reset(new filter_bc(*bc));
  }

  // Bloom filter to filter out low frequency k-mers. One pass
  // algorithm.
  std::unique_ptr<mer_dna_bloom_filter> bf;
  if(args.bf_size_given) {
    bf.reset(new mer_dna_bloom_filter(args.bf_fp_arg, args.bf_size_arg));
    mer_filter.reset(new filter_bf(*bf));
  }


  if(rank == 0) printf("start reading files\n");
  std::vector<std::string> input_dirs;
  for (file_vector::const_iterator iter = args.file_arg.begin(); 
       iter != args.file_arg.end(); iter++) {
      input_dirs.push_back((*iter));
  }

  int key_len = args.mer_len_arg * 2;
  param.filter = mer_filter.get();
  param.key_len = key_len / 8 + (key_len % 8 != 0);
  param.val_len = sizeof(CountType);
  param.filled  = 0;

  void (*map_fn)(Readable<char*, void> *,
                 Writable<char, CountType> *, void *) = NULL;

  if (args.min_qual_char_given)
      map_fn = parse_qual_sequence;
  else
      map_fn = parse_sequence;


  MimirContext<char, CountType, char*, void> *ctx =
      new MimirContext<char, CountType, char*, void>(
       input_dirs, std::string(args.output_arg),
       "text", "text",
       MPI_COMM_WORLD,
#ifndef USE_MIMIR_DB
       NULL,
#else
       combine_fn,
#endif
       NULL,
       repartition,
       param.key_len, 1, 1, 0, param.key_len, 1);
#ifndef USE_MIMIR_DB
  MerDatabase db(do_op, mer_filter.get());
  ctx->set_user_database(&db);
#endif
  ctx->map(map_fn);
  ctx->scan(analysis_fn);

  //if (args.text_flag) {
      MimirContext<const char*, CountType, char, CountType> *out_ctx =
          new MimirContext<const char*, CountType, char, CountType>(
                std::vector<std::string>(), std::string(args.output_arg),
                "null", "text",
                MPI_COMM_WORLD, NULL, NULL, NULL,
                1, 1, param.key_len, 1, 1, 1);
      out_ctx->insert_data_handle(ctx->get_data_handle());
      out_ctx->map(output_txt, NULL, false, true);
      delete out_ctx;
  //} else {
  //    ctx->output();
  //}
  delete ctx;

  auto after_count_time = system_clock::now();

  analysis();

  MPI_Finalize();

  auto after_dump_time = system_clock::now();

  if(args.timing_given) {
    std::ofstream timing_file(args.timing_arg);
    timing_file << "Init     " << as_seconds(after_init_time - start_time) << "\n"
                << "Counting " << as_seconds(after_count_time - after_init_time) << "\n"
                << "Writing  " << as_seconds(after_dump_time - after_count_time) << "\n";
  }

  return 0;
}
