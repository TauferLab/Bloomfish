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

#define MAX_ITEM_SIZE 1048576

#define COMBINE

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
            //printf("rank=%d, %s, %ld\n", rank, m.to_str().c_str(), v);
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

template<typename storage_t>
class MerDatabase : public BaseDatabase
{
  public:
    typedef typename storage_t::key_type key_type;
    typedef typename storage_t::region_iterator iterator;
    typedef typename jellyfish::mer_heap::heap<typename storage_t::key_type, iterator> heap_type;
    typedef typename heap_type::const_item_t    heap_item;

    uint64_t get_record_count() {
        return 0;
    }

    int key_len() {
        return ary->key_len();
    }

    MerDatabase(OPERATION op, filter* filter = new struct filter) {
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
    }

    ~MerDatabase() {
        delete ary;
    }

    bool open() {
        id = 0;
        it = NULL;
        buflen = 0;
        heap = new heap_type(ary->ary()->max_reprobe_offset());
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
        delete heap;
        heap = NULL;
    }

    BaseRecordFormat* read() {

        while (id * block_info.second < ary->ary()->size()) {
            while (heap->is_not_empty()) {
                heap_item item = heap->head();
                if(item->val_ >= min_val && item->val_ <= max_val) {
                    uniq  += item->val_ == 1;
                    total += item->val_;
                    maxcount    = std::max(maxcount, item->val_);
                    ++distinct;
                    if (args.text_flag) {
                        sprintf(buffer, "%s %ld\n", item->key_.to_str().c_str(), item->val_);
                        buflen = strlen(buffer);
                    } else {
                        memcpy(buffer, (const char*)(item->key_).data(), key_len_);
                        buflen += key_len_;
                        uint64_t v = std::min(max_val_, item->val_);
                        for (int ii = 0; ii < val_len_; ii++) {
                            buffer[buflen + ii] = (unsigned char)(v >> (ii * 8));
                        }
                        buflen += val_len_;
                    }
                    if (buflen > MAX_ITEM_SIZE)
                        report_error("The item size is too long!\n");
                    record.set_record(buffer, buflen);
                    buflen = 0;
                    heap->pop();
                    if(it->next()) heap->push(*it);
                    return &record;
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

        return NULL;
    }

    void write(BaseRecordFormat *record) {
        KVRecord *kv = (KVRecord*)record;

        mer_dna mer(mer_dna::k());
        memcpy(mer.data__(), kv->get_key(), kv->get_key_size());
        uint64_t count = 1;
#ifdef COMBINE
        if (kv->get_val_size() != 0) {
            count = decode_varint(kv->get_val());
        }
#endif
        if((*filter)(mer)){
            switch(op){
                case COUNT: {
                    //printf("rank=%d, %s, %ld, %ld\n", rank, mer.to_str().c_str(), count, nwrite);
                    ary->add(mer, count); break;
                 }
                case PRIME: ary->set(mer); break;
                case UPDATE: {
                    ary->update_add(mer, count);
                }; break;
            }
        }
    }

    void set_min_val(uint64_t min) {
        this->min_val = min;
    }

    void set_max_val(uint64_t max) {
        this->max_val = max;
    }

  private:
    mer_hash* ary;
    struct filter* filter;
    OPERATION op;

    iterator *it;
    heap_type*  heap;
    typename storage_t::key_type key;
    size_t id;

    std::pair<size_t, size_t> block_info;
    uint64_t min_val, max_val;
    int key_len_, val_len_;
    uint64_t max_val_;

    char  buffer[MAX_ITEM_SIZE];
    int   buflen;
    BaseRecordFormat record;

    //void write_header(storage_t* ary) {
    //  header->update_from_ary(*ary);
    //  if (args.text_flag) {
    //      header->format("text/sorted");
    //  } else {
    //      header->format("binary/sorted");
    //      header->counter_len(val_len_);
    //  }
    //  std::stringstream out;
      //out.open(filename.c_str());
    //  header->write(out);
    //  std::string str = out.str();
    //  if (str.size() > MAX_ITEM_SIZE) report_error("header too large!\n");
    //  char buf[MAX_ITEM_SIZE];
    //  int buflen = (int)(str.size());
    //  sprintf(buf, "%s", str.c_str());
    //  BaseRecordFormat record(buf, buflen);
    //  writer->write(&record);
      //out.close();
    //}
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

int repartition(const char *buffer, int len, bool islast) {
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

void add_mers(const char *seq, Writable *output, void *ptr)
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
#ifndef COMBINE
                    KVRecord output_record((char*)mer.data(), param.key_len, NULL, 0);
#else
                    char val[100];
                    int vallen = encode_varint(val, 1);
                    KVRecord output_record((char*)mer.data(), param.key_len, val, vallen);
#endif
                    //printf("write a mer %s\n", mer.to_str().c_str());
                    output->write(&output_record);
                }
            }
        } else {
            param.filled = 0;
        }
    }

}

void parse_sequence(Readable *input, Writable *output, void *ptr)
{
    StringRecord* record = NULL;
    char *seq = NULL;
    char *line = NULL;

    while ((record = (StringRecord*)(input->read())) != NULL ) {
        line = record->get_record();

        if (line[0] == '>'){
            param.filled = 0;
            param.ftype = FASTA_TYPE;
            continue;
        } else if (line[0] == '@') {
            param.ftype = FASTQ_TYPE;
            record = (StringRecord*)(input->read());
            if (!record) return;
            // sequence line
            line = record->get_record();
            seq = line;
            add_mers(seq, output, ptr);
            // '+' line
            record = (StringRecord*)(input->read());
            if (!record) report_error("Fastq file format incorrect!!!");
            // score line
            record = (StringRecord*)(input->read());
            if (!record) report_error("Fastq file format incorrect!!!");
            param.filled = 0;
        } else if (param.ftype == FASTA_TYPE) {
            seq = line;
            add_mers(seq, output, ptr);
        } else {
            report_error("File format incorrect!!!");
        }

    }
}

void add_qual_mers(std::string &seq, std::string &qual, Writable *output, void *ptr)
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
#ifndef COMBINE
                    KVRecord output_record((char*)mer.data(), param.key_len, NULL, 0);
#else
                    char val[100];
                    int vallen = encode_varint(val, 1);
                    KVRecord output_record((char*)mer.data(), param.key_len, val, vallen);
#endif
                    output->write(&output_record);
                }
            }
        }else{
            param.filled = 0;
        }
    }

    seq.clear();
    qual.clear();
}

void parse_qual_sequence(Readable *input, Writable *output, void *ptr)
{
    StringRecord* record = NULL;
    const char *line = NULL;
    std::string seq, qual;

    while ((record = (StringRecord*)(input->read())) != NULL) {
        line = record->get_record();

        if (line[0] == '@') {
            record = (StringRecord*)(input->read());
            if (!record) return;
            // sequence line
            line = record->get_record();
            seq += line;
            // '+' line
            record = (StringRecord*)(input->read());
            if (!record) report_error("Fastq file format incorrect!!!");
            // score line
            record = (StringRecord*)(input->read());
            if (!record) report_error("Fastq file format incorrect!!!");
            line = record->get_record();
            qual += line;
            param.filled = 0;
        } else {
            report_error("Fastq file format incorrect!!!");
        }

        add_qual_mers(seq, qual, output, ptr);
    }
}

void combine_mers(Combinable *combiner, KVRecord *kv1, KVRecord *kv2, void* ptr)
{
    uint64_t count = decode_varint(kv1->get_val()) 
        + decode_varint(kv2->get_val());
    char val[100] = {0};
    int vallen = encode_varint(val, count);
    KVRecord update_record(kv1->get_key(),
                           kv1->get_key_size(),
                           val, vallen);
    combiner->update(&update_record);
}

int mcount_main(int argc, char *argv[])
{
  auto start_time = system_clock::now();

  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  Mimir_Init(&argc, &argv, MPI_COMM_WORLD);
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

  MerDatabase<mer_array> merdata(do_op, mer_filter.get());

  int key_len = merdata.key_len();
  param.filter = mer_filter.get();
  param.key_len = key_len / 8 + (key_len % 8 != 0);
  param.val_len = sizeof(uint32_t);
  param.filled  = 0;

#ifdef COMBINE
  int vsize = KVVARINT;
  CombineCallback combine_fn = combine_mers;
#else
  int vsize = 0;
  CombineCallback combine_fn = NULL;
#endif

  MapCallback map_fn = NULL;
  if (args.min_qual_char_given)
      map_fn = parse_qual_sequence;
  else
      map_fn = parse_sequence;

  MimirContext<char*,char*> ctx(param.key_len, vsize,
                                map_fn, NULL, input_dirs,
                                std::string(args.output_arg),
                                repartition, combine_fn);
  ctx.set_user_database(&merdata);
  ctx.map();
  ctx.output();

  auto after_count_time = system_clock::now();

  Mimir_Finalize();

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
