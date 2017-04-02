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

//#define COMBINE

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

int rank, size;

void report_error(const char *msg) {
    fprintf(stderr, "%d[%d] Error: %s\n", rank, size, msg);
    MPI_Abort(MPI_COMM_WORLD, 1);
}

template<typename storage_t>
class mimir_dumper
{
  public:
    typedef typename storage_t::key_type key_type;
    typedef typename storage_t::region_iterator iterator;
    typedef typename jellyfish::mer_heap::heap<typename storage_t::key_type, iterator> heap_type;
    typedef typename heap_type::const_item_t    heap_item;

    mimir_dumper(const char* filename, jellyfish::file_header* header)
    {
        this->header = header;
        writer = FileWriter::getWriter(filename);
        this->filename = (writer->get_file_name()).c_str();
        min = 0;
        max = std::numeric_limits<uint64_t>::max();
        val_len_ = args.out_counter_len_arg;
        max_val_ = (((uint64_t)1 << (8 * val_len_)) - 1);
    }

    ~mimir_dumper() {
        //delete writer;
    }

    void write_header(storage_t* ary) {
      header->update_from_ary(*ary);
      if (args.text_flag) {
          header->format("text/sorted");
      } else {
          header->format("binary/sorted");
          header->counter_len(val_len_);
      }
      std::stringstream out;
      //out.open(filename.c_str());
      header->write(out);
      std::string str = out.str();
      if (str.size() > MAX_ITEM_SIZE) report_error("header too large!\n");
      char buf[MAX_ITEM_SIZE];
      int buflen = (int)(str.size());
      sprintf(buf, "%s", str.c_str());
      BaseRecordFormat record(buf, buflen);
      writer->write(&record);
      //out.close();

    }

    void dump(storage_t* ary) {
        char  buffer[MAX_ITEM_SIZE];
        int   buflen = 0;
        storage_t* ary_ = ary;
        std::pair<size_t, size_t> block_info 
            = ary_->blocks_for_records(5 * ary_->max_reprobe_offset());
        heap_type             heap(ary_->max_reprobe_offset());
        size_t                count = 0;
        typename storage_t::key_type key;

        key_len_ = ary->key_len() / 8 + (ary->key_len() % 8 != 0);

        //MPI_Barrier(MPI_COMM_WORLD);
        writer->open();

        if (writer->is_single_file()) {
            if (rank == 0) write_header(ary);
        } else {
            write_header(ary);
        }

        for(size_t id = 0; id * block_info.second < ary_->size(); id += 1) {
            // Fill buffer
            iterator it(ary_, id * block_info.second, (id + 1) * block_info.second, key);
            heap.fill(it);

            while(heap.is_not_empty()) {
                 heap_item item = heap.head();
                if(item->val_ >= min && item->val_ <= max) {
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
                    BaseRecordFormat record(buffer, buflen);
                    writer->write(&record);
                    buflen = 0;
                }
                ++count;
                heap.pop();
                if(it.next())
                    heap.push(it);
            }
        }

        writer->close();
     }

    void set_min(uint64_t min) {
        this->min = min;
    }

    void set_max(uint64_t max) {
        this->max = max;
    }

  private:
    std::string filename;
    FileWriter *writer;
    jellyfish::file_header *header;
    uint64_t                 min;
    uint64_t                 max;
    int val_len_;
    uint64_t max_val_;
    int key_len_;
};

extern mer_dna_bloom_counter* load_bloom_filter(const char* path);
enum OPERATION { COUNT, PRIME, UPDATE };

enum FileType { NULL_TYPE, FASTA_TYPE, FASTQ_TYPE };

struct MimirParam{
    mer_hash      *ary;
    struct filter *filter;
    OPERATION      op;
    int  key_len, val_len;
    mer_dna m, rcm;
    unsigned int filled;
    FileType ftype;
};

MimirParam param;

// add a mer to hash array
void add_mer(mer_dna &m, uint64_t v,
             mer_hash& ary, filter& filter, OPERATION op)
{
    if(filter(m)){
        switch(op){
        case COUNT: {
            //printf("rank=%d, %s\n", rank, m.to_str().c_str());
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

int fasta_padding(char *buffer, int len) {
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

class MerRecord : public StringRecord {
  public:
    virtual int get_border_size(char *buffer, uint64_t len, bool islast) {
        int padding = 0;
        int pre_pos = 0;
        char pre = 0;
        int pos = 0;
        int line_num = 5;

        for (int i=0; i < line_num; i++) {
            while ((uint64_t)pos < len && buffer[pos] != '\n') pos ++;
            if ((uint64_t)pos == len && !islast) {
                report_error("Cannot judge the border.");
            }
            if ((uint64_t)pos == len || (uint64_t)pos == (len - 1))
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
};

class MerCounter : public Writable {
  public:
    MerCounter() {
    }
    bool open() {
        record_count = 0;
        return true;
    }
    void close() {
    }
    void write(BaseRecordFormat *record) {
        KVRecord *kv = (KVRecord*)record;
        mer_dna mer(mer_dna::k());
        memcpy(mer.data__(), kv->get_key(), kv->get_key_size());
        uint64_t count = 1;
#ifdef COMBINE
        if (kv->get_val_size() != 0) {
            //count = strtoull((kv->get_val()), NULL, 0);
            count = decode_varint(kv->get_val());
        }
#endif
        add_mer(mer, count, *(param.ary), *(param.filter), param.op);
        record_count ++;
    }
    uint64_t get_record_count() { return record_count; }
  private:
    uint64_t record_count;
};

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
                    //char tmp[2] = {"1"};
                    char val[100];
                    int vallen = encode_varint(val, 1);
                    KVRecord output_record((char*)mer.data(), param.key_len, val, vallen);
#endif
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
    MerRecord* record = NULL;
    char *seq = NULL;
    char *line = NULL;

    while ((record = (MerRecord*)(input->read())) != NULL ) {
        line = record->get_record();

        if (line[0] == '>'){
            param.filled = 0;
            param.ftype = FASTA_TYPE;
            continue;
        } else if (line[0] == '@') {
            param.ftype = FASTQ_TYPE;
            record = (MerRecord*)(input->read());
            if (!record) return;
            // sequence line
            line = record->get_record();
            seq = line;
            add_mers(seq, output, ptr);
            // '+' line
            record = (MerRecord*)(input->read());
            if (!record) report_error("Fastq file format incorrect!!!");
            // score line
            record = (MerRecord*)(input->read());
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
                    //char tmp[2] = {"1"};
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
    MerRecord* record = NULL;
    const char *line = NULL;
    std::string seq, qual;

    while ((record = (MerRecord*)(input->read())) != NULL) {
        line = record->get_record();

        if (line[0] == '@') {
            record = (MerRecord*)(input->read());
            if (!record) return;
            // sequence line
            line = record->get_record();
            seq += line;
            // '+' line
            record = (MerRecord*)(input->read());
            if (!record) report_error("Fastq file format incorrect!!!");
            // score line
            record = (MerRecord*)(input->read());
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
  mer_hash ary(args.size_arg, args.mer_len_arg * 2, args.counter_len_arg, 1, args.reprobes_arg);
  ary.do_size_doubling(true);
  if(args.disk_flag){
    //ary.do_size_doubling(false);
    fprintf(stderr, "mcount currently does not support --disk parameter.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // do not support multithreading
  //char outfile[1024]={0};
  //sprintf(outfile, "%s.%d.%d", args.output_arg, size, rank);
  //std::unique_ptr<jellyfish::dumper_t<mer_array> > dumper;
  //if(args.text_flag)
  //  dumper.reset(new text_dumper(1, outfile, &header));
  //else
  //    dumper.reset(new binary_dumper(args.out_counter_len_arg, ary.key_len(), 
  //                                   1, outfile, &header));
  //ary.dumper(dumper.get());

  MimirContext mimir;

  int key_len = ary.key_len();
  param.key_len = key_len / 8 + (key_len % 8 != 0);
  param.val_len = sizeof(uint32_t);
  param.filled  = 0;
  mimir.set_key_length(param.key_len);
#ifdef COMBINE
  mimir.set_val_length(KVVARINT);
#else
  mimir.set_val_length(0);
#endif

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

  //if(args.min_qual_char_given) {
  //  fprintf(stderr, "mcount currently does not support -Q, --min-qual-char=string parameter.\n");
  //  MPI_Abort(MPI_COMM_WORLD, 1);
  //} else {
  param.ary = &ary;
  param.filter = mer_filter.get();
  param.op = do_op;
  //std::string path = *(args.file_arg.begin());

  if(rank == 0) printf("start reading files\n");
  // local counting of mers
  InputSplit input;
  for (file_vector::const_iterator iter = args.file_arg.begin(); 
       iter != args.file_arg.end(); iter++) {
      std::string file = *(iter);
      input.add(file.c_str());
  }
  InputSplit* splitinput = FileSplitter::getFileSplitter()->split(&input, BYSIZE);
  splitinput->print();
  StringRecord::set_whitespaces("\n");
  FileReader<MerRecord> *reader = FileReader<MerRecord>::getReader(splitinput);
  //read_sequence_files((FileReader<ByteRecordFormat>*)&reader, &param);
  MerCounter counter;
  if (args.min_qual_char_given)
      mimir.set_map_callback(parse_qual_sequence);
  else
      mimir.set_map_callback(parse_sequence);
  //mimit.set_shuffle_flag(false);
#ifdef COMBINE
  mimir.set_combine_callback(combine_mers);
#endif
  mimir.mapreduce(reader, &counter, NULL);

  auto after_count_time = system_clock::now();

  // If no intermediate files, dump directly into output file. If not, will do a round of merging
  if(!args.no_write_flag) {
    //if(dumper->nb_files() == 0) {
      //dumper->one_file(true);
      //MPIFileWriter writer(args.output_arg);
      mimir_dumper<mer_array> dumper(args.output_arg, &header);
      if(args.lower_count_given)
          dumper.set_min(args.lower_count_arg);
      if(args.upper_count_given)
          dumper.set_max(args.upper_count_arg);
      dumper.dump(ary.ary());
    //} else {
    //  fprintf(stderr, "mcount currently does not support intermediate files.\n");
    //  MPI_Abort(MPI_COMM_WORLD, 1);
      //dumper->dump(ary.ary());
      //if(!args.no_merge_flag) {
      //  std::vector<const char*> files = dumper->file_names_cstr();
      //  uint64_t min = args.lower_count_given ? args.lower_count_arg : 0;
      //  uint64_t max = args.upper_count_given ? args.upper_count_arg : std::numeric_limits<uint64_t>::max();
      //  try {
      //    merge_files(files, args.output_arg, header, min, max);
      //  } catch(MergeError e) {
      //    err::die(err::msg() << e.what());
      //  }
      //  if(!args.no_unlink_flag) {
      //    for(int i =0; i < dumper->nb_files(); ++i)
      //      unlink(files[i]);
      //  }
      //} // if(!args.no_merge_flag
    //} // if(!args.no_merge_flag
  }

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
    //sprintf(outfile, "%s.%d.%d", args.timing_arg, size, rank);
    std::ofstream timing_file(args.timing_arg);
    timing_file << "Init     " << as_seconds(after_init_time - start_time) << "\n"
                << "Counting " << as_seconds(after_count_time - after_init_time) << "\n"
                << "Writing  " << as_seconds(after_dump_time - after_count_time) << "\n";
    timing_file << "Array " << ary.size() << "\n";
  }

  return 0;
}
