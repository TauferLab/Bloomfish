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
#include <jellyfish/jellyfish.hpp>
#include <jellyfish/merge_files.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <sub_commands/count_main_cmdline.hpp>

#include "mimir.h"
using namespace MIMIR_NS;

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

extern mer_dna_bloom_counter* load_bloom_filter(const char* path);
enum OPERATION { COUNT, PRIME, UPDATE };

struct MimirParam{
    mer_hash      *ary;
    struct filter *filter;
    OPERATION      op;
    int  key_len, val_len;
    mer_dna m, rcm;
    std::string seq, qual;
    unsigned int filled;
};

MimirParam param;
int rank, size;

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

int fastq_padding(char *buffer, int len) {
    int padding = 0;
    while (padding < len && buffer[padding] != '\n') {
        param.seq += buffer[padding];
        padding ++;
    }
    padding ++;

    if (padding != '+') {
        fprintf(stderr, "Cannot find '+' in the padding for fastq file\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    while (padding < len && buffer[padding] != '\n')
        padding ++;
    padding ++;

    while (padding < len && buffer[padding] != '\n') {
        param.qual += buffer[padding];
        padding ++;
    }
    if (buffer[padding] == '\n') padding++;

    if (param.seq.size() != param.qual.size()) {
        fprintf(stderr, "The sequence and quality scores don't mactch!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    param.filled = 0;
    if (param.seq.size() > mer_dna::k() - 1) {
        param.seq = param.seq.substr(param.seq.size() - mer_dna::k() + 1);
        param.qual = param.qual.substr(param.qual.size() - mer_dna::k() + 1);
    }

    if (!args.min_qual_char_given) {
        for (size_t j = 0; j < param.seq.size(); j++) {
            char ch = param.seq[j];
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
        param.seq.clear();
        param.qual.clear();
    }

    return padding;
}

class MerRecord : public ByteRecord {
  public:
    virtual int get_left_border(char *buffer, uint64_t len, bool islast) {
        int count = 0;
        int padding = 0;
        int line_num = 5;
        bool has_qual = false;

        for (int i = 0; i < line_num; i++) {
            while ((uint64_t)count < len && buffer[count] != '\n') count ++;
            if ((uint64_t)count == len && !islast) {
                fprintf(stderr, "No enough line exsit in a partition\n");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            if ((uint64_t)count == len) return count;
            if ((uint64_t)count == len - 1) return count + 1;
            count ++;
            if (args.min_qual_char_given && buffer[count] == '>') {
                fprintf(stderr, "Fastq file format error with '>'!\n");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            if (buffer[count] == '>' || buffer[count] == '@')
                return count;
            if (buffer[count] == '+') {
                has_qual = true;
                line_num = i + 3;
            }
        }

        if (has_qual) padding = fastq_padding(buffer + count, (int)len - count);
        else padding = fasta_padding(buffer +count, (int)len - count);
        return count + padding;
    }
};

class MerCounter : public Writable {
  public:
    MerCounter(bool localcount = true) {
        this->localcount = localcount;
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
        uint32_t count = 1;
        if (!localcount) count = *(uint32_t*)(kv->get_val());
        add_mer(mer, count, *(param.ary), *(param.filter), param.op);
        record_count ++;
    }
    uint64_t get_record_count() { return record_count; }
  private:
    bool localcount;
    uint64_t record_count;
};

bool skip_line(Readable* input)
{
    bool ret = false;

    MerRecord* record = NULL;

    while ((record = (MerRecord*)(input->read())) != NULL) {
        char ch = *(record->get_record());
	if (ch == '\n') break;
	if (record->is_eof()) {
            ret = true;
            break;
        }
    }

    return ret;
}

void parse_sequence(Readable *input, Writable *output, void *ptr)
{
    MerRecord* record = NULL;

    char ch;

    while ((record = (MerRecord*)(input->read())) != NULL) {
        ch = *(record->get_record());
        if (ch == '>' || ch == '@') {
            param.filled = 0;
            skip_line(input);
            continue;
        }
        else if (ch == '+') {
            skip_line(input);
            skip_line(input);
            continue;
        }
        else if (ch == '\n') continue;
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
                    KVRecord output_record((char*)mer.data(), param.key_len, NULL, 0);
                    output->write(&output_record);
                }
            }
        } else {
            param.filled = 0;
        }

        if (record->is_eof()) param.filled = 0;
    }
}

void add_qual_mers(Writable *output, void *ptr)
{
    if (param.seq.size() != param.qual.size()) {
        fprintf(stderr, "the length of quality score (%ld) does not match sequence length (%ld)\n", 
                param.qual.size(), param.seq.size());
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::string::iterator seq_iter = param.seq.begin();
    std::string::iterator qual_iter = param.qual.begin();
    for (;seq_iter != param.seq.end(); seq_iter++, qual_iter++) {
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
                    KVRecord output_record((char*)mer.data(), param.key_len, NULL, 0);
                    output->write(&output_record);
                }
            }
        }else{
            param.filled = 0;
        }
    }

    param.seq.clear();
    param.qual.clear();
}

void parse_fastq_sequence(Readable *input, Writable *output, void *ptr)
{
    char ch;
    ByteRecord* record = NULL;

    while ((record = (ByteRecord*)(input->read())) != NULL) {
        ch = *(record->get_record());
        if (ch == '@') {
            add_qual_mers(output, ptr);
            param.filled = 0;
            skip_line(input);
            continue;
        }

        param.seq += ch;
        while ((record = (ByteRecord*)(input->read())) != NULL) {
            ch = *(record->get_record());
            if (ch == '\n') {
                break;
            }
            param.seq += ch;
        }

        record = (ByteRecord*)(input->read());
        ch = *(record->get_record());

        if (ch != '+') {
            fprintf(stderr, "The fastq file format error!\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        } else {
            skip_line(input);
        }

        while ((record = (ByteRecord*)(input->read())) != NULL) {
            ch = *(record->get_record());
            if (ch == '\n') {
                break;
            }
            param.qual += ch;
        }

        add_qual_mers(output, ptr);

        if (record->is_eof()) param.filled = 0;
    }
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
  char outfile[1024]={0};
  sprintf(outfile, "%s.%d.%d", args.output_arg, size, rank);
  std::unique_ptr<jellyfish::dumper_t<mer_array> > dumper;
  if(args.text_flag)
    dumper.reset(new text_dumper(1, outfile, &header));
  else
      dumper.reset(new binary_dumper(args.out_counter_len_arg, ary.key_len(), 
                                     1, outfile, &header));
  ary.dumper(dumper.get());

  MimirContext mimir;

  int key_len = ary.key_len();
  param.key_len = key_len / 8 + (key_len % 8 != 0);
  param.val_len = sizeof(uint32_t);
  param.filled  = 0;
  mimir.set_key_length(param.key_len);
  mimir.set_val_length(0);
  //mimir->set_combiner(combiner);

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
  std::string path = *(args.file_arg.begin());

  if(rank == 0) printf("start reading files\n");
  // local counting of mers
  InputSplit* splitinput = FileSplitter::getFileSplitter()->split(path.c_str(), BYSIZE);
  splitinput->print();
  //ByteRecord::set_separators(">@");
  FileReader<MerRecord> reader(splitinput);
  //read_sequence_files((FileReader<ByteRecordFormat>*)&reader, &param);
  MerCounter counter;
  if (args.min_qual_char_given)
      mimir.set_map_callback(parse_fastq_sequence);
  else
      mimir.set_map_callback(parse_sequence);
  //mimit.set_shuffle_flag(false);
  mimir.mapreduce(&reader, &counter, NULL);

  auto after_count_time = system_clock::now();

  // If no intermediate files, dump directly into output file. If not, will do a round of merging
  if(!args.no_write_flag) {
    if(dumper->nb_files() == 0) {
      dumper->one_file(true);
      if(args.lower_count_given)
        dumper->min(args.lower_count_arg);
      if(args.upper_count_given)
        dumper->max(args.upper_count_arg);
      dumper->dump(ary.ary());
    } else {
      fprintf(stderr, "mcount currently does not support intermediate files.\n");
      MPI_Abort(MPI_COMM_WORLD, 1);
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
    } // if(!args.no_merge_flag
  }

  Mimir_Finalize();
  MPI_Finalize();

  auto after_dump_time = system_clock::now();

  if(args.timing_given) {
    sprintf(outfile, "%s.%d.%d", args.timing_arg, size, rank);
    std::ofstream timing_file(outfile);
    timing_file << "Init     " << as_seconds(after_init_time - start_time) << "\n"
                << "Counting " << as_seconds(after_count_time - after_init_time) << "\n"
                << "Writing  " << as_seconds(after_dump_time - after_count_time) << "\n";
    timing_file << "Array " << ary.size() << "\n";
  }

  return 0;
}
