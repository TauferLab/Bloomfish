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
#include <mapreduce.h>

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
};

// add a mer to hash array
void add_mer(mer_dna &m, uint64_t v, mer_hash& ary, filter& filter, OPERATION op){
    if(filter(m)){
        switch(op){
        case COUNT: {
            std::string mer_str = m.to_str();
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

// split sequence files
bool split_sequence_files(MIMIR_NS::InputStream *in, void *ptr){
    unsigned char ch = in->get_byte();
    if(ch == '>') return true;

    return false;
}

// read sequence files
void read_sequence_files(MIMIR_NS::MapReduce* mr,
                    MIMIR_NS::InputStream* in,
                    void* ptr)
{
    MimirParam *param = (MimirParam*)ptr;

    unsigned int filled = 0;
    mer_dna m, rcm;

    char ch=0x00;
    while( in->is_empty() == false ){
        ch = in->get_byte();
        // go to another file
        if(in->is_eof()){
            filled = 0;
        }else if(ch == '>'){
            filled = 0;
            while(in->is_eof() == false && 
                  in->get_byte() != '\n')
                in->next();
            continue;
        }else if(ch != '\n'){
            int code = m.code(ch);
            m.shift_left(code);
            if(args.canonical_flag)
                rcm.shift_right(rcm.complement(code));
            filled = std::min(filled + 1, mer_dna::k());
            if(filled >= m.k()){
                mer_dna mer = args.canonical_flag || m < rcm ? m : rcm;
                add_mer(mer, 1, *(param->ary), *(param->filter), param->op);
            }
        }
        in->next();
    }

    (*(param->ary)).done();
}

// shuffle mers
void shuffle_mers(MIMIR_NS::MapReduce *mr, void *ptr){
    MimirParam *param = (MimirParam*)ptr;
    mer_array *ary = param->ary->ary();

    std::pair<size_t, size_t> block_info; // { nb blocks, nb records }
    block_info = ary->blocks_for_records(5 * ary->max_reprobe_offset());

    for(size_t id = 0; id * block_info.second < param->ary->size(); id ++) {
        mer_dna key;
        uint64_t val;
        mer_array::eager_iterator it(ary, 
                                     id * block_info.second, 
                                    (id + 1) * block_info.second);

        while(it.next()){
            key = it.key();
            val = it.val();
            mr->add_key_value((const char*)key.data(), param->key_len,
                              (const char*)&(val), param->val_len);
        }
    }
}

// add mers from Mimir
void add_kv_mers(MIMIR_NS::MapReduce* mr, 
              char *key, int keysize, 
              char *val, int valsize, 
              void *ptr){
    MimirParam *param = (MimirParam*)ptr;
    mer_dna mer(mer_dna::k());
    memcpy(mer.data__(), key, keysize);
    uint64_t count = *(uint64_t*)val;
    add_mer(mer, count, *(param->ary), *(param->filter), param->op);
}

int mcount_main(int argc, char *argv[])
{
  auto start_time = system_clock::now();

  int provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

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

  header.canonical(args.canonical_flag);
  mer_hash ary(args.size_arg, args.mer_len_arg * 2, args.counter_len_arg, args.threads_arg, args.reprobes_arg);
  ary.do_size_doubling(true);
  if(args.disk_flag){
    //ary.do_size_doubling(false);
    fprintf(stderr, "mcount currently does not support --disk parameter.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  std::unique_ptr<jellyfish::dumper_t<mer_array> > dumper;
  if(args.text_flag)
    dumper.reset(new text_dumper(args.threads_arg, args.output_arg, &header));
  else
      dumper.reset(new binary_dumper(args.out_counter_len_arg, ary.key_len(), args.threads_arg, args.output_arg, &header));
  ary.dumper(dumper.get());

  MimirParam param;
  std::unique_ptr<MIMIR_NS::MapReduce> mimir;
  mimir.reset(new MIMIR_NS::MapReduce(MPI_COMM_WORLD));

  int key_len = ary.key_len();
  param.key_len = key_len / 8 + (key_len % 8 != 0);
  param.val_len = sizeof(uint64_t);
  mimir->set_key_length(param.key_len);
  mimir->set_value_length(sizeof(uint64_t));

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

  if(args.min_qual_char_given) {
    fprintf(stderr, "mcount currently does not support -Q, --min-qual-char=string parameter.\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  } else {
    param.ary = &ary;
    param.filter = mer_filter.get();
    param.op = do_op;
    std::string path = *(args.file_arg.begin());
    // local counting of mers
    mimir->process_binary_file(path.c_str(), 1, 1, read_sequence_files, 
                               split_sequence_files, &param, 0);
    // Shuffle
    mimir->init_key_value(shuffle_mers, &param);
    ary.ary()->clear();
    // read back mers from Mimir
    mimir->map_key_value(mimir.get(), add_kv_mers, &param, 0);
  }

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
      dumper->dump(ary.ary());
      if(!args.no_merge_flag) {
        std::vector<const char*> files = dumper->file_names_cstr();
        uint64_t min = args.lower_count_given ? args.lower_count_arg : 0;
        uint64_t max = args.upper_count_given ? args.upper_count_arg : std::numeric_limits<uint64_t>::max();
        try {
          merge_files(files, args.output_arg, header, min, max);
        } catch(MergeError e) {
          err::die(err::msg() << e.what());
        }
        if(!args.no_unlink_flag) {
          for(int i =0; i < dumper->nb_files(); ++i)
            unlink(files[i]);
        }
      } // if(!args.no_merge_flag
    } // if(!args.no_merge_flag
  }

  mimir.release();
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
