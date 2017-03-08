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

#include <config.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <iostream>
#include <fstream>

#include <jellyfish/err.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/fstream_default.hpp>
#include <jellyfish/jellyfish.hpp>
#include <sub_commands/stats_main_cmdline.hpp>

#include "mimir.h"
using namespace MIMIR_NS;

namespace err = jellyfish::err;

template<typename reader_type>
void compute_stats(reader_type& reader, uint64_t low, uint64_t high,
                   uint64_t& uniq, uint64_t& distinct, uint64_t& total,
                   uint64_t& max) {
  //uniq = distinct = total = max = 0;

  while(reader.next()) {
    if(reader.val() < low || reader.val() > high) continue;
    uniq  += reader.val() == 1;
    total += reader.val();
    max    = std::max(max, reader.val());
    ++distinct;
  }
}

int mstats_main(int argc, char *argv[])
{
  stats_main_cmdline args(argc, argv);

  int provided, rank, size;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  Mimir_Init(&argc, &argv, MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  FileSplitter spliter;
  InputSplit*  input = spliter.split(args.db_arg, BYNAME);

  uint64_t uniq = 0, distinct = 0, total = 0, max = 0;

  FileSeg *seg = NULL;
  while ((seg = input->get_next_file()) != NULL) {
      printf("open file=%s\n", seg->filename.c_str());
      std::ifstream is(seg->filename.c_str());
      if(!is.good())
          err::die(err::msg() << "Failed to open input file '" << args.db_arg << "'");
      jellyfish::file_header header;
      header.read(is);
      jellyfish::mer_dna::k(header.key_len() / 2);
 
      if(!args.upper_count_given)
          args.upper_count_arg = std::numeric_limits<uint64_t>::max();
      if(!header.format().compare(binary_dumper::format)) {
          binary_reader reader(is, &header);
          compute_stats(reader, args.lower_count_arg, args.upper_count_arg, uniq, distinct, total, max);
      } else if(!header.format().compare(text_dumper::format)) {
          text_reader reader(is, &header);
          compute_stats(reader, args.lower_count_arg, args.upper_count_arg, uniq, distinct, total, max);
      } else {
          err::die(err::msg() << "Unknown format '" << header.format() << "'");
      } 
  }

  uint64_t global_uniq = 0, global_distinct = 0, global_total = 0, global_max = 0;
 
  MPI_Reduce(&uniq, &global_uniq, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&distinct, &global_distinct, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&total, &global_total, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&max, &global_max, 1, MPI_INT64_T, MPI_MAX, 0, MPI_COMM_WORLD);

  if (rank == 0) {

      ofstream_default out(args.output_given ? args.output_arg : 0, std::cout);
      if(!out.good())
      err::die(err::msg() << "Error opening output file '" << args.output_arg << "'");

      out << "Unique:    " << global_uniq << "\n"
          << "Distinct:  " << global_distinct << "\n"
          << "Total:     " << global_total << "\n"
          << "Max_count: " << global_max << "\n";
      out.close();
  }

  Mimir_Finalize();
  MPI_Finalize();

  return 0;
}
