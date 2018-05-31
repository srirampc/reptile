/*
 *  Author: Matthew Regennitter <mregenni@iastate.edu>
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdint.h>

#include <cassert>

#include <zlib.h>
#include "kseq.h"
#include "reptile_file.hpp"

KSEQ_INIT(gzFile, gzread)

char LOOKUP[4] = { 'A', 'C', 'G', 'T' };

int32_t main(int32_t argc, char* argv[]) {
   
   if(argc < 4) {
      std::cerr << "\nTool to produce an output file from a fasta file and its Reptile corrections" << std::endl;
      std::cerr << "\nUsage: " << argv[0] << " <source> <corrections> <output>" << std::endl << std::endl;
	  std::cerr << " <source> : the raw fasta file. The read names are 0 to (readNum - 1)" << std::endl;
	  std::cerr << " <corrections> : reptile correction of the raw fasta file" << std::endl;
	  std::cerr << " <output> : the user specified output fasta file name" << std::endl << std::endl;		
      return 1;
   }

   std::ifstream source(argv[1]);
	if(!source.good()) {
		std::cerr << "Could not open the source file " << argv[1] << ". Does the path exist?" << std::endl;
		return 1;
	}
    source.close();
	
   mlr::reptile_file corrections(argv[2]);
	if(!corrections.good()) {
		std::cerr << "Could not open the corrections file " << argv[2] << ". Does the path exist?" << std::endl;
		return 1;
	}

   std::ofstream out(argv[3]);

   if(!out.good()) {
      std::cerr << "Could not open output file " << argv[3] << ". Does the path exist?" << std::endl;
		out.close();
      return 1;
   }

    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "Input fasta/fastq file  : " << argv[1] << std::endl;
    std::cout << "Reptile corrections     : " << argv[2] << std::endl;
    std::cout << "Output fasta/fastq file : " << argv[3] << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;

   	gzFile fp;
	kseq_t *seq;
	int l;
	fp = gzopen(argv[1], "r");
	seq = kseq_init(fp);
   
   std::vector<mlr::fix_t> fixes;

   // Get the initial corrections
   std::string cindex = corrections.get_corrections(fixes);

   while ((l = kseq_read(seq)) >= 0) {
      std::string name(seq->name.s);
      std::string read(seq->seq.s);
      std::string comment;
      std::string qual;
      if (seq->comment.l) comment = seq->comment.s;
      if (seq->qual.l) qual = seq->qual.s;


      // Check the ID of this read
      std::string rindex = name;
      //std::cout << read << rindex << std::endl;

      if(cindex == rindex) {
         // Make the corrections
         for(uint32_t i = 0; i < fixes.size(); ++i) {
			   assert(read[fixes[i].pos] == LOOKUP[fixes[i].last]);
            read[fixes[i].pos] = LOOKUP[fixes[i].change];
         }
         // Get the next set of corrections
         cindex = corrections.get_corrections(fixes);         
      }
      //std::cout << read << std::endl;

      // Commmit the read (possibly with changes)
      if(seq->qual.l){
        out << "@" << name << std::endl;
        out << read << std::endl;
        out << "+" << name << std::endl;
        out << qual << std::endl;
      } else {
        out << ">" << name << std::endl;
        out << read << std::endl;
      }
   }

   out.close();
    kseq_destroy(seq);
    gzclose(fp);

   return 0;
}
