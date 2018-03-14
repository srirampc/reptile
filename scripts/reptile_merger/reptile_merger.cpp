/*
 *  Author: Matthew Regennitter <mregenni@iastate.edu>
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdint.h>

#include <cassert>

#include "fasta_file.hpp"
#include "reptile_file.hpp"

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

   mlr::fasta_file source(argv[1]);
	if(!source.good()) {
		std::cerr << "Could not open the source file " << argv[1] << ". Does the path exist?" << std::endl;
		return 1;
	}
	
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
   
   std::string read, name;
   std::vector<mlr::fix_t> fixes;

   // Get the initial corrections
   int64_t cindex = corrections.get_corrections(fixes);

   while(source.get_read(name, read)) {

      // Check the ID of this read
      int64_t rindex = (int64_t)atof(name.c_str());

      if(cindex == rindex) {
         // Make the corrections
         for(uint32_t i = 0; i < fixes.size(); ++i) {
			   assert(read[fixes[i].pos] == LOOKUP[fixes[i].last]);
            read[fixes[i].pos] = LOOKUP[fixes[i].change];
         }
         // Get the next set of corrections
         cindex = corrections.get_corrections(fixes);         
      }

      // Commmit the read (possibly with changes)
      out << ">" << name << std::endl;
      out << read << std::endl;
   }

   out.close();

   return 0;
}
