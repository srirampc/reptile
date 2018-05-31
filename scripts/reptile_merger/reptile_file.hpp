/*
 *  Author: Matthew Regennitter <mregenni@iastate.edu>
 */

#ifndef REPTILE_FILE_HPP
#define REPTILE_FILE_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <vector>
#include <stdint.h>

namespace mlr {

   // A structure representing a fix
   typedef struct fix_s {
      int32_t pos;
      int8_t change;
      int8_t last;
   } fix_t;

   class reptile_file 
   {
   public:

      // Constructor
      reptile_file(const std::string& file_path) {
         in_.open(file_path.c_str());
      }

      // Destructor
      ~reptile_file() {
         close();
      }

      // Method to get corrections for a particular read
      std::string get_corrections(std::vector<fix_t>& fixes) {
         std::string rc = "";
         fix_t fix;
         int32_t i = 0;

         // clear the vector
         fixes.clear();

         if(good()) {
            // Get the line
            in_.getline(buffer_, sizeof(buffer_));
            assert(in_.gcount()+1 < (int64_t)sizeof(buffer_));

            // Get the read number
            rc = strtok(buffer_, "\t");

            // Get the correction count
            const int64_t count = (int64_t)atof(strtok(NULL, "\t"));

            // Get the specified corrections
            while(i++ < count) {
               fix.pos = (int32_t)atof(strtok(NULL, "\t"));
               fix.change = (int8_t)atof(strtok(NULL, "\t"));
               fix.last = (int8_t)atof(strtok(NULL, "\t"));
               strtok(NULL, "\t");
               fixes.push_back(fix);
            }

            in_.peek();
         }
    
         return rc;
      }

      // Method to close the file (without waiting for destructor)
      void close() {
         if(in_.is_open()) {
            in_.close();
         }
      }

		// Method to determine whether the file is ready for reading
		bool good() const {
			return in_.is_open() && in_.good();
		}

   protected:

      // The input file
      std::ifstream in_;

   private:

      // Unused constructor
      reptile_file();

      // The data buffer
      char buffer_[4096];
   };

}

#endif // REPTILE_FILE_HPP

