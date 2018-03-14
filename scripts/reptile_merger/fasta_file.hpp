/*
 *  Author: Matthew Regennitter <mregenni@iastate.edu>
 */
#ifndef FASTA_FILE_HPP
#define FASTA_FILE_HPP

#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

namespace mlr {

   class fasta_file {

   public:

      // Constructor
      fasta_file(const std::string& file_path) {
         in_.open(file_path.c_str());
      }

      // Destructor
      ~fasta_file() {
         close();
      }

      // Method to retrieve a read and name from the file
      bool get_read(std::string& name, std::string& read) {
         bool rc = false;
         std::string buff;
         read.clear();
         name.clear();
         
         if(good()) {
            getline(in_, buff);
            if('>' == buff[0]) {
               name = &(buff[1]);
               while(in_.good() && '>' != in_.peek()) {
                  getline(in_, buff);
                  read += buff;
                  in_.peek();
               }
               rc = true;
            }
            in_.peek();
         }
         
         return rc;
      }

      // Method to retrieve a read from the file
      bool get_read(std::string& read) {
         std::string dummy;
         return get_read(dummy, read);
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
      fasta_file();

   };

}

#endif // FASTA_FILE_HPP

