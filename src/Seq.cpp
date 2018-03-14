
/***
 **
 *  File: Seq.cpp
 *  Created: Dec 05, 2009 8:27 PM
 *
 *  Author: Xiao Yang <isuxyang@gmail.com>
 * 
 *  Copyright 2009 Xiao Yang, Karin Dorman, Srinivas Aluru
 *
 *  This file is part of Reptile (version 1.1)
 *
 *  Reptile is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Reptile is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Libpnorm. If not, see <http://www.gnu.org/licenses/>.
 *
 */
 
#include "Seq.h"

Seq::Seq() {
}

Seq::Seq(const Seq& orig) {
}

Seq::~Seq() {
}

// flag = 0: short reads file; 1: quality file
void Seq::retrieve_seqs_from_fa(bIO::FASTA_input& handle, int max, bool flag){     

    typedef bIO::FASTA_input::value_type value_type;

    int position = 0;
    int num = 0;
    
    while (++ handle == true) {

            const value_type& v = *handle;

            if (!flag){

                headers_.push_back(v.first);

                SL_.push_back(position);

                S_.resize(S_.size() + v.second.length() + 1);

                memcpy(&S_[position], v.second.c_str(), v.second.length() + 1);

                position += v.second.length() + 1;
            }
            else {
                QL_.push_back(position);

                std::string quals = v.second;

                quals.erase(0, 1);

                Q_.resize(Q_.size() + quals.length() + 1);

                memcpy(&Q_[position], quals.c_str(), quals.length() + 1);

                position += quals.length() + 1;
            }       

            num ++;
            
            if (num >= max) break;
    }
    //converting everything to lower case
    if (!flag){
        for (int i = 0; i < (int) S_.size(); ++ i){
            if(isupper(S_[i])) {
                S_[i] = tolower(S_[i]);
            }
        }
    }

}


 