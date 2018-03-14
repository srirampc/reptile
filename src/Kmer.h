/***
 **
 *  File: Kmer.h
 *  Created: Dec 6, 2009
 *
 *  Author: Xiao Yang <isuxyang@gmail.com>
 *
 *  Copyright 2009 Xiao Yang, Karin Dorman, Srinivas Aluru
 *
 *  This file is part of Reptile (version 1.1).
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

#ifndef _KMER_H
#define	_KMER_H

#include "util.h"
#include "Seq.h"


// used by sorting KQ_array, compare length K cstr
struct KQcomp {
    bool operator() (const std::pair<uint64_t, char*>& e1,
                     const std::pair<uint64_t, char*>& e2) const {
        return e1.first < e2.first;
    }
};

typedef struct KC{
    uint64_t ID;
    int goodCnt;
    int cnt;
    KC(){};
    KC(uint64_t myID, int c1, int c2): ID(myID), goodCnt(c1), cnt(c2){};
}kc_t;

struct Knumcomp {
    bool operator() (const kc_t& e1, const kc_t& e2) const {
        return  (e1.ID < e2.ID);
    }
};

typedef std::vector<kc_t> kcvec_t;

class Kmer {

protected:
    kcvec_t elem_;
    bool goodQuality(char* qAddr,int kvalue, const Para& myPara);
public:
    Kmer();
    Kmer(const Kmer& orig);
    void getKC(const Seq& mySeq, int kvalue, const Para& myPara);
    const kcvec_t& getElem () const { return elem_;  }
    void clear(){
        elem_.clear();
    }
};




#endif	/* _KMER_H */

