/***
 **
 *  File: Kmer.cpp
 *  Created: December 6, 2009
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

#include "Kmer.h"
#include "Seq.h"

Kmer::Kmer() {
}

Kmer::Kmer(const Kmer& orig) {
}


bool Kmer::goodQuality(char* qAddr, int kvalue, const Para& myPara){
    if (!qAddr) return false;
    int badpos = 0;
    for (int i = 0; i < kvalue; ++ i){
        if ((int) qAddr[i] < myPara.qualThreshold){
            badpos ++;
        }
    }
    if (badpos > myPara.maxBadQPerKmer) return false;
    return true;
}

/*
 * Get kmer information, resulting sorted <ID, good_cnt, cnt> w.r.t ID field
 */
void Kmer::getKC(const Seq& mySeq, int kvalue, const Para& myPara) {

    std::vector<std::pair<uint64_t, char*> > KQ_array; // <ID, qualaddr>

    int numSeq = mySeq.getSL().size();

    for (int i = 0; i < numSeq; ++i) {

        int pos = mySeq.getSL()[i];

        char* addr = const_cast<char*> (&mySeq.getS()[pos]);
        char* qAddr = NULL;
        if (myPara.QFlag){
            qAddr = const_cast<char*> (&mySeq.getQ()[pos]);
        }

        int len = (int) strlen(addr) - kvalue + 1;

        for (int j = 0; j < len; ++j) {
            uint64_t ID;
            if (toID (ID, addr + j, kvalue)){
                if (qAddr){
                    KQ_array.push_back(std::pair<uint64_t, char*> (ID, qAddr + j));
                }
                else {
                    KQ_array.push_back(std::pair<uint64_t, char*> (ID, NULL));
                }
            }
        }
    }

    //std::cout << "\tsorting ...\n";
    std::sort(KQ_array.begin(), KQ_array.end(), KQcomp());

    // linear scan to remove the duplicates and calculating multiplicity
    int size = (int) KQ_array.size();

    int idx1 = 0;

    int card = 0;
    while (idx1 < size - 1) {

        card = 1;
        int goodCnt = 0;

        if (goodQuality(KQ_array[idx1].second, kvalue, myPara)) goodCnt++;

        int idx2 = idx1 + 1;

        while (idx2 < size) {

            if (KQ_array[idx1].first == KQ_array[idx2].first) {
                //two kmers non-overlapping on the same read
                //if(Kmer_SP[idx2] - Kmer_SP[idx1] >= K)
                if (goodQuality(KQ_array[idx2].second, kvalue, myPara)) goodCnt++;
                ++card;
                ++idx2;
            } else break;
        }

        elem_.push_back(kc_t(KQ_array[idx1].first, goodCnt, card));

        idx1 = idx2;
    }

    if (idx1 == size - 1) {
        int goodCnt = 0;

        if (goodQuality(KQ_array[idx1].second, kvalue, myPara)) goodCnt++;

        elem_.push_back(kc_t(KQ_array[idx1].first, goodCnt, card));
    }
}