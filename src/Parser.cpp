/***
 **
 *  File: Parser.cpp
 *  Created: Dec 12, 2009 4:05 PM
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
 

#include <algorithm>

#include "Parser.h"

bool file_exists(const char *fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}

/*
 * Merge all kmers in mySeq to TargetArray, which stores unique kmers
 */

void mergeK(uvec_t& TargetArray, const Seq& mySeq, int kvalue) {

    uvec_t tmpKArray;

    int numSeq = mySeq.getSL().size();

    for (int i = 0; i < numSeq; ++i) {

        int pos = mySeq.getSL()[i];

        char* addr = const_cast<char*> (&mySeq.getS()[pos]);
        int len = (int) strlen(addr) - kvalue + 1;

        for (int j = 0; j < len; ++j) {
            uint32_t ID;
            if (toID (ID, addr + j, kvalue)){
                tmpKArray.push_back(ID);
            }
        }
    }

    TargetArray.insert(TargetArray.end(), tmpKArray.begin(), tmpKArray.end());

    std::sort(TargetArray.begin(), TargetArray.end());
    uvec_t::iterator it = std::unique_copy(TargetArray.begin(),
            TargetArray.end(), TargetArray.begin());
    TargetArray.resize(it - TargetArray.begin());
}

/*
 * Merge <ID, good_cnt, cnt> tuples
 */
void mergeKInfo(kcvec_t& TargetArray, const kcvec_t& sourceArray){

    if (TargetArray.size() == 0){
        TargetArray = sourceArray;
        return;
    }

    uint64_t oriSize = TargetArray.size();
    /*
     * for every item in sourceArray, identify its position in TargetArray
     * Once identified, merge. Otherwise, add to the end of TargetArray to
     * be sorted later
     */
    for (unsigned int i = 0; i < sourceArray.size(); ++ i){
        std::pair<kcvec_t::iterator, kcvec_t::iterator> bounds =
        std::equal_range (TargetArray.begin(), TargetArray.begin() + oriSize,
                sourceArray[i], Knumcomp());
        if (bounds.first == bounds.second) { // didn't find
            TargetArray.push_back(sourceArray[i]);
        }
        else{
            (*bounds.first).cnt += sourceArray[i].cnt;
            (*bounds.first).goodCnt += sourceArray[i].goodCnt;
        }
    }
    std::sort(TargetArray.begin(), TargetArray.end(), Knumcomp());
}
 

/*
 * Load short reads and quality information
 */
void Parser::load(const Para& myPara) {

    if (!file_exists(myPara.iFaName.c_str())) {
        std::cerr << "open " << myPara.iFaName << "failed, does it exist? \n";
        exit(1);
    }

    std::cerr << "Loading reads (batch size: " << (double) myPara.batchSize / 1000
            << " K) \n\t Batch:";
    int batch = 0;
    double timing = get_time();
    int fetchRetCode = 0;

    Seq mySeq(myPara.iFaName, myPara.iQName);

    while (fetchRetCode >= 0) {

        std::cerr << batch << "\t";
        batch++;

        fetchRetCode = mySeq.retrieve_batch(myPara.batchSize);

        // kmer
        mergeK(kArray_, mySeq, myPara.K);
        // tile
        Kmer tmpTileArray;
        tmpTileArray.getKC(mySeq, myPara.step+myPara.K, myPara);
        mergeKInfo(tileArray_, tmpTileArray.getElem());
        
        mySeq.clear();
        print_time("", timing);
        
    }

    std::cerr << "Loading done! Approximately "
            << (double) batch * myPara.batchSize / 1000000 << " Million Reads. \n\n";

}


int bisearch(const kcvec_t& KArray, int size, uint64_t elem) {
    int lb = 0, ub = size - 1, mid;
    while (lb <= ub) {
        mid = (lb + ub) / 2;
        if (KArray[mid].ID == elem)
            return mid;
        else if (KArray[mid].ID < elem)
            lb = mid + 1;
        else if (KArray[mid].ID > elem)
            ub = mid - 1;
    }
    return -1;
}


void print(int i) {
    std::cerr << " " << i;
}

void printHex(int i) {
    std::cerr << " " << std::hex << i;
}

void Parser::tableMaker(const Para& myPara) {

    std::cerr << "Constructing Facility Tables ... \n";
    /*
     * 1. split K positions into random chunks based
     *    on myPara.k and myPara.eSearch
     */
    ivec_t indices(myPara.K, 0);
    for (int i = 0; i < myPara.K; ++i) indices[i] = i;
    std::random_shuffle(indices.begin(), indices.end());
    //std::for_each(indices.begin(), indices.end(), print);  std::cerr <<"\n";

    int unitSize = (log2(myPara.eSearch) / 2);
    int num = myPara.K / unitSize;
    
    for (int i = 0; i < num; ++i) {
        ivec_t tmp(unitSize, 0);
        std::copy(indices.begin() + i*unitSize, indices.begin() + (i + 1)*unitSize, tmp.begin());
        //tmp.insert(tmp.end(), indices.begin() + i*unitSize,
        //        indices.begin() + (i + 1) * unitSize);
        maskIdx_.push_back(tmp);
    }

    int lastUnitSize = myPara.K % unitSize;
    if (lastUnitSize != 0) {
        ivec_t tmp;
        tmp.insert(tmp.end(), indices.begin() + num*unitSize, indices.end());
        maskIdx_.push_back(tmp);
        num++;
    }

    //std::cerr << "\t" << num << " Tables to be created\n";

    /*
     * 2. Create Masks for tables
     */
    masks_.resize(num);
    for (int i = 0; i < num; masks_[i] = 0xFFFFFFFF, i++);
    for (int i = 0; i < num; ++i) {
        for (int j = 0; j < maskIdx_[i].size(); ++j) {
            masks_[i] &= ~(1 << (2 * maskIdx_[i][j]));
            masks_[i] &= ~(1 << (2 * maskIdx_[i][j] + 1));
        }
    }

    //std::for_each(masks_.begin(), masks_.end(), printHex);
    //std::cerr << "\n";

    /*
     * 3. Create Tables Then sort according to masks for each table
     */    
    for (int i = 0; i < num; ++ i){
        table_.insert(table_.end(), kArray_.begin(), kArray_.end());
    }
    /*
     * sort tables
     */
    int tableSize = kArray_.size();
    for (int i = 0; i < num; ++i) {
        std::sort(table_.begin() + i*tableSize,
                table_.begin() + (i + 1) * tableSize, TComp(masks_[i]));
    }

    std::cerr << "\tdone !\n\n";

}

void Parser::ec(const Para& myPara) {

    std::cerr << "Start Error Correction in batches...(batch size: "
            << (double) myPara.batchSize / 1000 << " K) \n\t Batch:";

    if (file_exists(myPara.iFaName.c_str()) == false) {
        std::cerr << "open " << myPara.iFaName << "failed :|\n";
        exit(1);
    }

    int batch = 0;
    double timing = get_time();
    int fetchRetCode = 0;

    std::ofstream *pHandle = 0;
    if(myPara.oErrName == "-")
        pHandle = new std::ofstream(myPara.oErrName.c_str());
    Seq mySeq(myPara.iFaName, myPara.iQName);

    while (fetchRetCode >= 0) {

        std::cerr << batch << "\t";
        batch++;

        fetchRetCode = mySeq.retrieve_batch(myPara.batchSize);

        //---------------------------------------------------------------------
        // now start correcting this batch -- could be multi-threaded
        //---------------------------------------------------------------------
        int numSeq = mySeq.getSL().size();

        for (int i = 0; i < numSeq; ++i) {
            int pos = mySeq.getSL()[i];
            readID_ = i; //
            char* addr = const_cast<char*> (&mySeq.getS()[pos]);
            char* qAddr = NULL;
            if (myPara.QFlag) {
                qAddr = const_cast<char*> (&mySeq.getQ()[pos]);
            }
            readEC(addr, qAddr, myPara);
            // readID_++;
        }// for
        //---------------------------------------------------------------------                
        output(mySeq, (pHandle == 0) ? std::cout : *pHandle);
        print_time("", timing);
    }//while (fasta_sr.operator bool())

    if(pHandle != 0){
        pHandle->close();
        delete pHandle;
    }
    std::cerr << "\tdone !\n\n";
}

void Parser::readEC(char* addr, char* qAddr, const Para& myPara) {


    ipair_t hdUb(myPara.hdMax, myPara.hdMax);

    ipair_t dPoints(0, myPara.K - myPara.step); //dividing points //a
    
    std::string myRead = addr;
    int readLen = myRead.length();
    int nextPos = 0, prePos = -1;
    bool crawlingFlag = false;
    int crawlTo = -1;

    while (1) {

        upair_t mosaic;
        uint64_t ID1, ID2;
        if (!(toID(ID1, &myRead[nextPos], myPara.K) &&
                toID(ID2, &myRead[nextPos + myPara.step], myPara.K)))
            return; /* containing non-acgt*/
        mosaic = upair_t(ID1, ID2);

        if (errorCall(mosaic, dPoints, hdUb, qAddr + nextPos, myPara)) {//a
 
            updateRead(myRead, nextPos);

            if (nextPos + myPara.step + myPara.K == readLen) return;

            prePos = nextPos;

            if (!crawlingFlag) {
                if (nextPos + myPara.step <= readLen - myPara.K - myPara.step) {
                    dPoints = ipair_t(myPara.K, 0);
                    nextPos += myPara.step;
                    hdUb = ipair_t(0, myPara.hdMax);
                }
                else {
                    dPoints = ipair_t(myPara.K, nextPos + myPara.K + 2 * myPara.step - readLen);
                    nextPos = readLen - myPara.K - myPara.step;
                    hdUb = ipair_t(0, myPara.hdMax);
                }
            }
        } 
        else {
            if (crawlingFlag == false){
                crawlingFlag = true;
                crawlTo = nextPos;
            }
        }

        if (crawlingFlag){

         // crawling --- to do reverse crawling
            if (nextPos == prePos + 1 || prePos == -1)  return;

            nextPos = prePos + 1;

            if(nextPos == crawlTo) crawlingFlag = false;
            hdUb = ipair_t(0, 1);
            dPoints = ipair_t(myPara.K, myPara.K - 3);
        }
    }
}

void Parser::updateRead(std::string& myRead, int shift) {

    if (readErr_.size()) {        
        // 1. adding shift to err positions w.r.t read & update read
        for (int i = 0; i < readErr_.size(); ++i) {
            readErr_[i].pos += shift;
            myRead[readErr_[i].pos] = bits_to_char(readErr_[i].from);
        }
        // 2. update records_
        if (records_.size() && records_.back().readID == readID_) {
            (records_.back()).evec.insert((records_.back()).evec.end(),
                    readErr_.begin(), readErr_.end());
        } else {
            records_.push_back(record_t(readID_, readErr_));
        }

        // 3. update snKArray_ and lnKArray_ --- todo

        readErr_.clear();
    }
}

void candidates(ivec_t& candies, const kcvec_t& tiles, int threshold) {
    for (int i = 1; i < tiles.size(); ++i) {
        if (tiles[i].goodCnt >= threshold)
            candies.push_back(i);
    }
}

void diff(evec_t& errs, uint64_t from, uint64_t to, char* qAddr, int len) {
    uint64_t tmp = from ^ to;
    int idx = 0;
    while (tmp) {
        if (tmp & 0x3) {
            e_t err;
            err.pos = len - 1 - idx;
            err.from = from & 0x3;
            err.to = to & 0x3;
            err.qual = (int) qAddr[(len - 1) - idx];
            errs.push_back(err);
        }
        tmp >>= 2;
        from >>= 2;
        to >>= 2;
        idx++;
    }
}

bool goodQuals(char* qAddr, int len, int threshold) {
    for (int i = 0; i < len; ++i) {
        if (qAddr[i] < threshold)
            return false;
    }
    return true;
}


bool Parser::errorCall(const upair_t& mosaic, const ipair_t& dPoints,
        const ipair_t& hdUb, char* qAddr, const Para& myPara) {

    upair_t rvMosaic(
            reverse_complementary<uint32_t, uint32_t > (mosaic.first, myPara.K),
            reverse_complementary<uint32_t, uint32_t > (mosaic.second, myPara.K));
    
    uvec_t heads, tails, rvHeads, rvTails;
    if (std::binary_search (kArray_.begin(), kArray_.end(), mosaic.first)){
        heads.push_back(mosaic.first);
    }
    if (std::binary_search (kArray_.begin(), kArray_.end(), mosaic.second)){
        tails.push_back(mosaic.second);
    }
    if (std::binary_search (kArray_.begin(), kArray_.end(), rvMosaic.second)){
        rvHeads.push_back(rvMosaic.second);
    }
    if (std::binary_search (kArray_.begin(), kArray_.end(), rvMosaic.first)){
        rvTails.push_back(rvMosaic.first);
    }

    ipair_t hd(0, 0);
    while (hd.first <= hdUb.first) {

        kcvec_t tiles, rvTiles;
        tiling(tiles, heads, tails, myPara);
  //      print_kcvec("tiles", tiles, myPara.K + myPara.step);
        tiling(rvTiles, rvHeads, rvTails, myPara);
  //      print_kcvec("rvtiles", rvTiles, myPara.K + myPara.step);
        mergeTiles(tiles, rvTiles, myPara);
  //      print_kcvec("After Merging", tiles, myPara.K + myPara.step);

        /*
         * Error Calling Current Tile
         * tiles.size() may equal to 0 due to error correction
         */
        if (tiles.size() > 0) {
            if (tiles[0].goodCnt >= myPara.tGoodTile) return true;
            else {
                if (tiles.size() == 1) { // try to correct
                    if (tiles[0].goodCnt >= myPara.tCard &&
                            goodQuals(qAddr, myPara.K + myPara.step, myPara.Qlb))
                        return true;
                } else {
                    ivec_t candies; // indices of tiles
                    candidates(candies, tiles, myPara.tCard);

                    if (tiles[0].goodCnt >= myPara.tCard) {
                        /*
                         * ec only if non-ambig correction of low qual bases
                         * could be identified and correcting to one of the
                         * high cardinality neighbors
                         */
                        int tGoodCnt = tiles[0].goodCnt / myPara.tRatio;

                        ivec_t highCardNbs;
                        for (int i = 0; i < candies.size(); ++i) {
                            if (tiles[candies[i]].goodCnt >= tGoodCnt)
                                highCardNbs.push_back(candies[i]);
                        }
                        int alterNum = 0;
                        for (int i = 0; i < highCardNbs.size(); ++i) {
                            /*
                             * check all candidates, if ambig, then do not ec
                             */
                            evec_t errs;
                            diff(errs, tiles[highCardNbs[i]].ID, tiles[0].ID,
                                    qAddr, myPara.step + myPara.K);
                            for (int j = 0; j < errs.size(); ++j) {
                                if (errs[j].qual < myPara.Qlb) {
                                    readErr_ = errs;
                                    alterNum++;
                                    break;
                                }
                            }
                        }
                        if (alterNum > 1) { // ambig
                            readErr_.clear();
                        } else return true;
                    } else {
                        /*
                         * ec only if non-ambig correction to one of [candies]
                         * could be identified
                         */
                        if (candies.size() == 1) {
                            diff(readErr_, tiles[candies[0]].ID, tiles[0].ID,
                                    qAddr, myPara.step + myPara.K);
                            return true;
                        }
                    }
                }
            }
        }
        /*
         * increase Hamming Distance to search for more neighbors
         */
        hd.second++;
        if (hd.second <= hdUb.second) {
            unitNeighbor(tails, dPoints.second, myPara);
            unitNeighbor(rvHeads, dPoints.second, myPara);
        } else {
            hd.first++;
            if (hd.first <= hdUb.first) {
                unitNeighbor(heads, dPoints.first, myPara);
                unitNeighbor(rvTails, dPoints.first, myPara);
            }
            else{ // all possible searches have been done, till this stage
                  // if tile[0] is the maximum and > tCard then, it is considered
                  // as correct due to low coverage region
                int maxCnt = 0;
                bool tflag = true;
                if (tiles.size() > 0) maxCnt = tiles[0].goodCnt;
                if (maxCnt >= myPara.tCard) {
                    for (int i = 1; i < tiles.size(); i ++){
                        if (maxCnt < tiles[i].goodCnt) {
                            tflag = false;
                            break;
                        }
                    }
                    if (tflag) return true;
                }                
            }
        }
    }    
    return false;
}

bool checkPoint(const ivec_t& indices, int dPoint) {
    for (int i = 0; i < indices.size(); ++i) {
        if (indices[i] >= dPoint)
            return true;
    }
    return false;
}

void Parser::unitNeighbor(uvec_t& myNB, int dPoint, const Para& myPara) {

    if (myNB.size() == 0) return;

    uvec_t inputIDs(myNB.size(), 0);
    for (int i = 0; i < myNB.size(); inputIDs.push_back(myNB[i]), ++ i);
    std::sort(inputIDs.begin(), inputIDs.end());
    
    uvec_t nbIDs; // neigbhoring IDs

    for (int i = 0; i < myNB.size(); ++i) {

        uint32_t myID = myNB[i];

        uvec_t myNeighbor; //store ID's neighbors by searching in all tables

        /*
         *  search in every table
         */

        int tableSize = table_.size() / masks_.size();
        for (int j = 0; j < masks_.size(); ++j) {

            if (!checkPoint(maskIdx_[j], dPoint)) continue; //a

            std::pair<uvec_t::iterator, uvec_t::iterator> bounds;

            bounds = std::equal_range(table_.begin() + j*tableSize,
                    table_.begin() + (j + 1) * tableSize, myID, TComp(masks_[j]));

            // identify if HD = 1
            for (uvec_t::iterator it = bounds.first; it != bounds.second; ++it) {
                // check if differed position satisfies dPoint
                int idx = myPara.K - 1;
                bool dPointFlag = true;

                //if (*it == myID) continue;
                if (std::binary_search(inputIDs.begin(), inputIDs.end(), *it)) continue;
                
                uint32_t tmp = (*it) ^ myID;
                int cnt = 0;
                while (tmp) {
                    
                    if (idx < dPoint) {
                        dPointFlag = false;
                        break;
                    }
                    
                    if (tmp & 0x3) cnt++;
                    tmp >>= 2;
                }
                if (cnt == 1 && dPointFlag) {
                    myNeighbor.push_back(*it);
                }
            }
        }

        nbIDs.insert(nbIDs.begin(), myNeighbor.begin(), myNeighbor.end());
    }

    /*
     * remove duplicates
     */
    std::sort(nbIDs.begin(), nbIDs.end());
    uvec_t::iterator it
            = std::unique_copy(nbIDs.begin(), nbIDs.end(), nbIDs.begin());
    nbIDs.resize(it - nbIDs.begin());

    /*
     * now have all nbIDs from table, then search in [snKArray_]
     * to get complete multiplicity and quality score information
     */

    for (int i = 0; i < nbIDs.size(); ++i) {

        if (std::binary_search(kArray_.begin(), kArray_.end(), nbIDs[i])){
            myNB.push_back(nbIDs[i]);
        }
    }
}

/*
 * Keep elements of N only if the ID equals to the last 2*len bits in tiles
 */
void updateNodes(kcvec_t& N, kcvec_t& N_rv, const kcvec_t& tiles, int len) {

    std::set<uint64_t> IDs, IDs_rv;
    kcvec_t tmpVec;
    for (int i = 0; i < tiles.size(); ++i) {
        uint64_t last2k = tiles[i].ID & ((0x1 << 2 * len) - 1);
        IDs.insert(last2k);
        IDs_rv.insert(reverse_complementary <uint32_t, uint32_t > (last2k, len));
    }
    for (int i = 0; i < N.size(); ++i) {
        if (IDs.count(N[i].ID)) tmpVec.push_back(N[i]);
    }
    N = tmpVec;
    tmpVec.clear();
    for (int i = 0; i < N_rv.size(); ++i) {
        if (IDs_rv.count(N_rv[i].ID)) tmpVec.push_back(N_rv[i]);
    }
    N_rv = tmpVec;
}


bool unicpy(const kc_t& e1, const kc_t& e2) {
    return (e1.ID == e2.ID);
}

void Parser::mergeTiles(kcvec_t& tileTo, kcvec_t& tileFrom, const Para& myPara) {
    /*
     *  convert IDs in [tileFrom] to revcompl IDs and insert to [tileTo]
     *  then sort and merge (remove duplicates); keep the repID in place
     */
    if (tileFrom.size() == 0) return;

    /* convert to reverse compl IDs
     */
    for (int i = 0; i < tileFrom.size(); ++i) {
        tileFrom[i].ID = reverse_complementary <uint64_t, uint64_t >
                (tileFrom[i].ID, myPara.K + myPara.step);

    }

    if (tileTo.size() == 0){
        tileTo = tileFrom;
        tileFrom.clear();
        return;
    }

    uint64_t repID = tileTo[0].ID;
    tileTo.insert(tileTo.end(), tileFrom.begin(), tileFrom.end());
    tileFrom.clear();

    /* sort
     */
    std::sort(tileTo.begin(), tileTo.end(), Knumcomp());

    /* merge
     */
    int idx1 = 0;

    for (int idx2 = idx1 + 1; idx2 < tileTo.size(); ++idx2) {
        if (tileTo[idx2].ID == tileTo[idx1].ID) {
            tileTo[idx1].goodCnt += tileTo[idx2].goodCnt;
            tileTo[idx1].cnt += tileTo[idx2].cnt;
        } else {
            idx1 = idx2;
        }
    }

    kcvec_t::iterator it =
            std::unique_copy(tileTo.begin(), tileTo.end(), tileTo.begin(), unicpy);
    tileTo.resize(it - tileTo.begin());

    /* keep repID in the first pos
     */
    if (tileTo[0].ID != repID) {
        for (int i = 1; i < tileTo.size(); ++i) {
            if (tileTo[i].ID == repID) {
                kc_t tmp = tileTo[0];
                tileTo[0] = tileTo[i];
                tileTo[i] = tmp;
                break;
            }
        }
    }
}

/*
 * flag: true -- forward    flase--reverse 
 */
void Parser::tiling(kcvec_t& tiles, const uvec_t& N1,
        const uvec_t& N2, const Para& myPara) {

    if (N1.size() == 0 || N2.size() == 0) return;
    uint64_t reptile;
    if (!overlay(reptile, N1[0], N2[0], myPara)) {
        std::cerr << "Err: errCall, reptile construction fail\n";
        exit(1);
    }
    int idx = bisearch(tileArray_, tileArray_.size(), reptile);
    if (idx != -1) {
        tiles.push_back(tileArray_[idx]);
    }
    else {
        tiles.push_back(kc_t(reptile, 0, 0));
    }
    for (int i = 0; i < N1.size(); ++i) {
        for (int j = 0; j < N2.size(); ++j) {
            if (i == 0 && j == 0) continue;
            uint64_t tmptile;
            if (overlay(tmptile, N1[i], N2[j], myPara)) {
                //binary search 
                idx = bisearch(tileArray_, tileArray_.size(), tmptile);

                if (idx != -1) {
                    tiles.push_back(tileArray_[idx]);
                }
            }
        }
    }
}

bool Parser::overlay(uint64_t& rslt, uint32_t n1, uint32_t n2, const Para& myPara) {

    int shift = 2 * myPara.step;
    uint64_t a = n1, b = n2;
    // comparison suffix - prefix if necessary
    if (myPara.step >= myPara.K || 
            (b >> shift) == (a & ((1 << 2*(myPara.K - myPara.step)) - 1))) {
        rslt = (a << shift) | b;
        return true;
    }
    return false;
}

void Parser::output(Seq& mySeq, std::ostream& oHandle) {

    /*
     * Resulting file format:
     * ReadID (sorted)  ErrNum  [pos from to qual] [pos from to qual]...
     * from: reference; to: read (numercial value Aa:0 Cc:1 Gg:2 Tt:3 others 4)
     * qual: quality (numerical value)
     */
    for (int i = 0; i < records_.size(); ++i) {

        oHandle << mySeq.getHeaders()[records_[i].readID] << "\t" << records_[i].evec.size();

        for (int j = 0; j < records_[i].evec.size(); ++j) {
            oHandle << "\t" << records_[i].evec[j].pos << "\t"
                    << records_[i].evec[j].from << "\t" << records_[i].evec[j].to
                    << "\t" << records_[i].evec[j].qual;
        }
        oHandle << "\n";
    }
    records_.clear();
}

