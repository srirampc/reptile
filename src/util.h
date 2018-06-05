/***
 **
 *  File: util.h
 *  Created: Dec 12, 2009 4:05 PM
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

#ifndef _UTIL_H
#define	_UTIL_H

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <complex> //pow
#include <ctype.h>
#include <cmath>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <iomanip>
#include <stdint.h>

#if defined (_MSC_VER)
  #include "getWinTime.hpp"
#else
  #include <sys/time.h>
#endif

typedef std::vector<int> ivec_t;
typedef std::vector<ivec_t> iivec_t;
typedef std::vector<uint32_t> uvec_t;
typedef std::vector<uvec_t> uuvec_t;
typedef std::vector<char> cvec_t;
typedef std::vector<bool> bvec_t;
typedef std::vector<std::string> strvec_t;
typedef std::map<int, int> imap_t;
typedef std::set<int> iset_t;
typedef std::pair<int, int> ipair_t;
typedef std::pair<uint32_t, uint32_t> upair_t;

class Para {
public:
    std::string iFaName;
    int QFlag;
    std::string iQName;
    int batchSize;  //specify number of reads to be loaded
    int K;          
    std::string oErrName; 
    int step;    
    int qualThreshold;
    int Qlb;    
    int maxBadQPerKmer;
    int eSearch;
    int tGoodTile;
    int tCard;
    //int tConst;
    int hdMax; 
    double tRatio;
    
    Para(int argc, char** argv) : argnum(argc), arg(argv) {
        setPara();
    };

private:
    int argnum;
    char** arg;

    void setPara() {

        if (argnum != 2) { std::cerr << "Specify Config file \n"; exit(1); }

        std::string line, s1;
        std::ifstream input(arg[1]);
        if (!input) { std::cerr << "cant open " << arg[1] << "\n"; exit(1);}

        std::istringstream buf;

        while(getline(input, line)){
            buf.clear();
            buf.str(line);

            if (buf >> s1){
                if (s1 == "InFaFile")
                    buf >> iFaName;
                else if (s1 == "QFlag")
                    buf >> QFlag;
                else if (s1 == "IQFile")
                    buf >> iQName;
                else if (s1 == "OErrFile")
                    buf >> oErrName;
                else if (s1 == "BatchSize")
                    buf >> batchSize;
                else if (s1 == "KmerLen")
                    buf >> K;
                else if (s1 == "hd_max")
                    buf >> hdMax;                
                else if (s1 == "Step")
                    buf >> step;
                else if (s1 == "ExpectSearch")
                    buf >> eSearch;
                else if (s1 == "T_ratio")
                    buf >> tRatio;
                else if (s1 == "QThreshold")
                    buf >> qualThreshold;
                else if (s1 == "MaxBadQPerKmer")
                    buf >> maxBadQPerKmer;
                else if (s1 == "Qlb")
                    buf >> Qlb;
                else if (s1 == "T_expGoodCnt")
                    buf >> tGoodTile;
                else if (s1 == "T_card")
                    buf >> tCard;                            
            }
        }
        if (iFaName.empty()) {
            std::cerr << "Err: InFaFile is not specified!\n";
            exit (1);
        }
        std::cerr << "Input Parameters:\n----------------------------------\n";
        std::cerr << "Short Reads file: " << iFaName << "\n";
        if (QFlag){
            std::cerr << "I/QualFile = " << iQName << "\n"
                      << "Qthreshold = " << qualThreshold << "\n";
        }
        if (oErrName.empty()) {
            std::cerr << "Err: Output file is not specified!\n";
            exit (1);
        }
        std::ofstream oHandle(oErrName.c_str());
        if (!oHandle.good()) {
            std::cerr << "open " << oErrName << " failed, correct path?\n";
            exit(1);
        }
        oHandle.close();
        std::cerr << "O/ErrFile = " << oErrName << "\n";
        std::cerr << "(K, step, tile) = " << "(" << K << "," << step << "," << K + step << ")\n"
                  << "BatchSize = " << batchSize << "\n"
                  << "Max Hamming Distance Allowed = " << hdMax << "\n"
                  << "ExpectSearch = " << eSearch << "\n"
                  << "T_ratio = " << tRatio << "\n"
                  << "QThreshold = " << qualThreshold << "\n"
                  << "MaxBadQPerKmer = " << maxBadQPerKmer << "\n"
                  << "Qlb = " << Qlb << "\n"
                  << "T_expGoodCnt = " << tGoodTile << "\n"
                  << "T_card = " << tCard << "\n";
        std::cerr << "----------------------------------\n";

        if (K > 16 || (K + step) > 32) {
            std::cerr << "Set K in the range of (0, 16] and K+step in the range of (2, 32]\n";
            exit(1);
        }
    }
};



inline int char_to_bits(char c) {
    static bool flag = false;
    int cvalue = -1;

    switch (c) {
        case 'A':
        case 'a':
            cvalue = 0;
            break;
        case 'C':
        case 'c':
            cvalue = 1;
            break;
        case 'G':
        case 'g':
            cvalue = 2;
            break;
        case 'T':
        case 't':
            cvalue = 3;
            break;
    }

    if (cvalue != -1) {
        return cvalue;
    } else {
        if(flag == false){
                std::cerr << "Note: non-ACGT character found: " << c
                          << " , all will be ignored\n";
                flag = true;
        }
        return -1;
    }
}
inline char bits_to_char(int value) {

    switch (value) {
        case 0:
            return 'a';
        case 1:
            return 'c';
        case 2:
            return 'g';
        case 3:
            return 't';
        default:
            return 'n';
    }
}

/*
 *  reverse_complementary test code:
 *  uint32_t test = 314324;
 *  std::cerr << toString(test, 10) << "\n";
 *  uint32_t result = reverse_complementary<uint32_t, uint32_t> (test, 10);
 *  std::cerr << toString(result, 10) << "\n";
 *
*/
template <typename Tin, typename Tout>
Tout reverse_complementary (Tin num, int len){
    Tout rv = 0;
    for (int i = 0; i < len; ++ i){
        int cvalue = num & 0x3;
        switch (cvalue){
            case 0:
                rv = (rv << 2) | 0x3;
                break;
            case 1:
                rv = (rv << 2) | 0x2;
                break;
            case 2:
                rv = (rv << 2) | 0x1;
                break;
            case 3:
                rv = rv << 2;
                break;
        }
        num >>= 2;
    }
    return rv;
}

template <typename T>
inline bool toID(T& ID, char* addr, int len) {
    ID = 0;
    for (int i = 0; i < len; ++ i){
        int c = char_to_bits(addr[i]);
        if (c == -1) return false;
	ID  = (ID << 2 | c);
    }
    return true;
}

inline std::string toString(uint64_t ID, int len){

    std::string kmer = "";
    for (int i = 0; i < len; ++ i){
        int last = (ID & 0x3);
        char c;
        switch (last){
            case 0: c = 'a';
            break;
            case 1: c = 'c';
            break;
            case 2: c = 'g';
            break;
            case 3: c = 't';
            break;
        }
        kmer += c;
        ID = ID >> 2;
    }
    std::reverse(kmer.begin(), kmer.end());

    return kmer;
}

inline double get_time() {
      timeval t;
      gettimeofday(&t, 0);
      return t.tv_sec + (0.000001 * t.tv_usec);
} // get_time

inline void print_time (const std::string& msg, double& timing){
    double cur_time = get_time();
    std::cerr << msg << "(" << cur_time - timing << " secs)\n\n";
    timing = cur_time;
}

#endif	/* _UTIL_H */
 

