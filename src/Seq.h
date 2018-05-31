
/***
 **
 *  File: Seq.h
 *  Created: December 5, 2009
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
 

#ifndef _SEQ_H
#define	_SEQ_H

#include "util.h"
#include "fasta_file.hpp"

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

class Seq {

protected:
    cvec_t S_;              // full string
    ivec_t SL_;             // location of each read
    cvec_t Q_;              // quality string
    ivec_t QL_;
    strvec_t headers_;  

	gzFile srfp_;
	kseq_t *srfseq_;
    bIO::FASTA_input *fasta_sr, *fasta_q;
    bool qflag_;

    int retrieve_batch_fa(int max);
    int retrieve_batch_fa(bIO::FASTA_input& handle, int max, bool flag);
    int retrieve_batch_fastq(int max);
    
public:
    Seq(const std::string& srfname, const std::string& qfname);
    virtual ~Seq();
    
    int retrieve_batch(int max);
    void qualHist(const Para& myPara);
    const cvec_t& getS () const {
        return S_;
    }
    const ivec_t& getSL() const {
        return SL_;
    }
    const cvec_t& getQ () const {
        return Q_;
    }
    const ivec_t& getQL() const {
        return QL_;
    }
    const strvec_t& getHeaders() const {
        return headers_;
    }
    void clear() {
        S_.clear(); SL_.clear(); Q_.clear(); QL_.clear(); headers_.clear();
    }
};

#endif	/* _SEQ_H */

