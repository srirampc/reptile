
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
#include <zlib.h>
#include "kseq.h"
#include "fasta_file.hpp"

Seq::Seq(const std::string& srfname, const std::string& qfname) {
    if(qfname.length() == 0){
        qflag_ = false;
        srfp_ = gzopen(srfname.c_str(), "r");
        srfseq_ = kseq_init(srfp_);
    } else {
        qflag_ = true;
        fasta_sr = new bIO::FASTA_input(srfname);
        fasta_q = new bIO::FASTA_input(qfname);
    }

}

Seq::~Seq() {
    if(qflag_){
        delete fasta_sr;
        delete fasta_q;
    } else {
        kseq_destroy(srfseq_);
        gzclose(srfp_);
    } 
}

int Seq::retrieve_batch(int max){

    clear();

    if(qflag_)
        return retrieve_batch_fa(max);
    else
        return retrieve_batch_fastq(max);
}

int Seq::retrieve_batch_fastq(int max){
    int l = 0;
    int position = 0;
    int num = 0;

	while ((l = kseq_read(srfseq_)) >= 0) {
        headers_.push_back(srfseq_->name.s);
        SL_.push_back(position);
        S_.resize(S_.size() + srfseq_->seq.l + 1);
        memcpy(&S_[position], srfseq_->seq.s, srfseq_->seq.l + 1);

        if (srfseq_->qual.l > 0){
            QL_.push_back(position);
            Q_.resize(Q_.size() + srfseq_->qual.l + 1);
            memcpy(&Q_[position], srfseq_->qual.s, srfseq_->qual.l + 1);
        } else {
            std::string tmpx(srfseq_->seq.l, 'A');
            QL_.push_back(position);
            Q_.resize(Q_.size() + tmpx.length() + 1);
            memcpy(&Q_[position], tmpx.c_str(), tmpx.length() + 1);
        }

        position += srfseq_->seq.l + 1;
        num++;
        if (num >= max) break;
    }
    for (int i = 0; i < (int) S_.size(); ++ i)
        if(isupper(S_[i])) S_[i] = tolower(S_[i]);
    return l;
}

int Seq::retrieve_batch_fa(int max){
    int sf = retrieve_batch_fa(*fasta_sr, max, 0);
    int qf = retrieve_batch_fa(*fasta_q, max, 1);
    
    return (sf || qf) ? 0 : -1;
}

// flag = 0: short reads file; 1: quality file
int Seq::retrieve_batch_fa(bIO::FASTA_input& handle, int max, bool flag){     

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
    return (handle == true);
}

 