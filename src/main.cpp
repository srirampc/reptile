
/***
 **
 *  File: main.cpp
 *  Created: Dec 12, 2009 4:05 PM
 *
 *  Author: Xiao Yang <isuxyang@gmail.com>
 * 
 *  Copyright 2009 Xiao Yang, Karin Dorman, Srinivas Aluru
 *
 *  This file is part of Reptile (version 1.1)
 *
 *  Changes made compared to the former version:
 *  -- more efficient ways to record kmer and tiles
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
 *  along with Reptile. If not, see <http://www.gnu.org/licenses/>.
 *
 */
 

#include "util.h"
#include "Parser.h"
/*
 * Note Although we ingore N , we need remove reads containing too many Ns
 */
int main(int argc, char** argv) {

    Para myPara(argc, argv);

    Parser myParser;

    double timing = get_time();
  
    myParser.load(myPara);

    myParser.tableMaker(myPara);
    
    myParser.ec(myPara);
    
    myParser.output(myPara.oErrName);

    print_time("Program Finished Successfully!!\n", timing);
    
    return (EXIT_SUCCESS);
}

