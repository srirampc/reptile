
  - Requirements

  * Perl – for running some preprocessing steps.
  * GNU make – to build Reptile we use Make build system. We tested GNU make.
  * C++ compiler – GNU g++ compiler is recommended.

  - Preprocessing

	  * fastaq-converter.pl
		preprocess all the raw fastq files in an input directory. 

		Usage:

		./fastaq-converter.pl [InDIR] [OutDIR] [Flag]

		Output: 
	
		Each input .fastq file, two output files .fa (reads file) and .q (qual
		file)
 		
  - Reptile (c++)

	  * To compile

		> make all 	 
		
		Alternatively, in the src/ folder, you can use the command to generate exec file: reptile-1.1
		
		> g++ -O3 *.cpp -o reptile-1.1

	  * Configuration file -- "config" need be placed in the same folder as the exe file

	  * Usage (all the parameters used are explained in the config file):
		./reptile config

	  * Output File Format //

		ReadID (sorted)  ErrNum  [pos from to qual] [pos from to qual] ... //
		... ...	//
		
		fields: 
		  *	from: reference; 
		  * to: read (numercial value Aa:0 Cc:1 Gg:2 Tt:3 others 4) 
		  * qual: quality score (ASCII value)
 
  - Post-processing

	  * Currently, a stand-alone script called reptile_merger is provided to generate the fasta format
		output after Reptile correction.

  - Experiments

	  * Mapping short reads to reference genome using RMAP2.04 (http://rulai.cshl.edu/rmap/): 

		./rmap_v2.04/bin/rmap -v -o output.bed -a output.ambig -c ref-genome.fa -m 5 short-reads.fa

	  * Reptile has been applied to the following datasets downloaded from NCBI Short Read Archive

		  * Ref Genome: A. sp (NC_005966)  SRA: SRR006332
          * Ref Genome: E. coli 		   SRA: SRX000429
		  * Ref Genome: E. coli            SRA: SRR001665_1

  - Contact
	  If you have any questions regarding Reptile, please contact Xiao Yang at isuxyang@gmail.com

  - When using Reptile please cite:

	  X. Yang, K. Dorman and S. Aluru, “Reptile: Representative tiling for short read error correction”, Bioinformatics, 2010 (In print). 
