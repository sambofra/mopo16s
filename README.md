## mopo16s: Multi-Objective Primer Optimisation for 16S experiments
Copyright (c) 2015 Francesco Sambo, Dept. of Information Engineering,
University of Padova

Version 1.0

### Installation instructions

Either use git to directly retrieve the source

	git clone https://github.com/sambofra/mopo16s

or download and decompress the archive.

mopo16s requires the SeqAn c++ library: set the correct include directory 
of your SeqAn installation in the Makefile.rules file and then compile with

	make

The executable can then be found in the release folder.
 
To build the html documentation (requires the Doxygen software):

	make doc

mopo16s requires:
- a reference set of 16s sequences for which to design the optimal set 
  of primers and
- a set of (possibly degenerate) candidate primer pairs from which to 
  begin the search for the optimal primer set.

Both files should be in .fasta format, with the latter saved alternating 
forward and corresponding reverse primers.

To generate the two files, rep_set.fasta and good_pairs.fasta, 
starting from the same repositories used in the mopo16s paper, 
first set the required GreenGeenes cluster similarity level for the 
reference set in the shell script download_data.sh, then run:

	cd data	
	./download_data.sh

Then edit the R script process_data.R to set the required GreenGeenes 
cluster similarity level for the reference set and the minimum and maximum 
desired amplicon length, minimum and maximum primer length and maximum 
total degeneracy for the candidate primer pairs.
Then, run the R script with

	Rscript process_data.R

to generate the files rep_set.fasta and good_pairs.fasta (requires R and 
the additional R packages Biostrings and xlsx).
