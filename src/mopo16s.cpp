/*
 *  This file is part of the mopo16s program.
 *  Copyright (c) Francesco Sambo <sambofra@dei.unipd.it>
 *
 *  mopo16s is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  mopo16s is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See 
 *  the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>
#include <getopt.h>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/find.h>
#include <seqan/index.h>

#include "referenceset.hpp"
#include "efficiency.hpp"
#include "pairset.hpp"
#include "parameters.hpp"
#include "localsearch.hpp"

using namespace seqan;
using std::string;
using std::cout;
using std::endl;

void usage();

int main( int argc, char * argv[] )
{
	Parameters params; // initialized with defaults
	
	// options descriptor 
	static struct option longopts[] = {
		{ "seed",					required_argument,	NULL,	's' },
		{ "restarts", 				required_argument, 	NULL, 'r' },
		{ "outFileName",			required_argument,	NULL, 'o' },
		{ "maxMismatches",		required_argument,	NULL,	'M' },
		{ "minPrimerLen",			required_argument,	NULL,	'l' },
		{ "maxPrimerLen",			required_argument,	NULL,	'L' },
		{ "minTm",			  		required_argument,	NULL,	'm' },
		{ "minGCCont",				required_argument,	NULL,	'c' },
		{ "maxGCCont", 			required_argument, 	NULL, 'C' },
		{ "maxDimers", 			required_argument, 	NULL, 'D' },
		{ "maxHomopLen",			required_argument,	NULL,	'p' },
		{ "maxDeltaTm",			required_argument,	NULL,	'd' },
		{ "maxALenSpanC",			required_argument,	NULL,	'S' },
		{ "maxALenSpanE",  		required_argument,	NULL,	'e' },
		{ "maxALenSpanEQ",		required_argument,	NULL,	'q' },
		{ "minTmInterv", 			required_argument, 	NULL, 't' },
		{ "minGCContInt", 		required_argument, 	NULL, 'g' },
		{ "maxDimersInt",			required_argument,	NULL,	'i' },
		{ "deltaTmInt",			required_argument,	NULL,	'T' },
		{ "maxHLenInt",			required_argument,	NULL,	'P' },
		{ "maxALenSpanEI",		required_argument,	NULL,	'E' },
		{ "help",					no_argument,			NULL,	'h' }
	};

	// process options
	int ch; 
	while (
		(ch = getopt_long(argc, argv, "s:r:o:M:l:L:m:c:C:D:p:d:S:e:q:t:g:i:T:P:E:h", 
		longopts, NULL)) != -1
		)
		switch (ch) 
		{
			case 's': 
				params.seed = atol(optarg);
				break;
			case 'r': 
				params.rest = atoi(optarg);
				break;
			case 'o':
				params.outFName = optarg;
				break;
			case 'M': 
				params.maxMism = atoi(optarg);
				break;
			case 'l': 
				params.minPLen = atoi(optarg);
				break;
			case 'L': 
				params.maxPLen = atoi(optarg);
				break;
			case 'm':
				params.minTm = atof(optarg);
				break;
			case 'c':
				params.minGC = atof(optarg);
				break;
			case 'C':
				params.maxGC = atof(optarg);
				break;
			case 'D':
				params.maxDim = atoi(optarg);
				break;
			case 'p':
				params.maxHomLen = atoi(optarg);
				break;
			case 'd':
				params.dTm = atof(optarg);
				break;
			case 'S':
				params.maxALenSpanCov = atoi(optarg);
				break;
			case 'e':
				params.maxALenSpanEff = atoi(optarg);	
				break;
			case 'q':
				params.maxALenSpanEffQtile = atof(optarg);
				break;
			case 't':
				params.minTmInt = atof(optarg);	
				break;
			case 'g':
				params.minGCInt = atof(optarg);
				break;
			case 'i':
				params.maxDimInt = atoi(optarg);	
				break;
			case 'T':
				params.dTmInt = atof(optarg);
				break;
			case 'P':
				params.maxALenSpanEffInt = atoi(optarg);	
				break;
			case 'h':
				usage();
				exit( EXIT_SUCCESS );
			default:
				usage();
				exit( EXIT_FAILURE ); 
		}
		
	argc -= (optind-1);
	argv += (optind-1);
		
	string baseName;
		
	if( argc == 1 )
	{  
		std::cerr << "No reference set file name provided" << endl;
		usage();
		exit( EXIT_FAILURE );
	}
	else if ( argc == 2 )
	{
		std::cerr << "No initial primer pairs file name provided" << endl;
		usage();
		exit( EXIT_FAILURE );
	} 

	// read inputs from files	
	uint t0 = clock();
	cout << "Read reference set" << endl;
	ReferenceSet refSet( argv[1], params );
	
	cout << "Read initial primer pairs" << endl;
	SeqFileIn goodFileIn(argv[2]);
	StringSet<IupacString> goodPairs;
	StringSet<CharString> goodPairsNames;
	readRecords( goodPairsNames, goodPairs, goodFileIn );
	
	cout << "Optimisation" << endl;
	Archive ar = multiObjSearch( refSet, goodPairs, params );
	
	vector<uint> pf = ar.paretoFront();
	
	cout << "Saving primer pairs from the final Pareto Front in "; 
	cout << params.outFName << ".primers" << endl;
	std::ofstream outPrimers( params.outFName + ".primers" ); 
	PairSet ps = ar.get( 0 );
	for( uint s = 0; s < pf.size(); s++ )
	{
		ps = ar.get( pf[s] );
		// print forward primers
		for( uint p = 0; p < ps.setLength(true); p++ )
			outPrimers << ps.getPrimer(true, p) << "\t";
		// separator
		outPrimers << "x";
		// print reverse primers
		for( uint p = 0; p < ps.setLength(false); p++ )
			outPrimers << "\t" << ps.getPrimer(false, p);
		outPrimers << endl;
	}
	outPrimers.close();
	
	cout << "Saving scores of the primer pairs in "; 
	cout << params.outFName << ".scores" << endl;
	std::ofstream outScores( params.outFName + ".scores" );
	outScores << "Efficiency\tCoverage\tVariability" << endl;
	for( uint s = 0; s < pf.size(); s++ )
	{
		ps = ar.get( pf[s] );
		outScores << ps.getEffSc() << "\t" << ps.getCovSc() << "\t";
		outScores << ps.getVarSc() << endl;
	} 
	outScores.close();
	
	uint t1 = clock();
	std::cout << "Processing time: " << (t1 - t0)/(double)CLOCKS_PER_SEC;
	std::cout << " seconds" << endl;
}

void usage()
{
	cout << "Copyright (c) 2015 Francesco Sambo, Dept. of Information Engineering," << endl;
	cout << "University of Padova, Italy" << endl << endl;
	cout << "mopo16s V1.0: optimal multi-objective design of forward and reverse primer" << endl; 
	cout << "sets for metagenomics studies." << endl << endl;
	cout << "Usage: mopo16s [OPTIONS] reference_set_file initial_primer_pairs_file" << endl << endl;
	cout << "reference_set_file is a .fasta file containing the reference set of" << endl ;
	cout << "sequences for which the primer are designed." << endl << endl;
	cout << "initial_primer_pairs_file is a .fasta file containing a set of (possibly" << endl;
	cout << "degenerate) primer pairs from which to start the optimisation, saved" << endl; 
	cout << "alterning the forward and its corresponding reverse primer." << endl << endl;
	cout << "Common options:" << endl << endl;
	cout << "  -s, --seed=LONG             Seed of the random number generator (default 0)" << endl << endl;
	cout << "  -r, --restarts=INT          Number of restarts of the multi-objective" << endl;
	cout << "                              optimisation algorithm (default 20)" << endl << endl;
	cout << "  -o, --outFileName=FNAME     Root name of the output files (default \"out\")" << endl << endl;
	cout << "  -h, --help                  Print this help and exit" << endl << endl;
	cout << "Coverage-related options:" << endl << endl;
	cout << "  -M, --maxMismatches=INT     Maximum number of mismatches between the" << endl;
	cout << "                              non-3\'-end of the primer and a 16S sequence to" << endl;
	cout << "                              consider the latter covered by the primer, in" << endl;
	cout << "                              case also the 3\'-end perfectly matches" << endl;
	cout << "                              (default 2)" << endl << endl;
	cout << "  -S, --maxALenSpanC=INT      Maximum amplicon length span considered when" << endl;
	cout << "                              computing coverage (half above, half below " << endl;
	cout << "                              median) (default 200)" << endl << endl;
	cout << "Efficiency-related options:" << endl << endl;
	cout << "  -l, --minPrimerLen=INT      Minimum primer length (default 17)" << endl << endl;
	cout << "  -L, --maxPrimerLen=INT      Maximum primer length (default 21)" << endl << endl;
	cout << "  -m, --minTm=INT             Minimum primer melting temperature (default 52)" << endl << endl; 
	cout << "  -c, --minGCCont=DOUBLE      Minimum primer GC content (default 0.5)" << endl << endl; 
	cout << "  -C, --maxGCCont=DOUBLE      Maximum primer GC content (default 0.7)" << endl << endl; 
	cout << "  -D, --maxDimers=INT         Maximum number of self-dimers, ie of dimers" << endl; 
	cout << "                              between all possible gap-less alignments of the" << endl;
	cout << "                              primer with its reverse complement (default 8)" << endl << endl;
	cout << "  -p, --maxHomopLen=INT       Maximum homopolymer length (default 4)" << endl << endl;
	cout << "  -d, --maxDeltaTm=INT        Maximum span of melting temparatures for the" << endl;
	cout << "                              primer sets (default 3)" << endl << endl;
	cout << "  -e, --maxALenSpanE=INT      Maximum span (maxALenSpanE) between median and" << endl;
	cout << "  -q, --maxALenSpanEQ=DOUBLE  given quantile (maxALenSpanEQ) of amplicon" << endl;
	cout << "                              length (default 50 and 0.01, respectively)" << endl << endl;
	cout << "Fuzzy tolerance intervals for efficiency-related options:" << endl << endl;
	cout << "  -t, --minTmInterv=INT       Fuzzy tolerance interval for minimum melting" << endl;
	cout << "                              temperature (default 2)" << endl << endl;
	cout << "  -g, --minGCContInt=DOUBLE   Fuzzy tolerance interval for minimum GC" << endl;
	cout << "                              content (default 0.1)" << endl << endl;
	cout << "  -i, --maxDimersInt=INT      Fuzzy tolerance interval for maximum number of" << endl;
	cout << "                              self dimers (default 3)" << endl << endl;
	cout << "  -T, --deltaTmInt=INT        Fuzzy tolerance interval for span of melting" << endl;
	cout << "                              temperatures of the primer set (default 2)" << endl << endl;
	cout << "  -P, --maxHLenInt=INT        Fuzzy tolerance interval for maximum" << endl;
	cout << "                              homopolymer length (default 2)" << endl << endl;
	cout << "  -E, --maxALenSpanEI=INT     Fuzzy tolerance interval for maximum span" << endl;
	cout << "                              between median and given quantile amplicon" << endl; 
	cout << "                              length (default 50)"<< endl << endl;
	cout << "Mandatory arguments to long options are also mandatory for any corresponding" << endl;
	cout << "short options." << endl << endl;
}


/*int main( int argc, char* argv[] )
{
	// argv[1] = Reference Set
	// argv[2] = Initial solutions
	// argv[3] = Number of restarts
	// argv[4] = Seed
	
	// Read inputs an create data structures
	int i,j,k,mism,plen,bp,n,p,f,l;
	Parameters pars;
	pars.rest = atoi(argv[3]);
	pars.seed = atoi(argv[4]);
	
	uint t0 = clock();
	std::cout << "Read reference set\n";
	ReferenceSet refSet( argv[1], pars );
	
	std::cout << "Read initial solutions\n";
	SeqFileIn goodFileIn(argv[2]);
	StringSet<IupacString> goodPairs;
	StringSet<CharString> goodPairsNames;
	readRecords( goodPairsNames, goodPairs, goodFileIn );
	
	std::cout << "Optimisation\n";
	Archive ar = multiObjSearch( refSet, goodPairs, pars );
	
	std::cout << "Save results\n";
	
	uint t1 = clock();
	std::cout << "Time: " << (t1 - t0)/(double)CLOCKS_PER_SEC << "s\n";
	return 0;
}*/

/*int main( int argc, char* argv[] )
{
	// Read inputs an create data structures	
	IupacString s("ACGTWN");
	std::vector< IupacString > fwd;
	fwd.push_back(s);
	std::cout << fwd[0] << "\n";
	
	StringSet<CharString> id;
	StringSet<IupacString> allIn;

	SeqFileIn seqFileIn("data/good_pairs.fa");
	readRecords(id, allIn, seqFileIn);
	
	int i = 0;
	for( i = 0; i < length(id); i++ )
	{
		std::cout << id[i] << "\t" << allIn[i] << "\n";
	}
	int l = length(id);
	std::cout << id[l-1] << "\t" << allIn[l-1] << "\n";
	
	StringSet<IupacString> nonDegen = degenToSet(allIn[l-1]);
	for( i = 0; i < length(nonDegen); i++ )
		std::cout << nonDegen[i] << "\n";
	
	ReferenceSet refSet( argv[1], 2 );
	Archive ar;
	
	// load solutions and plug them into archive
	string baseName = "/Users/sambofra/dottorato/data/16s/res/multi.sols.85.quant.";
	string name;
	std::ifstream ifs;
	for( i = 1; i < 21; i++ )
	{
		std::cout << baseName + std::to_string(i) + ".txt" << "\n";
		ifs.open(baseName + std::to_string(i) + ".txt", std::ifstream::in);
		string line, fact;
		while( getline(ifs,line) )
		{
			//std::cout << line << "\n";
			std::istringstream iss(line);
			StringSet<IupacString> fwd;
			iss >> fact;
			//std::cout << fact << "\n";
			while( fact.compare("x") != 0 ) // fwd primers
			{
				seqan::appendValue(fwd, IupacString(fact));
				iss >> fact;
				//std::cout << fact << "\n";
			}
			StringSet<IupacString> rev; // rev primers
			iss >> fact;
			while( iss )
			{
				seqan::appendValue(rev, IupacString(fact));
				iss >> fact;
			}
			// add PairSet to Archive
			ar.add( PairSet( fwd, rev, refSet, pars ) );
		}
		ifs.close();
	}
	
	std::cout << ar.length() << "\n";
	
	// print all
	std::cout << "\nEntire archive\n";
	for( i = 0; i < ar.length(); i++ )
	{
		(ar.get(i)).printSet();
		std::cout << (ar.get(i)).getCovSc() << "\t" << (ar.get(i)).getEffSc();
		std::cout << "\t" << (ar.get(i)).getVarSc() << "\n";
		std::cout << i << "\n";
	}
	
	// print pfront
	// std::cout << "\nPareto front\n";
	// vector<uint> pFront = ar.paretoFront();
	for( i = 0; i < pFront.size(); i++ )
	{
		(ar.get(pFront[i])).printSet();
		std::cout << (ar.get(pFront[i])).getCovSc() << "\t" << (ar.get(pFront[i])).getEffSc();
		std::cout << "\t" << (ar.get(pFront[i])).getVarSc() << "\n";
		std::cout << i << "\n";
	}
	
	// Pareto front from original archive

	// read good primer pairs

	std::cout << "Read good pairs\n";
	SeqFileIn goodFileIn(argv[2]);
	StringSet<IupacString> goodPairs;
	StringSet<CharString> goodPairsNames;
	readRecords( goodPairsNames, goodPairs, goodFileIn );
	
	// expand degenerate primers to non degenerate primer sets
	std::cout << "Create and populate archive\n";
	Archive arOrig;
	PairSet ps(goodPairs[0], goodPairs[1], refSet, pars);
	
	for( i = 0; i < length(goodPairs); i+=2 )
	{
		arOrig.add( PairSet(goodPairs[i], goodPairs[i+1], refSet, pars) );
		// ps = ar.get( ar.length() - 1 ); 
	}
	std::cout << "\nOriginal archive length: " << arOrig.length() << "\n"; 
	vector<uint> pFrontOrig = arOrig.paretoFront();
	std::cout << "\nOriginal pareto front length: " << pFrontOrig.size() << "\n"; 
	for( i = 0; i < pFrontOrig.size(); i++ )
	{
		ps = arOrig.get( pFrontOrig[i] );
		ps.printSet();
		std::cout << "F: " << ps.getEffSc() << ", C: " << ps.getCovSc() 
			<< ", V: " << ps.getVarSc() << "\n"; 
	}
	
	// save to file
	std::ofstream origFS("orig.archive.txt");
	origFS << "fscore\tcscore\tvscore\n";
	for( i = 0; i < pFrontOrig.size(); i++ )
	{
		ps = arOrig.get( pFrontOrig[i] );
		origFS << ps.getEffSc() << "\t" << ps.getCovSc() 
			<< "\t" << ps.getVarSc() << "\n"; 
	}
	origFS.close();
	
	// save new pFront to file
	vector<uint> pFront = ar.paretoFront();
	std::ofstream newFS("pfront_true.txt");
	newFS << "fscore\tcscore\tvscore\n";
	for( i = 0; i < pFront.size(); i++ )
	{
		ps = ar.get( pFront[i] );
		newFS << ps.getEffSc() << "\t" << ps.getCovSc() 
			<< "\t" << ps.getVarSc() << "\n"; 
	}
	newFS.close();
	
	//std::cout << "\n" << goodPairsNames[44] << " " << goodPairsNames[45] << "\n";
	//std::cout << goodPairs[44] << "\n" << goodPairs[45] << "\n";
	//pSet.push_back( PairSet(goodPairs[44], goodPairs[45], refSet, pars) );
	//fwd.push_back( degenToSet(goodPairs[44]) );
	//rev.push_back( degenToSet(goodPairs[45]) );
	
	// happy local search (best improvement)
	vector<double> alpha = {1.0/3, 1.0/3, 1.0/3};
	Bounds b( pSet[0] );
	std::cout.precision(4);
	PairSet ps = localSearch( pSet[0], alpha, pars, b );
	ps.printSet();
	std::cout << "First pair, forward:\n";
	for(i = 0; i < length(fwd[0]); i++)
		std::cout << fwd[0][i] << "\n";
	std::cout << "\n";
	std::cout << "Reverse complements:\n";
	for(i = 0; i < length(fwd[0]); i++)
	{
		IupacString str = fwd[0][i];
		reverseComplement(str);
		std::cout << str << "\n";
	}
	std::cout << "\n";
	std::cout << "Melting temperature: " << meltTemp(fwd[0][0]) << "\n";
	std::cout << "GC content: " << fracGCCont(fwd[0][0]) << "\n";
	std::cout << "Last 3 ATs: " << last3AT(fwd[0][0]) << "\n";
	std::cout << "Last 4 CGs: " << last5CG(fwd[0][0]) << "\n";
	std::cout << "Longest homopolymer: " << maxHomLength(fwd[0][0]) << "\n";
	std::cout << "Maximum self-dimers: " << dimers(fwd[0][0],fwd[0][0],8) << "\n";
	std::cout << "Hairpin: " << noHairpins(fwd[0][0]) << "\n";
	vector<double> fuzScore = effTestsPrimer(fwd[0][0],pars);
	std::cout << "Fuzzy efficiency score:";
	for (i = 0; i < 7; i++)
		std::cout << " " << fuzScore[i];  
	std::cout << "\n";
	std::cout << fwd[22][0] << "\n";
	for (i = 0; i < length(fwd[22]); i++)
		std::cout << " " << noHairpins(fwd[22][i]);  
	std::cout << "\n";
	for (i = 0; i < length(rev[22]); i++)
		std::cout << " " << noHairpins(rev[22][i]);  
	std::cout << "\n";
	
	std::cout << "First pair, reverse:\n";
	for(i = 0; i < length(rev[0]); i++)
		std::cout << rev[0][i] << "\n";
	std::cout << "\n";
	
	std::cout << "Reverse complements:\n";
	for(i = 0; i < length(rev[0]); i++)
	{
		IupacString str = rev[0][i];
		reverseComplement(str);
		std::cout << str << "\n";
	}
	std::cout << "\n";
	
	std::cout << "Melting temperature: " << meltTemp(rev[0][0]) << "\n";
	std::cout << "GC content: " << fracGCCont(rev[0][0]) << "\n";
	std::cout << "Last 3 ATs: " << last3AT(rev[0][0]) << "\n";
	std::cout << "Last 4 CGs: " << last5CG(rev[0][0]) << "\n";
	std::cout << "Longest homopolymer: " << maxHomLength(rev[0][0]) << "\n";	
	std::cout << "Maximum self-dimers: " << maxDimers(rev[0][0],rev[0][0],8) << "\n";
	std::cout << "Hairpin: " << hairpin(rev[0][0]) << "\n";	
	std::cout << "\n";
	
	std::cout << "Dimers fwd/rev: " << maxDimers(fwd[0][0],rev[0][0],8) << "\n\n";
	// try some search
	uint t0 = clock();
	for( i = 0; i < length(fwd); i++ )
	{
		std::cout << "fwd " << i << "\n";
		refSet.matchPosFwd( fwd[i][0] );
	}
	
	uint t1 = clock();
	for( i = 0; i < length(rev); i++ )
	{
		std::cout << "rev " << i << "\n";
		refSet.matchPosRev( rev[i][0] );
	}
	
	uint t2 = clock();
	std::cout << "Time for fwd search was: " << (t1 - t0)/(double)CLOCKS_PER_SEC << "\n";
	std::cout << "Time for rev search was: " << (t2 - t1)/(double)CLOCKS_PER_SEC << "\n";
	
	// launch multi-objective optimisation
	return 0;
}*/
