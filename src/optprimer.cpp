#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <ctime>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/find.h>
#include <seqan/index.h>

#include "referenceset.hpp"
#include "feasibility.hpp"
#include "pairset.hpp"
#include "parameters.hpp"
#include "localsearch.hpp"

using namespace seqan;
using std::string;
/*
void usage();

int main( int argc, char * argv )
{
	Parameters params; // initialized with defaults
	
	// options descriptor 
	static struct option longopts[] = {
		{ "randSeed",						required_argument,	NULL,	's' },
		{ "minPrimerLen",					required_argument,	NULL,	'l' },
		{ "maxPrimerLen",					required_argument,	NULL,	'L' },
		{ "minMeltTemp",  				required_argument,	NULL,	't' },
		{ "minGCCont",						required_argument,	NULL,	'c' },
		{ "maxGCCont", 					required_argument, 	NULL, 'C' },
		{ "maxDimers", 					required_argument, 	NULL, 'D' },
		{ "maxDeltaTm",					required_argument,	NULL,	's' },
		{ "maxHomopLen",					required_argument,	NULL,	'b' },
		{ "maxAmplLenSpanCov",			required_argument,	NULL,	't' },
		{ "maxAmplLenSpanFeas",  		required_argument,	NULL,	'c' },
		{ "maxAlenSpanFeasQtile",		required_argument,	NULL,	'o' },
		{ "minMeltTempInterv", 			required_argument, 	NULL, 'C' },
		{ "minGCContInterv", 			required_argument, 	NULL, 'i' },
		{ "maxDimersInterv",				required_argument,	NULL,	's' },
		{ "deltaTmInterv",				required_argument,	NULL,	'b' },
		{ "maxHomopLenInterv",			required_argument,	NULL,	't' },
		{ "maxAmplLenSpanFeasInterv",	required_argument,	NULL,	'c' },
	};
	
	
	uint randSeed = 0, uint minPrimerLen = 17, uint maxPrimerLen = 21,
			uint minMeltTemp = 52, double minGCCont = 0.5, double maxGCCont = 0.7, 
			uint maxDimers = 8, double maxDeltaTm = 3, uint maxHomopLen = 4,
			uint maxAmplLenSpanCov = 200, uint maxAmplLenSpanFeas = 100,
		  	double maxAlenSpanFeasQtile = 0.01, uint minMeltTempInterv = 2, 
			double minGCContInterv = 0.1, double maxDimersInterv = 3, 
			double deltaTmInterv = 2, uint maxHomopLenInterv = 2,
	    	uint maxAmplLenSpanFeasInterv = 50
	
	
	// process options 
	while ((ch = getopt_long(argc, argv, "s:b:t:c:o:C:i:Svhm", longopts, NULL)) != -1)
		switch (ch) 
		{
			case 's': 
				seed = atol(optarg);
				break;
			case 'b': 
				B = atol(optarg);
				break;
			case 't': 
				corrTh = atof(optarg);
				break;
			case 'c': 
				crossVal = atoi(optarg);
				break;
			case 'o':
				outClassFile = optarg;
				break;
			case 'C':
				covDataFile = optarg;
				break;
			case 'i':
				covInfoFile = optarg;
				break;
			case 'S':
				method = SIM;
				break;
			case 'v':
				verbose = true;
				break;
			case 'm':
				multiple = true;
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
		cerr << "No root file name provided" << endl;
		usage();
		exit( EXIT_FAILURE );
	}
	else
		baseName = argv[1];
	
}

void usage()
{}
	*/

int main( int argc, char* argv[] )
{
	// argv[1] = Reference Set
	// argv[2] = Initial solutions
	// argv[3] = Number of restarts
	// argv[4] = Seed
	
	// Read inputs an create data structures
	int i,j,k,mism,plen,bp,n,p,f,l;
	Parameters pars;
	
	uint t0 = clock();
	std::cout << "Read reference set\n";
	ReferenceSet refSet( argv[1], pars );
	
	std::cout << "Read initial solutions\n";
	SeqFileIn goodFileIn(argv[2]);
	StringSet<IupacString> goodPairs;
	StringSet<CharString> goodPairsNames;
	readRecords( goodPairsNames, goodPairs, goodFileIn );
	
	uint nRest = atoi(argv[3]);
	uint seed = atoi(argv[4]);
	
	std::cout << "Optimisation\n";
	Archive ar = multiObjSearch( refSet, goodPairs, pars, nRest, seed );
	
	std::cout << "Save results\n";
	
	uint t1 = clock();
	std::cout << "Time: " << (t1 - t0)/(double)CLOCKS_PER_SEC << "s\n";
	return 0;
}

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
		std::cout << (ar.get(i)).getCovSc() << "\t" << (ar.get(i)).getFeaSc();
		std::cout << "\t" << (ar.get(i)).getVarSc() << "\n";
		std::cout << i << "\n";
	}
	
	// print pfront
	// std::cout << "\nPareto front\n";
	// vector<uint> pFront = ar.paretoFront();
	for( i = 0; i < pFront.size(); i++ )
	{
		(ar.get(pFront[i])).printSet();
		std::cout << (ar.get(pFront[i])).getCovSc() << "\t" << (ar.get(pFront[i])).getFeaSc();
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
		std::cout << "F: " << ps.getFeaSc() << ", C: " << ps.getCovSc() 
			<< ", V: " << ps.getVarSc() << "\n"; 
	}
	
	// save to file
	std::ofstream origFS("orig.archive.txt");
	origFS << "fscore\tcscore\tvscore\n";
	for( i = 0; i < pFrontOrig.size(); i++ )
	{
		ps = arOrig.get( pFrontOrig[i] );
		origFS << ps.getFeaSc() << "\t" << ps.getCovSc() 
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
		newFS << ps.getFeaSc() << "\t" << ps.getCovSc() 
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
	vector<double> fuzScore = feasTestsPrimer(fwd[0][0],pars);
	std::cout << "Fuzzy feasibility score:";
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
