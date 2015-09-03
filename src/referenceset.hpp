#ifndef REFERENCESET_H
#define REFERENCESET_H

#include <vector>
#include <cmath>
#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include "parameters.hpp"

using std::vector;

class ReferenceSet 
{
	seqan::StringSet<seqan::IupacString> refSet;
	seqan::StringSet<seqan::CharString> refSetNames;
	vector<vector<vector<uint> > > cache;
	uint mMism;
	uint len;
	
	// compute an index to the cache from a 5 nt string
	uint hash( seqan::IupacString& tpEnd ); 

public:
	// constructor, read fasta file
	ReferenceSet( char* fastaName, Parameters& pars );

	// easy getter
	uint length();
	
	// returns a vector of matching positions for the given 
	// forward primer, or -1 in case of no match
	vector<int> matchPosFwd( seqan::IupacString& primer );
	
	// returns a vector of matching positions for the given 
	// reverse primer, or -1 in case of no match
	vector<int> matchPosRev( seqan::IupacString& primer );
	
	seqan::IupacString& getReference( uint index );
};

#endif