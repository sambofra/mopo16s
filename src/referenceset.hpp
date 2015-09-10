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

/*! \file */

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

/// Object containing all reference sequences to be covered by the primer pair sets
class ReferenceSet 
{
	seqan::StringSet<seqan::IupacString> refSet;
	seqan::StringSet<seqan::CharString> refSetNames;
	vector<vector<vector<uint> > > cache;
	uint mMism;
	uint len;
	
	// computes an index to the cache from a 5 nt string
	uint hash( seqan::IupacString& tpEnd ); 

public:
	/// Constructor, reads fasta file
	ReferenceSet( char* fastaName, Parameters& pars );

	/// easy getter
	uint length();
	
	/// returns a vector of matching positions for the given 
	/// forward primer, or -1 in case of no match
	vector<int> matchPosFwd( seqan::IupacString& primer );
	
	/// returns a vector of matching positions for the given 
	/// reverse primer, or -1 in case of no match
	vector<int> matchPosRev( seqan::IupacString& primer );
	
	//seqan::IupacString& getReference( uint index );
};

#endif