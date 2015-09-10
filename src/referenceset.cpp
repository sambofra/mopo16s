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

#include "referenceset.hpp"

ReferenceSet::ReferenceSet( char* fastaName, Parameters& pars )
	:  mMism(pars.maxMism), cache( 1024 ) // 5 ^ 4
{	
	// read from fasta file
	seqan::SeqFileIn repFileIn(fastaName);
	seqan::readRecords( refSetNames, refSet, repFileIn );
	seqan::close( repFileIn );
	len = seqan::length(refSet);
}

uint ReferenceSet::length()
{
	return len;
}

uint ReferenceSet::hash( seqan::IupacString& tpEnd )
{
	uint i, ctr, out = 0;
	
	for( i = 0; i < 5; i++)
		// ordValues for Iupac are powers of 2
		out += pow( 4, i ) * log2( seqan::ordValue( tpEnd[i]) );
	// std::cout << tpEnd << " " << out << "\n";
	return( out );
}

/*
seqan::IupacString& ReferenceSet::getReference( uint index )
{
	return( refSet[index] );
}*/

vector<int> ReferenceSet::matchPosFwd( seqan::IupacString& primer )
{
	uint i, j, k, mism;
	
	// extract three prime end
	uint pLen = seqan::length(primer);
	seqan::IupacString tpEnd = seqan::suffix( primer, pLen-5 );
	
	// search the cache for the match vector corresponding to the tpEnd
	vector<vector<uint> > matchVect;
	uint hashInd = hash( tpEnd ); 
	if( (cache[ hashInd ]).size() > 0 ) // found in the cache
		matchVect = cache[hashInd];
	else
	{
		// search for the three prime end into the entire reference set
		matchVect.resize(length());
		seqan::Pattern<seqan::IupacString, seqan::Horspool> pattern( tpEnd );
		seqan::Finder<seqan::IupacString> finder;
		
		for( i = 0; i < length(); i++ )
		{
			finder = seqan::Finder<seqan::IupacString>(refSet[i]);
			while (seqan::find(finder, pattern))
				matchVect[i].push_back( seqan::beginPosition(finder) );
		}
		
		// save the search to cache
		cache[ hashInd ] = matchVect;
	}
	
	// check the remainder of the primer vs the matches of the tpEnd
	// and save the position of the complete matches, -1 if no match
	vector<int> pos( length(), -1 );
	for( i = 0; i < length(); i++ )
	{
		for( j = 0; j < matchVect[i].size(); j++ )
		{
			// the remainder of the primer should fit
			if( matchVect[i][j] >= pLen - 5 )
			{
				mism = 0;
				for( k = 1; k <= pLen-5; k++ )
					if( primer[pLen-5-k] != refSet[i][matchVect[i][j]-k] )
						mism++;
				if( mism <= mMism )
				{
					pos[i] =  matchVect[i][j] + 5; // first nt after the match
					break;
				}
			}
		}	
	}
	
	return pos;
}

vector<int> ReferenceSet::matchPosRev( seqan::IupacString& primer )
{
	// reverse complement of the primer
	seqan::IupacString rcPrimer = primer;
	seqan::reverseComplement( rcPrimer );
	
	uint i, j, k, mism;
	
	// extract three prime end (left!)
	uint pLen = seqan::length(rcPrimer);
	seqan::IupacString tpEnd = seqan::prefix( rcPrimer, 5 );
	
	// search the cache for the match vector corresponding to the tpEnd
	vector<vector<uint> > matchVect;
	uint hashInd = hash( tpEnd ); 
	if( (cache[ hashInd ]).size()  > 0 ) // found in the cache
		matchVect = cache[hashInd];
	else
	{
		// search for the three prime end into the entire reference set
		matchVect.resize(length());
		seqan::Pattern<seqan::IupacString, seqan::Horspool> pattern( tpEnd );
		seqan::Finder<seqan::IupacString> finder;
		
		for( i = 0; i < length(); i++ )
		{
			finder = seqan::Finder<seqan::IupacString>(refSet[i]);
			while (seqan::find(finder, pattern))
				matchVect[i].push_back( seqan::beginPosition(finder) );
		}
		
		// save the search to cache
		cache[ hashInd ] = matchVect;
	}
	
	// check the remainder of the primer vs the matches of the tpEnd
	// and save the position of the complete matches, -1 if no match
	vector<int> pos( length(), -1 );
	for( i = 0; i < length(); i++ )
	{
		for( j = 0; j < matchVect[i].size(); j++ )
		{
			//std::cout << matchVect[i][j] << "\n";
			// the remainder of the primer should fit
			if( matchVect[i][j] + pLen <= seqan::length( refSet[i] ) )
			{
				mism = 0;
				for( k = 1; k <= pLen-5; k++ )
					if( rcPrimer[4+k] != refSet[i][matchVect[i][j]+4+k] )
						mism++;
				if( mism <= mMism )
				{
					pos[i] =  matchVect[i][j]; // first nt of the match
					break;
				}
			}
		}	
	}
	
	return pos;
	
}