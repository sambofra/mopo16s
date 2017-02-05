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

#include "localsearch.hpp"

uint Archive::length()
{
	return pSets.size();
}

bool Archive::add( const PairSet& pSet )
{
	uint i;
	// check if already present in the archive, based on the scores
	for( i = 0; i < pSets.size(); i++ )
		if( pSets[i].getCovSc() == pSet.getCovSc() && 
			pSets[i].getEffSc() == pSet.getEffSc() &&
			pSets[i].getVarSc() == pSet.getVarSc())
			return false;
	// add
	pSets.push_back(pSet);
	return true;
}

PairSet Archive::get( uint ind )
{
	return( pSets[ind] );
}

vector<uint> Archive::paretoFront()
{
	vector<uint> res;
	uint i,j;
	bool worse;
	for( i = 0; i < pSets.size(); i++ )
	{
		worse = false;
		// search if equal of worse according to all scores
		for( j = 0; j < pSets.size(); j++ )
			if( j != i && pSets[i].getCovSc() <= pSets[j].getCovSc() && 
				pSets[i].getEffSc() <= pSets[j].getEffSc() &&
				pSets[i].getVarSc() >= pSets[j].getVarSc() )
			{
				worse = true;
				break;
			}
		if( !worse )
			res.push_back(i);
	}
	return res;
}

Bounds Archive::getBounds()
{
	uint i;
	Bounds b;
	b.maxCovSc = 0; b.minCovSc = 1;
	b.maxFeaSc = 0; b.minFeaSc = 1;
	b.maxVarSc = 0; b.minVarSc = 1;
	for( i = 0; i < this->length(); i++ )
	{
		if( pSets[i].getCovSc() > b.maxCovSc )
			b.maxCovSc = pSets[i].getCovSc();
		if( pSets[i].getCovSc() < b.minCovSc )
			b.minCovSc = pSets[i].getCovSc();
		
		if( pSets[i].getEffSc() > b.maxFeaSc )
			b.maxFeaSc = pSets[i].getEffSc();
		if( pSets[i].getEffSc() < b.minFeaSc )
			b.minFeaSc = pSets[i].getEffSc();
		
		if( pSets[i].getVarSc() > b.maxVarSc )
			b.maxVarSc = pSets[i].getVarSc();
		if( pSets[i].getVarSc() < b.minVarSc )
			b.minVarSc = pSets[i].getVarSc();
	}
	return b;
}

PairSet localSearch( PairSet& init, vector<double>& alpha, Parameters& pars, Bounds& b )
{
	uint p,n,f,l;
	PairSet curPair = init;
	PairSet bestPair = init;
	PairSet tmpPair = init;
	
	double curWSc = curPair.getWSc( alpha, b );
	double bestWSc = curWSc;
	double tmpWSc;
	
	std::cout.precision(4);
	std::cout << bestPair.getCovSc() << "\t" << bestPair.getEffSc() << "\t";
	std::cout << bestPair.getVarSc() << "\t" << bestWSc << "\n";
		
	bool improve = true;
	while( improve )
	{
		improve = false;
		for( int fwdSet = 1; fwdSet >= 0; fwdSet-- ) // first forward then reverse
		{
			for( p = 0; p < curPair.setLength( fwdSet ); p++ ) // scan primers
			{
				// flip
				for( n = 0; n < curPair.prLength( fwdSet, p ); n++ )
				{
					tmpPair = curPair;
					for( f = 0; f < 3; f++ )
					{
						tmpPair.flip( fwdSet, n, p, 1 );
						tmpWSc = tmpPair.getWSc( alpha, b );
						if( tmpWSc < bestWSc )
						{
							bestPair = tmpPair;
							bestWSc = tmpWSc;
							improve = true;
							std::cout << bestPair.getCovSc() << "\t" << bestPair.getEffSc() << "\t";
							std::cout << bestPair.getVarSc() << "\t" << bestWSc << "\t";
							std::cout << (fwdSet ? "fwd " : "rev ") << p << ": flip\t";
							std::cout << bestPair.getPrimer( fwdSet, p ) << "\n";
						}
					}
				}
				// add and trim
				for( uint head = 0; head < 2; head++ )
				{
					// add
					if( curPair.prLength( fwdSet, p ) < pars.maxPLen )
					{
						for( l = 0; l < 4; l++ )
						{
							tmpPair = curPair;
							tmpPair.add( fwdSet, head, p, l );
							tmpWSc = tmpPair.getWSc( alpha, b );
							if( tmpWSc < bestWSc )
							{
								bestPair = tmpPair;
								bestWSc = tmpWSc;
								improve = true;
								std::cout << bestPair.getCovSc() << "\t" << bestPair.getEffSc() << "\t";
								std::cout << bestPair.getVarSc() << "\t" << bestWSc << "\t";
								std::cout << (fwdSet ? "fwd " : "rev ") << p << ": add\t";
								std::cout << bestPair.getPrimer( fwdSet, p ) << "\n";
							}
						}
					}
					// trim
					if( curPair.prLength( fwdSet, p ) > pars.minPLen )
					{
						tmpPair = curPair;
						tmpPair.trim( fwdSet, head, p );
						tmpWSc = tmpPair.getWSc( alpha, b );
						if( tmpWSc < bestWSc )
						{
							bestPair = tmpPair;
							bestWSc = tmpWSc;
							improve = true;
							std::cout << bestPair.getCovSc() << "\t" << bestPair.getEffSc() << "\t";
							std::cout << bestPair.getVarSc() << "\t" << bestWSc << "\t";
							std::cout << (fwdSet ? "fwd " : "rev ") << p << ": trim\t";
							std::cout << bestPair.getPrimer( fwdSet, p ) << "\n";
						}
					}
				}
				curPair = bestPair;
				curWSc = bestWSc;
			}
		}
	}
	return bestPair;
}

Archive initArchive( ReferenceSet& ref, StringSet<IupacString>& goodPairs, Parameters& pars )
{
   uint i;
   
	// build and populate archive
	Archive ar;
	for( i = 0; i < seqan::length(goodPairs); i+=2 )
		ar.add( PairSet(goodPairs[i], goodPairs[i+1], ref, pars) );
   
   return ar;
}

Archive multiObjSearch( Archive& ar, Parameters& pars )
{
	uint i;
	// initialize rng
	std::mt19937 rng( pars.seed );
	
	// compute Pareto front
	vector<uint> pFront = ar.paretoFront();
	vector<double> alpha = sampleAlpha( rng );
	Bounds b = ar.getBounds();
	PairSet init = ar.get( 0 );
	
	// Cycle, sampling from Pareto front and running local search
	for( i = 0; i < pars.rest; i++ )
	{
		std::cout << "Restart " << i << "\n";
		init = ar.get( pFront[ rng() % pFront.size()] );
		ar.add( localSearch( init, alpha, pars, b ) );
		// update state
		pFront = ar.paretoFront();
		alpha = sampleAlpha( rng );
		b = ar.getBounds();
	}
	
	return ar;
}

vector<double> sampleAlpha( std::mt19937& rng )
{
	vector<double> res(3);
	res[0] = (double)rng()/rng.max();
	res[1] = (double)rng()/rng.max();
	std::sort(res.begin(), res.begin() + 2);
	res[2] = 1.0 - res[1];
	res[1] = res[1] - res[0];
	res[0] = res[0] - 0.0;
	return res;
}