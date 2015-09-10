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

#ifndef PAIRSET_H
#define PAIRSET_H

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include "referenceset.hpp"
#include "efficiency.hpp"
#include "parameters.hpp"

using std::vector;
using seqan::IupacString;
using seqan::StringSet;

class Bounds; // forward declaration, defined afterwards

struct Stats{
	uint median;
	uint quantile;
};

/// Object containing a pair of sets of primer, one of forward and one of revers primers 
class PairSet
{
	ReferenceSet * ref;
	
	// parameters for computing the scores
	Parameters * pars;
	
	vector<vector< int > > mpFwd; // matching position of fwd primers
	vector<vector< int > > mpRev; // matching position of rev primers
	// coverage of all pairs fwd x rev, with amplicon length 
	// if covered or -1 otherwise
	vector<vector< int > > covPairs; 
	
	seqan::StringSet<IupacString> fwd; // forward primers set
	seqan::StringSet<IupacString> rev; // revers primers set
	vector<double> tmFwd; // melting temperatures
	vector<double> tmRev; // melting temperatures
	
	// scores
	double covScore, varScore, feaScore;
	
	// convert a degenerate primer to a primer set
	StringSet<IupacString> degenToSet( IupacString& degen );
	
	// compute initial scores 
	void computeScores();
	
	// update covPairs
	void updateCovPairs();
	
	// compute median and given quantile of positive amplicon lenghts from covPairs
	Stats statsAmpLen( double qtile );
	
public:	
	
	/// constructor from degenerate forward and reverse primers
	PairSet( IupacString degPrFwd, IupacString degPrRev, ReferenceSet& refSet,
	Parameters& params);
	
	/// constructor from non-degenerate forward and reverse primer sets
	PairSet( seqan::StringSet<IupacString> prFwd, seqan::StringSet<IupacString> prRev, 
	ReferenceSet& refSet, Parameters& params);	
	
	// getters for the scores
	double getCovSc() const; /*!< Getter for the coverage score */
	double getVarSc() const; /*!< Getter for the variability score */
	double getEffSc() const; /*!< Getter for the efficiency score */
	
	/// getter for the scaled weighted score, to be minimized
	double getWSc( vector<double>& alpha, Bounds& b );
	
	/// Length of the primer sets
	/*! \param fwdSet true if forward set, false if reverse
	*/
	uint setLength( bool fwdSet );
	
	/// Length of the p-th forward or reverse primer
	/*! \param fwdSet true if forward set, false if reverse
	*/
	uint prLength( bool fwdSet, uint p );
	
	/// Returns the p-th forward or reverse primer
	/*! \param fwdSet true if forward set, false if reverse
	*/
	IupacString getPrimer( bool fwdSet, uint p );
	
	/// Nice printout of the primer set
	void printSet();
	
	/// Flip a nucleotide
	/*! 
		Flip the n-th nucleotide of the p-th primer from either 
		fwd or rev to its l-th subsequent letter in lexicographic order 
		and recompute scores 
	*/
	void flip( bool fwdSet, uint n, uint p, uint l );
	
	/// Add a nucleotide 
	/*! 
		Add the l-th nucleotide letter in lexicographic order to the 
		p-th primer from either fwd or rev and to either its head or tail  
		and recompute scores 
	*/
	void add( bool fwdSet, bool toHead, uint p, uint l );
	
	/// Trim a nucleotide
	/*! 
		Trim either the head or the tail nucleotide from the p-th primer 
		of either fwd or rev and recompute scores 
	*/
	void trim( bool fwdSet, bool fromHead, uint p );

};

/// To be used by PairSet::getWSc, to compute the weighted score
class Bounds
{
public:
	/// empty default constructor
	Bounds(){} 
	
	Bounds(PairSet& ps) : // set with ps as worse solution ever
		maxCovSc(1), minCovSc( ps.getCovSc() ),
		maxFeaSc(10), minFeaSc( ps.getEffSc() ),
		maxVarSc( ps.getVarSc() ), minVarSc(0) 
	{} 
		
	double maxCovSc;
	double minCovSc;
	double maxFeaSc;
	double minFeaSc;
	double maxVarSc;
	double minVarSc;
};

#endif