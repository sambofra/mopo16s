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

/*! \file 
	 \brief Classes and functions for single and multi-objective local search.
*/

#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H

#include <random>
#include <vector>

#include "referenceset.hpp"
#include "efficiency.hpp"
#include "pairset.hpp"
#include "parameters.hpp"

using seqan::IupacString;
using seqan::StringSet;


/// Class for storing pair sets encountered duirng local search
class Archive
{
	vector<PairSet> pSets;
	
public:
	/// get archive length
	uint length();
	
	/// add one element, if not duplicated, and returns true if added
	bool add( const PairSet& pSet );
	
	/// retrieve the element with the corresponding index
	PairSet get( uint ind );
	
	/// returns the indices of the elements forming the Pareto Front
	vector<uint> paretoFront(); 
	
	/// get bounds on the three scores from the solutions in the archive
	Bounds getBounds();
};

/*! 
	\brief Single-objective local search
	\param init Initial PairSet
	\param alpha Vector of weights for the three objectives
	\param pars Parameters object containing all required parameters
	\param b Bounds object containing upper and lower bounds for the three scores

	Starting from the PairSet init, performs the sequence of local
	moves (addition, removal or flip of a nucleotide in a primer) each leading
	to the largest decrease in the weighted sum of the three objectives,
	computed with weight vector alpha and score bounds b
*/
PairSet localSearch( PairSet& init, vector<double>& alpha, Parameters& pars, Bounds& b );

/*! 
	\brief Multi-objective local search
	\param ref ReferenceSet object with the reference sequences that should be covered
	\param goodPairs StringSet with good degenerate primer pairs, 
			 alternating forward and reverse primers 
	\param pars Parameters object containing all required parameters
*/
Archive multiObjSearch( ReferenceSet& ref, StringSet<IupacString>& goodPairs, Parameters& pars );

/// Sample uniformly at random three weights that sum to 1
vector<double> sampleAlpha( std::mt19937& rng );

#endif