#ifndef LOCALSEARCH_H
#define LOCALSEARCH_H

#include <random>
#include <vector>

#include "referenceset.hpp"
#include "feasibility.hpp"
#include "pairset.hpp"
#include "parameters.hpp"

using seqan::IupacString;
using seqan::StringSet;


class Archive
{
	vector<PairSet> pSets;
	
public:
	// get archive length
	uint length();
	
	// add one element, if not duplicated, and returns true if added
	bool add( const PairSet& pSet );
	
	// retrieve the element with the corresponding index
	PairSet get( uint ind );
	
	// returns the indices of the elements forming the Pareto Front
	vector<uint> paretoFront(); 
	
	// get bounds on the three scores from the solutions in the archive
	Bounds getBounds();
};

PairSet localSearch( PairSet& init, vector<double>& alpha, Parameters& pars, Bounds& b );

Archive multiObjSearch( ReferenceSet& ref, StringSet<IupacString>& goodPairs, Parameters& pars,
	 uint nRest, uint seed );

vector<double> sampleAlpha( std::mt19937& rng );

#endif