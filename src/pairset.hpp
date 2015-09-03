#ifndef PAIRSET_H
#define PAIRSET_H

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include "referenceset.hpp"
#include "feasibility.hpp"
#include "parameters.hpp"

using std::vector;
using seqan::IupacString;
using seqan::StringSet;

class Bounds; // forward declaration, defined afterwards

struct Stats{
	uint median;
	uint quantile;
};

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
	
	// constructor from degenerate primers
	PairSet( IupacString degPrFwd, IupacString degPrRev, ReferenceSet& refSet,
	Parameters& params);
	
	// constructor from non-degenerate primer sets
	PairSet( seqan::StringSet<IupacString> prFwd, seqan::StringSet<IupacString> prRev, 
	ReferenceSet& refSet, Parameters& params);	
	
	// getters for the scores
	double getCovSc() const;
	double getVarSc() const;
	double getFeaSc() const;
	
	// getter for the scaled weighted score, to be minimized
	double getWSc( vector<double>& alpha, Bounds& b );
	
	// Length of the primer sets
	uint setLength( bool fwdSet );
	// Length of the p-th forward or reverse primer
	uint prLength( bool fwdSet, uint p );
	
	// Returns the p-th forward or reverse primer
	IupacString getPrimer( bool fwdSet, uint p );
	
	// Nice printout of the primer set
	void printSet();
	
	/* flip the n-th nucleotide of the p-th primer from either 
	fwd or rev to its l-th subsequent letter in lexicographic order 
	and recompute scores */
	void flip( bool fwdSet, uint n, uint p, uint l );
	
	/* add the l-th nucleotide letter in lexicographic order to the 
	p-th primer from either fwd or rev and to either its head or tail  
	and recompute scores */
	void add( bool fwdSet, bool toHead, uint p, uint l );
	
	/* trim either the head or the tail nucleotide from the p-th primer 
	of either fwd or rev and recompute scores */
	void trim( bool fwdSet, bool fromHead, uint p );

};

// To be used by getWSc, to compute the weighted score
class Bounds
{
public:
	Bounds(){} // empty default constructor
	
	Bounds(PairSet& ps) : // set with ps as worse solution ever
		maxCovSc(1), minCovSc( ps.getCovSc() ),
		maxFeaSc(10), minFeaSc( ps.getFeaSc() ),
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