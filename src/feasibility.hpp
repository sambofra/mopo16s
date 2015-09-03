/* Set of primer feasibility functions */

#ifndef FEASIBILITY_H
#define FEASIBILITY_H

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/modifier.h>
#include <seqan/stream.h>

#include <cmath>
#include <vector>
#include <algorithm>

#include "parameters.hpp"

using std::vector;
using seqan::IupacString;
using seqan::ordValue;

// Melting temperature
double meltTemp( IupacString& primer );

// Fraction of GC content
double fracGCCont( IupacString& primer );

// Not all As and Ts in the last three nucleotides
bool last3AT( IupacString& primer );

// Less than 4 Cs or Gs in the last five nts
bool last5CG( IupacString& primer );

// Length of the longest homopolymer
uint maxHomLength( IupacString& primer );
	
// Maximum number of dimers between all possible alignments
// of two primers (with a lower threshold threshold)
uint dimers( IupacString& p1, IupacString& p2, uint th );
	
// No hairpins, i.e. self alignments of 
// the last nt + at least 3 out of the 4 previous nts
bool noHairpins( IupacString& primer );

// Fuzzy feasibility tests of a single primer
vector<double> feasTestsPrimer( IupacString& primer, Parameters& pars );
	
#endif