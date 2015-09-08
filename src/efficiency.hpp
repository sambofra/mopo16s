/*
 *  This file is part of the optprimer program.
 *  Copyright (c) Francesco Sambo <sambofra@dei.unipd.it>
 *
 *  optprimer is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  optprimer is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See 
 *  the GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/*! 
	\file
	\brief Set of efficiency functions for a single primer
*/

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

/// Melting temperature
double meltTemp( IupacString& primer );

/// Fraction of GC content
double fracGCCont( IupacString& primer );

/// Not all As and Ts in the last three nucleotides
bool last3AT( IupacString& primer );

/// Less than 4 Cs or Gs in the last five nts
bool last5CG( IupacString& primer );

/// Length of the longest homopolymer
uint maxHomLength( IupacString& primer );
	
/// Maximum number of dimers between all possible alignments of two primers (with a lower threshold th)
uint dimers( IupacString& p1, IupacString& p2, uint th );
	
/// No hairpins, ie self alignments of the last nt + at least 3 out of the 4 previous nts
bool noHairpins( IupacString& primer );

/// Fuzzy efficiency tests of a single primer
vector<double> effTestsPrimer( IupacString& primer, Parameters& pars );
	
#endif