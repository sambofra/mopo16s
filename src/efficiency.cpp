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

#include "efficiency.hpp"

// matrix for calculating melting temperature
const vector<vector<double> > dHMat = {
	{-7.9,-8.4,-7.8,-7.2},
	{-8.5,-8.0,-10.6,-7.8},
	{-8.2,-9.8,-8.0,-8.4},
	{-7.2,-8.2,-8.5,-7.9}};

const vector<vector<double> > dSMat = {
	{-22.2,-22.4,-21.0,-20.4},
	{-22.7,-19.9,-27.2,-21.0},
	{-22.2,-24.4,-19.9,-22.4},
	{-21.3,-22.2,-22.7,-22.2}};
	
double meltTemp( IupacString& primer )
{
   double salt = 50;
   double c_dna = 50;
   double dh, ds;
	
	uint i, oldInd, curInd;
	uint l = seqan::length(primer);
	
	// end
   if( ordValue(primer[l-1]) == 1 || ordValue(primer[l-1]) == 8 )
   {
		dh = 2.3;
		ds = 4.1;
   }
   else
   {
		dh = 0.1;
		ds = -2.8;
   }
	// start
	oldInd = (uint)log2(ordValue(primer[0]));
	if( oldInd == 0 || oldInd == 3 )
	{
      dh += 2.3;
      ds += 4.1;
	}
	else
	{
		dh += 0.1;
		ds += -2.8;
	}
	// neighbours
   for( i = 1; i < l; i++ )
   {
		curInd = (uint)log2(ordValue(primer[i]));
		dh += dHMat[oldInd][curInd];
		ds += dSMat[oldInd][curInd];
		oldInd = curInd;
   }
   double sc = 0.368 * ( l - 1 ) * log ( salt / 1000 );
   double Rk = 1.987 * ( log( 1e-9 * c_dna/4 ) );
   return( 1000*dh / (ds + sc + Rk) - 273.15 );
}


// Fraction of GC content
double fracGCCont( IupacString& primer )
{
	uint count = 0;
	uint l = seqan::length(primer);
	uint i;
	
	for( i = 0; i < l; i++ )
		if( ordValue(primer[i]) == 2 || ordValue(primer[i]) == 4 ) // C or G
			count++;
	
	return( (double)count / l );
}

// Not all As and Ts in the last three nucleotides
bool last3AT( IupacString& primer )
{
	uint l = seqan::length(primer);
	if( 
		(ordValue(primer[l-1]) == 1 || ordValue(primer[l-1]) == 8) &&
		(ordValue(primer[l-2]) == 1 || ordValue(primer[l-2]) == 8) &&
		(ordValue(primer[l-3]) == 1 || ordValue(primer[l-3]) == 8))
		return(false);
	
	return(true);
}

// Less than 4 Cs or Gs in the last five nts
bool last5CG( IupacString& primer )
{
	uint l = seqan::length(primer);
	uint i, count = 0;
	
	for( i = l-5; i < l; i++ )
		if( ordValue(primer[i]) == 2 || ordValue(primer[i]) == 4 )
			count++; 
	return( count < 4 );
}

// Length of the longest homopolymer
uint maxHomLength( IupacString& primer )
{
	uint i, l = seqan::length(primer);
	uint lMaxHom = 1;
	uint lCurHom = 1;
	for( i = 1; i < l ; i++ )
	{
		if(primer[i] == primer[i-1])
			lCurHom++;
		else
		{
			if( lCurHom > lMaxHom )
				lMaxHom = lCurHom;
			lCurHom = 1;
		}
	}
	
	return( lMaxHom );
}

// Maximum number of dimers between all possible alignments
// of two primers (with a lower threshold)
uint dimers( IupacString& p1, IupacString& p2, uint th )
{
	uint i, j, lRC, lBase, curMat, maxMat = 0;
	IupacString base;
	IupacString rc;
	if( seqan::length(p1) < seqan::length(p2) )
	{
		base = p2;
		rc = p1;
		seqan::reverseComplement(rc);
		lRC = seqan::length(p1);
		lBase = seqan::length(p2);
	}
	else
	{
		base = p1;
		rc = p2;
		seqan::reverseComplement(rc);
		lRC = seqan::length(p2);
		lBase = seqan::length(p1);
	}
	
	/* have rc slide over base and count the number of matches */
	// left
	for( i = lRC-th; i > 0; i-- ) // how much of rc is out of base
	{
		curMat = 0;
		for( j = 0; j < lRC - i; j++ )
			if( base[j] == rc[j+i] )
				curMat++;
		if( curMat > maxMat ) // update maximum 
			maxMat = curMat;
	}
	// inside
	for( i = 0; i <= lBase - lRC; i++ ) // right shift of rc
	{
		curMat = 0;
		for( j = 0; j < lRC; j++ )
			if( base[j+i] == rc[j] )
				curMat++;
		if( curMat > maxMat ) // update maximum 
			maxMat = curMat;
	}
	// right
	for( i = lRC - 1; i >= th; i-- ) // how much of rc is in base
	{
		curMat = 0;
		for( j = 0; j < i; j++ )
			if( base[ lBase - i + j] == rc[j] )
				curMat++;
		if( curMat > maxMat ) // update maximum 
			maxMat = curMat;
	}
	return( maxMat );
}

bool noHairpins( IupacString& primer )
{
	uint i, j, count, l = seqan::length(primer);
	IupacString rc = primer;
	reverseComplement(rc);
	
	// 3 prime end of the primer
	for( i = l-1; i > 8; i-- ) // end position 
	{
		// check 3 prime point
		if( rc[i] == primer[l-1] )
		{
			// check the 4 preceding nucleotides
			count = 0;
			for( j = 1; j < 5; j++ ) // left stride 
				if( rc[i-j] == primer[l-1-j] )
					count++;
			if( count > 2 )
				return(false);
		}
	}
	
	// 3 prime end of the rc
	for( i = l-1; i > 8; i-- ) // end position 
	{
		// check 3 prime point
		if( primer[i] == rc[l-1] )
		{
			// check the 4 preceding nucleotides
			count = 0;
			for( j = 1; j < 5; j++ ) // left stride 
				if( primer[i-j] == rc[l-1-j] )
					count++;
			if( count > 2 )
				return(false);
		}
	}
	
	return(true);
}

vector<double> effTestsPrimer( IupacString& primer, Parameters& pars )
{
	/* compute efficiency scores */
	vector<double > res( 7, 1 );
	
	// melting temperature
	double tm = meltTemp(primer);
	if( tm < pars.minTm )
		res[0] = std::max((tm + pars.minTmInt - pars.minTm)/pars.minTmInt, 0.0);
    
	// minGC% < GC < minGC%
	double frGC = fracGCCont( primer );
	if( frGC > pars.maxGC )
		res[1] = 0;
	else if( frGC < pars.minGC )
		res[1] = std::max((frGC + pars.minGCInt - pars.minGC) / 
			(double)pars.minGCInt, 0.0);
    
	// not all As and Ts in last 3 bases and no more than 3 Cs or Gs in last 5
	res[2] = last3AT( primer );
	res[3] = last5CG( primer );
  
	// no homopolymers longer than 4
	uint mHL = maxHomLength( primer );
	if ( mHL > 4 )
		res[4] = std::max((4.0 + pars.maxHomLenInt - mHL)/pars.maxHomLenInt, 0.0);
		
	// no self dimers longer than maxDim
	uint selfDim = dimers( primer, primer, pars.maxDim );
	if( selfDim > pars.maxDim )
		res[5] = std::max((pars.maxDim + pars.maxDimInt - selfDim) / 
			(double)pars.maxDimInt, 0.0);		
  
	// no hairpins  
	res[6] = noHairpins( primer );
	
	return( res );
}
