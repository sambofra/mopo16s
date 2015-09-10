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

#include "pairset.hpp"

PairSet::PairSet( IupacString degPrFwd, IupacString degPrRev, ReferenceSet& refSet,
	Parameters& params )
{
	ref = &refSet;
	pars = &params;
	
	// std::cout << "\ndegen\n";
	uint i,j;
	// initialize sets with non-degenerate primers
	fwd = degenToSet( degPrFwd );
	rev = degenToSet( degPrRev );
	/*std::cout << "fwd\n";
	for( i = 0; i < fwdLength(); i++ )
		std::cout << fwd[i] << "\n";
	std::cout << "rev\n";
	for( i = 0; i < revLength(); i++ )
		std::cout << rev[i] << "\n";
	*/	
	// compute coverage of fwd and rev primers
	// std::cout << "coverage ";
	// std::cout << "\n";
	mpFwd.resize( setLength( true ) );
	for( i = 0; i < setLength( true ); i++ )
	{
		mpFwd[i] = ref->matchPosFwd( fwd[i] );
		//for( j = 0; j < ref->length(); j++ )
		//j = 185;
		//	std::cout << mpFwd[i][j] << " ";
		//std::cout << "\n";
	}
	//std::cout << "\n";
	mpRev.resize( setLength( false ) );
	for( i = 0; i < setLength( false ); i++ )
	{
		mpRev[i] = ref->matchPosRev( rev[i] );
		//for( j = 0; j < ref->length(); j++ )
		//j = 185;
		//	std::cout << mpRev[i][j] << " ";
		//std::cout << "\n";
	}
	//std::cout << "\n";
	
	// resize coverage structure
	covPairs.resize(  setLength( true ) * setLength( false ) );
	for( i = 0; i < covPairs.size(); i++ )
		covPairs[i].resize( ref->length() );
	
	// compute melting temperatures
	tmFwd.resize( setLength( true ) );
	for( i = 0; i < setLength( true ); i++ )
		tmFwd[i] = meltTemp( fwd[i] );
	
	tmRev.resize( setLength( false ) );
	for( i = 0; i < setLength( false ); i++ )
		tmRev[i] = meltTemp( rev[i] );	
	
	// compute initial scores
	// std::cout << "scores ";
	computeScores();
}

PairSet::PairSet( seqan::StringSet<IupacString> prFwd, seqan::StringSet<IupacString> prRev, 
ReferenceSet& refSet, Parameters& params) :
fwd( prFwd ), rev( prRev )
{
	ref = &refSet;
	pars = &params;
	
	uint i,j;
	// compute coverage of fwd and rev primers
	// std::cout << "coverage ";
	// std::cout << "\n";
	mpFwd.resize( setLength( true ) );
	for( i = 0; i < setLength( true ); i++ )
	{
		mpFwd[i] = ref->matchPosFwd( fwd[i] );
		//for( j = 0; j < ref->length(); j++ )
		//j = 185;
		//	std::cout << mpFwd[i][j] << " ";
		//std::cout << "\n";
	}
	//std::cout << "\n";
	mpRev.resize( setLength( false ) );
	for( i = 0; i < setLength( false ); i++ )
	{
		mpRev[i] = ref->matchPosRev( rev[i] );
		//for( j = 0; j < ref->length(); j++ )
		//j = 185;
		//	std::cout << mpRev[i][j] << " ";
		//std::cout << "\n";
	}
	//std::cout << "\n";
	
	// resize coverage structure
	covPairs.resize(  setLength( true ) * setLength( false ) );
	for( i = 0; i < covPairs.size(); i++ )
		covPairs[i].resize( ref->length() );
	
	// compute melting temperatures
	tmFwd.resize( setLength( true ) );
	for( i = 0; i < setLength( true ); i++ )
		tmFwd[i] = meltTemp( fwd[i] );
	
	tmRev.resize( setLength( false ) );
	for( i = 0; i < setLength( false ); i++ )
		tmRev[i] = meltTemp( rev[i] );	
	
	// compute initial scores
	// std::cout << "scores ";
	computeScores();
}

StringSet<IupacString> PairSet::degenToSet( IupacString& degen )
{
	// conversion table
	StringSet<IupacString> table;
	seqan::appendValue(table, IupacString("T"));
	seqan::appendValue(table, IupacString("A"));
	seqan::appendValue(table, IupacString("C"));
	seqan::appendValue(table, IupacString("AC"));
	seqan::appendValue(table, IupacString("G"));
	seqan::appendValue(table, IupacString("AG"));
	seqan::appendValue(table, IupacString("CG"));
	seqan::appendValue(table, IupacString("ACG"));	
	seqan::appendValue(table, IupacString("T"));
	seqan::appendValue(table, IupacString("AT"));
	seqan::appendValue(table, IupacString("CT"));
	seqan::appendValue(table, IupacString("ACT"));
	seqan::appendValue(table, IupacString("GT"));
	seqan::appendValue(table, IupacString("AGT"));
	seqan::appendValue(table, IupacString("CGT"));
	seqan::appendValue(table, IupacString("ACGT"));
	
	uint i, j, k, len = length(degen);
	// count the number of non degenerate nucleotides for each base
	uint degenBases[len]; 
	for( i = 0; i < len; i++ )
		degenBases[i] = seqan::length(table[seqan::ordValue(degen[i])]);
	
	// generate all possible combinations of degenerate nucleotides
	StringSet<IupacString> out;
	uint outLen = 1;
	for( i = 0; i < len; i++ )
		outLen *= degenBases[i];
	
	uint indices[len];
	for( i = 0; i < len; i++ )
		indices[i] = 0;
	
	uint incrPtr = 0;
	IupacString primer = degen; 
	for( i = 0; i < outLen; i++ )
	{
		// copy primer indicized by _indices_
		for( j = 0; j < len; j++ )
			primer[j] = table[ seqan::ordValue(degen[j]) ][ indices[j] ];
		seqan::appendValue(out, primer);
		
		// update indices
		indices[incrPtr]++;
		while( (indices[incrPtr] == degenBases[incrPtr]) & (incrPtr <= len) )
			indices[++incrPtr]++;
		
		for( k = 0; k < incrPtr; k++ )
			indices[k] = 0;
		incrPtr = 0;
	}
	return out;			
}

void PairSet::computeScores()
{
	uint i,j;
	double m,sqr;

	// update coverage of each combination of primers
	updateCovPairs();
	
	// std::cout << "median\n";
	
	Stats stats = statsAmpLen( pars->maxALenSpanEffQtile );
	// std::cout << "Median amplicon length: " << medALen << "\n";
	
	covScore = 0;
	varScore = 1;
	if( stats.median > 0 ) // otherwise no coverage at all
	{
		/* compute coverage and variability scores */
		vector<uint> coverage( ref->length(), 0 );
		for( i = 0; i < covPairs.size(); i++ )
			for( j = 0; j < coverage.size(); j++ )
				// check that ampl len lies in the correct span
				if( covPairs[i][j] <= stats.median + pars->maxALenSpanCov/2 && 
					covPairs[i][j] >= stats.median - pars->maxALenSpanCov/2 )
						coverage[j]++;
		
		//for( i = 0; i < ref->length(); i++ )
		//	std::cout << coverage[i] << " ";
		//std::cout << "\n";
	
		// coverage score, fraction of covered references		
		for ( i = 0; i < coverage.size(); i++ )
			if( coverage[i] > 0)
				covScore = covScore + 1.0;
		covScore /= coverage.size();
		
		// variation score, coefficient of variation (sd/m)
		m = 0.0;
		for ( i = 0; i < coverage.size(); i++ )
			m += coverage[i];
		m /= coverage.size();
		sqr = 0.0;
		for ( i = 0; i < coverage.size(); i++ )
			sqr += pow( coverage[i] - m, 2 );
		varScore = sqrt( sqr / (coverage.size() - 1) ) / m;
	}
		
	//std::cout << "Coverage score: " << covScore << "\n";
	//std::cout << "Variation score: " << varScore << "\n";
	
	/* compute efficiency score */
	
	// Fuzzy efficiency tests of the primer pair
	vector<double > eff( 10, 0 );
	
	// melting temperature
	double dTm = std::max( *(std::max_element(tmFwd.begin(), tmFwd.end())), 
			*(std::max_element(tmRev.begin(), tmRev.end())) ) -
		std::min( *(std::min_element(tmFwd.begin(), tmFwd.end())), 
			*(std::min_element(tmRev.begin(), tmRev.end())) );
	
	// std::cout << dTm << "\n";
	 
	if( dTm > pars->dTm )
		eff[7] = std::max((pars->dTm + pars->dTmInt - dTm)/pars->dTmInt, 0.0);
	else
		eff[7] = 1;
	
	// amplicon length span
	if( stats.median > 0 ) // otherwise no coverage at all
	{
		int span = stats.median - stats.quantile;
		if( span > pars->maxALenSpanEff )
			eff[8] = std::max( (double)(pars->maxALenSpanEff + pars->maxALenSpanEffInt - span) /
				 pars->maxALenSpanEffInt, 0.0 );
		else eff[8] = 1;
	}
	else eff[8] = 0;
	
	// dimers
	uint curDim, maxDim = 0;
	for( i = 0; i < setLength( true ); i++ )
		for( j = 0; j < setLength( false ); j++ )
		{
			curDim = dimers( fwd[i], rev[j], pars->maxDim );
			//std::cout << curDim << " ";
			if( curDim > maxDim )
				maxDim = curDim;
		}
	//std::cout << "\n";
	if (maxDim > pars->maxDim)
		eff[9] = std::max( ((double)pars->maxDim + pars->maxDimInt - maxDim) / 
			pars->maxDimInt, 0.0 );
	else eff[9] = 1;
	
	// Average efficiency of the single primers
	vector<double> singleEff;
	for (i = 0; i < setLength( true ); i++ )
	{
		singleEff = effTestsPrimer( fwd[i], *pars );
		for( j = 0; j < 7; j++ )
			eff[j] += singleEff[j];
	}
	for (i = 0; i < setLength( false ); i++ )
	{
		singleEff = effTestsPrimer( rev[i], *pars );
		for( j = 0; j < 7; j++ )
			eff[j] += singleEff[j];
	}
	for( j = 0; j < 7; j++ )
		eff[j] /= (setLength( true ) + setLength( false ));
	
	// Effibiliy score as sum of eff
	feaScore = std::accumulate(eff.begin(), eff.end(), 0.0);
	//for( i = 0; i < 10; i++ )
	//	std::cout << eff[i] << " ";
	//std::cout << "\n";
	//std::cout << "Effibility score: " << feaScore << "\n";
}

void PairSet::updateCovPairs()
{
	uint i, j, k;
	for( i = 0; i < setLength( true ); i++ )
		for( j = 0; j < setLength( false ); j++ )
		{
			for( k = 0; k < ref->length(); k++ )
			{
				if( mpFwd[i][k] > 0 && mpRev[j][k] > 0 )
					covPairs[i*setLength( false ) + j][k] = mpRev[j][k] - mpFwd[i][k];
				else
					covPairs[i*setLength( false ) + j][k] = -1;
				// std::cout << covPairs[i*revLength() + j][k] << " ";
			}
			// std::cout << covPairs[i*revLength() + j][2105] << " ";
		}
		// std::cout << "\n";
}

Stats PairSet::statsAmpLen( double qtile )
{
	uint i,j;
	vector<uint> ampLen;
	ampLen.reserve( covPairs.size() * ref->length() ); //no larger than that
	for( i = 0; i < covPairs.size(); i++ )
		for( j = 0; j < ref->length(); j++ )
			if( covPairs[i][j] > 0 )
				ampLen.push_back( covPairs[i][j] );
	
	Stats stats;
	if( ampLen.size() > 0 )
	{
		std::sort( ampLen.begin(), ampLen.end() );
		// std::cout << "amp len dist" << ampLen[0] << " " << ampLen[(ampLen.size())/2] << " " << ampLen[ampLen.size()-1] << "\n";
		stats.median = ampLen[(ampLen.size())/2];
		stats.quantile = ampLen[(int)(ampLen.size() * qtile)];
		return( stats );
	}
	// in case of zero coverage
	stats.median = -1;
	stats.quantile = -1;
	return( stats );
}

/*int PairSet::medianAmpLen()
{
	uint i,j;
	vector<int> ampLen;
	ampLen.reserve( covPairs.size() * ref->length() ); // it would be this long
	// copy all in a vector
	for( i = 0; i < covPairs.size(); i++ )
		ampLen.insert( ampLen.end(), covPairs[i].begin(), covPairs[i].end() );
	
	// sort and find the last -1, i.e NaN
	std::sort( ampLen.begin(), ampLen.end() );
	int noVal = -1;
	std::vector<int>::iterator it = std::find_end(ampLen.begin(), ampLen.end(), &noVal, &noVal);
	
	if( it == (ampLen.end() - 1) ) // all NaN
		return( -1 );
	else
		// median of the remainder
		return *( it + std::distance(it,ampLen.end())/2 );	
}*/


/*int PairSet::quantileAmpLen( double qtile )
{
	uint i,j;
	vector<uint> ampLen;
	for( i = 0; i < covPairs.size(); i++ )
		for( j = 0; j < ref->length(); j++ )
			if( covPairs[i][j] > 0 )
				ampLen.push_back( covPairs[i][j] );
	if( ampLen.size() > 0 )
	{
		std::sort( ampLen.begin(), ampLen.end() );
		return( ampLen[(int)(ampLen.size() * qtile)] );
	}
	// in case of zero coverage
	return( -1 );
}*/

/* Public methods */

uint PairSet::setLength( bool fwdSet )
{
	if( fwdSet )
		return seqan::length(fwd);
	return seqan::length(rev);
}

uint PairSet::prLength( bool fwdSet, uint p )
{
	if( fwdSet )
		return seqan::length(fwd[p]);
	return seqan::length(rev[p]);
}

IupacString PairSet::getPrimer( bool fwdSet, uint p )
{
	if( fwdSet )
		return(fwd[p]);
	return(rev[p]);
}

void PairSet::printSet()
{
	uint count = 0;
	//uint oldw = std::cout.width();
	//std::cout.width(22); 
	std::cout << "\nForward primers  \tReverse primers\n";
	while( count < setLength(true) && count < setLength(false) )
	{
		std::cout << fwd[count] << "\t" << rev[count] << "\n";
		count++;
	}
	// any of the two is over
	while( count < setLength(true)  )
	{
		std::cout << fwd[count] << "\t                 " << "\n";
		count++;
	}
	// any of the two is over
	while( count < setLength(false)  )
	{
		std::cout << "                 \t" << rev[count] << "\n";
		count++;
	}
}

double PairSet::getCovSc() const
{
	return( covScore );
}

double PairSet::getVarSc() const
{
	return( varScore );
}

double PairSet::getEffSc() const
{
	return( feaScore );
}

// getter for the scaled weighted score, to be minimized
double PairSet::getWSc( vector<double>& alpha, Bounds& b )
{
	return( 
		alpha[0] * (b.maxCovSc - covScore ) / (b.maxCovSc - b.minCovSc ) + 
		alpha[1] * (b.maxFeaSc - feaScore ) / (b.maxFeaSc - b.minFeaSc ) +
		alpha[2] * (varScore - b.minVarSc ) / (b.maxVarSc - b.minVarSc ));		
}

/* flip the n-th nucleotide of the p-th primer from either 
fwd or rev to its l-th subsequent letter in lexicographic order 
and recompute scores */
void PairSet::flip( bool fwdSet, uint n, uint p, uint l )
{
	char bases[] = {'A','C','G','T'};
	uint orig,j;
	
	if( fwdSet )
	{
		orig = log2(seqan::ordValue(fwd[p][n])); // 0..3 value
		fwd[p][n] = bases[(orig + l)%4];
		// std::cout << fwd[p] << "\n";
		// update matching positions and melting temperature
		mpFwd[p] = ref->matchPosFwd( fwd[p] );
		tmFwd[p] = meltTemp( fwd[p] );
	}
	else
	{
		orig = log2(seqan::ordValue(rev[p][n])); // 0..3 value
		rev[p][n] = bases[(orig + l)%4];
		// update matching positions  and melting temperature
		// std::cout << rev[p] << "\n";
		mpRev[p] = ref->matchPosRev( rev[p] );
		tmRev[p] = meltTemp( rev[p] );
	}
	
	// update scores
	computeScores();
}

/* add the l-th nucleotide letter in lexicographic order to the 
p-th primer from either fwd or rev and to either its head or tail  
and recompute scores */
void PairSet::add( bool fwdSet, bool toHead, uint p, uint l )
{
	char bases[] = {'A','C','G','T'};
	if( fwdSet )
	{
		if( toHead )
			insertValue(fwd[p],0,bases[l]);
		else
			fwd[p] += bases[l];
		
		// update matching positions and melting temperature
		mpFwd[p] = ref->matchPosFwd( fwd[p] );
		tmFwd[p] = meltTemp( fwd[p] );	
	}
	else 
	{
		if( toHead )
			insertValue(rev[p],0,bases[l]);
		else
			rev[p] += bases[l];
		
		// update matching positions and melting temperature
		mpRev[p] = ref->matchPosFwd( rev[p] );
		tmRev[p] = meltTemp( rev[p] );	
	}
	
	// update scores
	computeScores();
}

/* trim either the head or the tail nucleotide from the p-th primer 
of either fwd or rev and recompute scores */
void PairSet::trim( bool fwdSet, bool fromHead, uint p )
{
	if( fwdSet )
	{
		if( fromHead )
			seqan::erase(fwd[p],0);
		else
			seqan::eraseBack(fwd[p]);
		
		// update matching positions and melting temperature
		mpFwd[p] = ref->matchPosFwd( fwd[p] );
		tmFwd[p] = meltTemp( fwd[p] );	
	}
	else 
	{
		if( fromHead )
			seqan::erase(rev[p],0);
		else
			seqan::eraseBack(rev[p]);
		
		// update matching positions and melting temperature
		mpRev[p] = ref->matchPosFwd( rev[p] );
		tmRev[p] = meltTemp( rev[p] );	
	}
	
	// update scores
	computeScores();	
}

