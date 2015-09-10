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

#ifndef PARAMETERS_H
#define PARAMETERS_H

/// parameters of the various pieces of the algorithm, with defaults
class Parameters
{
public:
	
	Parameters( uint randSeed = 0, uint restarts = 20, uint maxMismatches = 2, 
		uint minPrimerLen = 17, uint maxPrimerLen = 21, uint minMeltTemp = 52, 
		double minGCCont = 0.5, double maxGCCont = 0.7, uint maxDimers = 8, 
		double maxDeltaTm = 3, uint maxHomopLen = 4,
		uint maxAmplLenSpanCov = 200, uint maxAmplLenSpanEff = 50,
		double maxAlenSpanEffQtile = 0.01, uint minMeltTempInterv = 2, 
		double minGCContInterv = 0.1, double maxDimersInterv = 3, 
		double deltaTmInterv = 2, uint maxHomopLenInterv = 2,
   	uint maxAmplLenSpanEffInterv = 50, std::string outFileName = "out" ) :
		seed(randSeed), rest(restarts), minPLen(minPrimerLen), maxPLen(maxPrimerLen),
		minTm(minMeltTemp), minGC(minGCCont), maxGC(maxGCCont), maxDim(maxDimers),
		maxHomLen(maxHomopLen),dTm(maxDeltaTm),maxALenSpanCov(maxAmplLenSpanCov), 
		maxALenSpanEff(maxAmplLenSpanEff), maxALenSpanEffQtile(maxAlenSpanEffQtile),
		minTmInt(minMeltTempInterv), minGCInt(minGCContInterv), 
		maxDimInt(maxDimersInterv), maxHomLenInt(maxHomopLenInterv),  
		dTmInt(deltaTmInterv), maxALenSpanEffInt(maxAmplLenSpanEffInterv),
		outFName(outFileName)
	{}
		
	// seed of the random number generator
	uint seed;
	
	// number of restarts of the local search algorithm
	uint rest;
	
	// primer alignment parameters
	uint maxMism;
	
	// primer features
	uint minPLen;
	uint maxPLen;
	uint minTm;
	double minGC; 
	double maxGC; 
	uint maxDim;
	uint maxHomLen;
	
	// primer pair features 
	double dTm;	
	uint maxALenSpanCov;
	uint maxALenSpanEff;
	double maxALenSpanEffQtile;
	
	// fuzzy intervals
	uint minTmInt; 
	double minGCInt;
	double maxDimInt;
	uint maxHomLenInt; 
	double dTmInt;
	uint maxALenSpanEffInt;
	
	// output file name
	std::string outFName;
};

	
#endif