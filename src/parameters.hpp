#ifndef PARAMETERS_H
#define PARAMETERS_H

/* Parameters of the various pieces of the algorithm, with defaults */
class Parameters
{
public:
	
	Parameters( uint randSeed = 0, uint maxMismatches = 2, uint minPrimerLen = 17, 
		uint maxPrimerLen = 21, uint minMeltTemp = 52, double minGCCont = 0.5, 
		double maxGCCont = 0.7, uint maxDimers = 8, double maxDeltaTm = 3, uint maxHomopLen = 4,
		uint maxAmplLenSpanCov = 200, uint maxAmplLenSpanFeas = 100,
		double maxAlenSpanFeasQtile = 0.01, uint minMeltTempInterv = 2, 
		double minGCContInterv = 0.1, double maxDimersInterv = 3, 
		double deltaTmInterv = 2, uint maxHomopLenInterv = 2,
   	uint maxAmplLenSpanFeasInterv = 50 ) :
		seed(randSeed), minPLen(minPrimerLen), maxPLen(maxPrimerLen),
		minTm(minMeltTemp), minGC(minGCCont), maxGC(maxGCCont), maxDim(maxDimers),
		maxHomLen(maxHomopLen),dTm(maxDeltaTm),maxALenSpanCov(maxAmplLenSpanCov), 
		maxALenSpanFeas(maxAmplLenSpanFeas), maxALenSpanFeasQtile(maxAlenSpanFeasQtile),
		minTmInt(minMeltTempInterv), minGCInt(minGCContInterv), 
		maxDimInt(maxDimersInterv), maxHomLenInt(maxHomopLenInterv),  
		dTmInt(deltaTmInterv), maxALenSpanFeasInt(maxAmplLenSpanFeasInterv) 
	{}
		
	// seed of the random number generator
	uint seed;
	
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
	uint maxALenSpanFeas;
	double maxALenSpanFeasQtile;
	
	// fuzzy intervals
	uint minTmInt; 
	double minGCInt;
	double maxDimInt;
	uint maxHomLenInt; 
	double dTmInt;
	uint maxALenSpanFeasInt;	
};

	
#endif