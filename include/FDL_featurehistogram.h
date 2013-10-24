/*
 * FDL_featurehistogram.h
 *
 *  Created on: Nov 9, 2012
 *      Author: abon
 */

#ifndef FDL_FEATUREHISTOGRAM_H_
#define FDL_FEATUREHISTOGRAM_H_

#include "FDL_header.h"
#include "FDL_util.h"

class FDL_featurehistogram
{
public:

	int nData;
	int nResolution;
	int nBin;
	double binWidth;
	double* freqList;

public:

	// Constructor
	FDL_featurehistogram();
	FDL_featurehistogram( int nr, int nb );
	FDL_featurehistogram( const FDL_featurehistogram& that );
	FDL_featurehistogram& operator= ( const FDL_featurehistogram& that );

	// Destructor
	~FDL_featurehistogram();

	// Initialize
	void initialize( int nr, int nb );

	// Get-set functions
	int getNumResolutions();
	int getNumBin();
	void getFrequencyList( double* fList );

	// Construct and/or update histogram
	void update( list<SLSegment> segmentList, int res );
	void normalize();

	int index2DTo1D( int r, int b );

	void printFeatureHisotgram();
	void saveFeatureHisotgram( FILE* ptr );
};


#endif /* FDL_FEATUREHISTOGRAM_H_ */
