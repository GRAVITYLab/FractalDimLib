/*
 * FDL_featurehistogram.cpp

 *
 *  Created on: Nov 9, 2012
 *      Author: abon
 */

#include "FDL_featurehistogram.h"

FDL_featurehistogram::FDL_featurehistogram()
{
	nData = 0;
	nResolution = 0;
	nBin = 0;
	binWidth = 0;
	freqList = NULL;
}

FDL_featurehistogram::FDL_featurehistogram( int nr, int nb )
{
	nData = 0;
	nResolution = nr;
	nBin = nb;
	binWidth = 3.0 / nBin;
	freqList = NULL;

	// Allocate
	freqList = new double[nResolution*nBin];
	memset( freqList, 0, sizeof(double)*nResolution*nBin );
}

FDL_featurehistogram::FDL_featurehistogram( const FDL_featurehistogram& that )
{
	this->nData = 0;
	this->nResolution = that.nResolution;
	this->nBin = that.nBin;
	this->binWidth = that.binWidth;
	this->freqList = NULL;

	if( that.freqList != NULL )
		memcpy( this->freqList, that.freqList, sizeof(double)*this->nResolution*this->nBin );
}

FDL_featurehistogram& FDL_featurehistogram::operator= ( const FDL_featurehistogram& that )
{
	if (this != &that ) // protect against invalid self-assignment
	{
		this->nData = 0;
		this->nResolution = that.nResolution;
		this->nBin = that.nBin;
		this->freqList = NULL;

		if( that.freqList != NULL )
			memcpy( this->freqList, that.freqList, sizeof(double)*this->nResolution*this->nBin );
	}
	// by convention, always return *this
	return *this;
}

FDL_featurehistogram::~FDL_featurehistogram()
{
	if( freqList != NULL )	delete [] freqList;
}

void
FDL_featurehistogram::initialize( int nr, int nb )
{
	nData = 0;
	nResolution = nr;
	nBin = nb;
	binWidth = 3.0 / nBin;
	freqList = NULL;

	// Allocate
	freqList = new double[nResolution*nBin];
	memset( freqList, 0, sizeof(double)*nResolution*nBin );
}

int
FDL_featurehistogram::getNumResolutions()
{
	return nResolution;
}

int
FDL_featurehistogram::getNumBin()
{
	return nBin;
}

void
FDL_featurehistogram::getFrequencyList( double* fList )
{
	assert( fList != NULL );
	memcpy( fList, freqList, sizeof(double)*nResolution*nBin );
}

void
FDL_featurehistogram::update( list<SLSegment> segmentList, int res )
{
	int binId = 0, binIdOverall = 0;
	list<SLSegment>::const_iterator pIter;

	// Update histogram
	for( pIter = segmentList.begin();
		 pIter != segmentList.end();
		 pIter++ )
	{
		// Determine bin Id for current resolution pair
		binId = (int) floor( (*pIter).score / binWidth );
		if( binId<0 ) binId = 0;
		if( binId>=nBin ) binId = nBin;

		// Determine bin Id including all resolution pairs
		binIdOverall = index2DTo1D( res, binId );

		// Update frequency
		freqList[binIdOverall] ++;

	}// end for

	// Update number of data points
	nData += segmentList.size();
}

void
FDL_featurehistogram::normalize()
{
	assert( nData > 0 );

	for( int i=0; i<(nResolution*nBin); i++ )
		freqList[i] /= nData;
}

int
FDL_featurehistogram::index2DTo1D( int r, int b )
{
	return (r*nBin + b);
}

void
FDL_featurehistogram::printFeatureHisotgram()
{
}

void
FDL_featurehistogram::saveFeatureHisotgram( FILE* ptr )
{
	for( int i=0; i<(nResolution*nBin)-1; i++ )
		fprintf( ptr, "%g, ", freqList[i] );
	fprintf( ptr, "%g\n", freqList[(nResolution*nBin)-1] );
}


