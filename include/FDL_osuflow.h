/**
 * FDL_osuflow.h
 *
 *  Created on: Jan 26, 2010
 *  Author: Abon
 */

#ifndef FDL_OSUFLOW_H_
#define FDL_OSUFLOW_H_

#include "FDL_header.h"
#include "FDL_util.h"
#include "FDL_streamline.h"

class FDL_osuflow
{
public:

	OSUFlow *osuflow;                        /**< Pointer to OSUFlow object */
	int nSeed;								 /**< Total number of seeds generated */
	VECTOR3 *seeds;							 /**< List of all seeds generated */
	int *streamlineLengthList;				 /**< List of lengths of all streamlines generated */
	list<VECTOR3*> streamlineList;	         /**< List of streamlines generated */

	int dataDim[3];
	VECTOR3 dataLow;
	VECTOR3 dataHigh;

	FILE *seedFile;
	FILE *streamlineFile;

public:

	/**
	 * Constructor.
	 */
	FDL_osuflow( OSUFlow *osuflow );
	/**
	 * Seed Generation (Regular) Function.
	 */
	void generateSeedAtEachVoxel();
	/**
	 * Seed Generation (Regular within a subvolume) Function.
	 */
	void generateSeedAtEachVoxel( float* low, float* high );
	/**
	 * Seed Generation (along a line) Function.
	 */
	void generateSeedAlongLine( float* start, float* end, int nseed );
	/**
	 * Seed Generation (on an arbitrarily oriented plane) Function.
	 */
	void generateSeedOnPlane( VECTOR3 *corners );
	/**
	 * Seed Generation (in a arbitrarily oriented box) Function.
	 */
	void generateSeedInBox( VECTOR3 *corners );	
	/**
	 * Seed Generation (Random within a subvolume) Function.
	 */
	void generateSeedInBox( float* low, float* high, int nseed );
	/**
	 * Seed Generation (Random) Function.
	 */
	/**
	 * Seed Generation (Random) Function.
	 */
	void generateSeedsRandom( int nSeed );
	/**
	 * Streamline computation function.
	 */
	void advect( int nStep, float p1, float p2 );
	/**
	 * Streamline computation (forward direction only) function.
	 */
	void advect_unidirection( int nStep, float p1, float p2 );
	/**
	 * Streamline computation helper function.
	 */
	VECTOR3* joinTwoTraceDirections( VECTOR3 *fwdTraceData, int fwdTraceLength,
			   					     VECTOR3 *backTraceData, int backTraceLength );

	/**
	 * File I/O functions
	 */
	void readSeeds( int nSeed, const char *seedFileName );
	void writeSeeds( const char *seedFileName );
	void writeStreamlines( const char *streamlineFileName );

	/**
	 * Destructor
	 */
	~FDL_osuflow();
};

#endif /* FDL_OSUFLOW_H_ */
