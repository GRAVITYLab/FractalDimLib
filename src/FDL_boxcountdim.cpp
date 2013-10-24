/*
 * FDL_boxcountdim.cpp
 *
 *  Created on: Jan 31, 2011
 *      Author: abon
 */

#include "FDL_boxcountdim.h"

FDL_boxcountdim::FDL_boxcountdim( VECTOR3 *lowgrid, VECTOR3 *highgrid,
									  bool transformflag )
{
	// Set extents of the field
	lowGrid.Set( lowgrid->x(), lowgrid->y(), lowgrid->z() );
	highGrid.Set( highgrid->x(), highgrid->y(), highgrid->z() );

	// Set grid resolution pair to default
	setGridResolutionPair( 1, 2 );

	// Set geometry transformation related flag
	this->toTransformGeometry = transformflag;

	// NULL initialize rest of variables

}// end constructor

/**
 * boxLenLow: Lower box size
 * boxLenHigh: Higher box size
 */
void
FDL_boxcountdim::setGridResolutionPair( double boxLenLow, double boxLenHigh )
{
	// Set box length for two resolutions
	boxLen[0].Set( boxLenLow, boxLenLow, boxLenLow );
	boxLen[1].Set( boxLenHigh, boxLenHigh, boxLenHigh );
	
	// Compute grid properties for both resolutions
	for( int i=0; i<2;  i++ )
	{
		resolutions[i][0] = ceil( ( highGrid[0]-lowGrid[0] ) / boxLen[i].x() );
		resolutions[i][1] = ceil( ( highGrid[1]-lowGrid[1] ) / boxLen[i].y() );
		resolutions[i][2] = ceil( ( highGrid[2]-lowGrid[2] ) / boxLen[i].z() );

		this->nBox[i] = this->resolutions[i].x() * this->resolutions[i].y() * this->resolutions[i].z();

		#ifdef DEBUG_MODE
		fprintf( stderr, "Box length at resolution %d: %g %g %g\n", i, this->boxLen[i][0], this->boxLen[i][1], this->boxLen[i][2] );
		fprintf( stderr, "Number of boxes at resolution %d: %g %g %g\n", i, this->resolutions[i].x(), this->resolutions[i].y(), this->resolutions[i].z() );
		fprintf( stderr, "Total number of boxes at resolution %d: %d\n", i, this->nBox[i] );
		#endif
	}

}// end function

double
FDL_boxcountdim::computeBCD_SingleStreamline_SingleWindow( VECTOR3* trace, int traceLength )
{
	// Allocate temporary memory for storing the trace and the box count for each resolution
	VECTOR3 traceCopy[traceLength];
	double boxcountArray[2];

	// Transform geometry, if needed OR just create a copy of the data
	memcpy( traceCopy, trace, sizeof(VECTOR3) * traceLength );
	//fprintf( stderr, "Trace length: %d\n", traceLength );
	if( toTransformGeometry == true && traceLength >= 3 )
		transformGeometry( traceCopy, traceLength );

	#ifdef DEBUG_MODE
	for( int i=0; i<traceLength; i++ )
		fprintf( stderr, "%d: %g, %g, %g, %g, %g, %g\n", i, trace[i][0], trace[i][1], trace[i][2], traceCopy[i][0], traceCopy[i][1], traceCopy[i][2] );

	float mX = traceCopy[0][0];
	float MX = traceCopy[0][0];
	float mY = traceCopy[0][1];
	float MY = traceCopy[0][1];
	float mZ = traceCopy[0][2];
	float MZ = traceCopy[0][2];
	for( int i=0; i<traceLength; i++ )
	{
		if( traceCopy[i][0] < mX ) mX = traceCopy[i][0];
		if( traceCopy[i][0] > MX ) MX = traceCopy[i][0];

		if( traceCopy[i][1] < mY ) mY = traceCopy[i][1];
		if( traceCopy[i][1] > MY ) MY = traceCopy[i][1];

		if( traceCopy[i][2] < mZ ) mZ = traceCopy[i][2];
		if( traceCopy[i][2] > MZ ) MZ = traceCopy[i][2];
	}
	fprintf( stderr, "New Limits: %g, %g, %g, %g, %g, %g\n", mX, MX, mY, MY, mZ, MZ );
	#endif

	#ifdef DEBUG_MODE
	fprintf( stderr, "Trace segment:\n" );
	for(int i=0; i<traceLength; i++ )
		fprintf( stderr, "%d: %g %g %g\n", i, traceCopy[i][0], traceCopy[i][1], traceCopy[i][2] );
	#endif

	// Compute the box count for each resolution
	for( int iRes=0; iRes<2; iRes++ )
		boxcountArray[iRes] = computeBC_SingleResolution( traceCopy, traceLength, iRes );
		//boxcountArray[iRes] = computeBC_SingleResolution2( traceCopy, traceLength, iRes );

	#ifdef DEBUG_MODE
	for( int i=0; i<2; i++ )
		fprintf( stderr, "Box count for %d-th resolution %g\n", i, boxcountArray[i] );
	#endif

	// compute and return final value of box count dimension
	double countRatio = boxcountArray[0] / boxcountArray[1];
	double fracDim = log( countRatio ) / log ( 2.0f );

	#ifdef DEBUG_MODE
	fprintf( stderr, "Computed FD of trace segment: %g %g %g\n", boxcountArray[0], boxcountArray[1], fracDim );
	#endif

	// Return computed fractal dimension
	return fracDim;
}// end function

double
FDL_boxcountdim::computeBCD_SingleStreamline_Full( FDL_streamline* streamline )
{
	// Get streamline trace
	int traceLen =  streamline->getLength();
	VECTOR3* totalTrace = new VECTOR3[traceLen];
	streamline->getTotalTrace( totalTrace );

	// Compute
	double score = computeBCD_SingleStreamline_SingleWindow( totalTrace, traceLen );

	// Clear
	delete [] totalTrace;

	return score;
}

double
FDL_boxcountdim::computeBCD_SingleStreamline_InAxisAlignedBox( VECTOR3* trace, int traceLength,
																	   double *boxMin, double *boxMax )
{
	set<VECTOR3,Comparator2> *boxidset = NULL;
	//fprintf( stderr, "%g %g %g %g %g %g\n", boxMin[0], boxMin[1], boxMin[2], boxMax[0], boxMax[1], boxMax[2]  );

	// Allocate memory for storing box count for each resolution
	double boxcountArray[2];

	for( int iRes=0; iRes<2; iRes++ )
	{
		// Initialize set for current resolution
		boxidset = new set<VECTOR3,Comparator2>();

		// Add up box count for each contributing window
		int validTraceCount = 0;
		for( int iT=0; iT<traceLength-1; iT++ )
		{
			// The next step of the streamline
			VECTOR3* nextTrace = (trace+iT);
			//fprintf( stderr, "%g %g %g\n", trace[iT].x(), trace[iT].y(), trace[iT].z() );
			//fprintf( stderr, "%g %g %g\n", trace[iT+1].x(), trace[iT+1].y(), trace[iT+1].z() );
			
			// Skip trace step if any of the two ends of the trace are outside the box			
			if( trace[iT].x() < boxMin[0] || trace[iT].x() > boxMax[0] ||
				trace[iT].y() < boxMin[1] || trace[iT].y() > boxMax[1] ||
				trace[iT].z() < boxMin[2] || trace[iT].z() > boxMax[2] ||
				trace[iT+1].x() < boxMin[0] || trace[iT+1].x() > boxMax[0] ||
				trace[iT+1].y() < boxMin[1] || trace[iT+1].y() > boxMax[1] ||
				trace[iT+1].z() < boxMin[2] || trace[iT+1].z() > boxMax[2] )
				continue;
			else
				validTraceCount ++;
			
			// Compute box count ONLY if both end points of the trace are inside the box
			computeBC_SingleResolution2( nextTrace, 2, iRes, boxidset );
		}

		// Handle special cases		
		if( validTraceCount == 0 )return 0;
		if( validTraceCount < 10 )return 1;

		// Get total box count for this resolution
		boxcountArray[iRes] = boxidset->size();

		// Clear set
		delete boxidset;

	}// end for
	
	#ifdef DEBUG_MODE
	for( int i=0; i<this->nResolution; i++ )
		fprintf( stderr, "Box count for %d-th resolution %g\n", i, boxcountArray[i] );
	#endif
	
	// Return 0 if the streamline does not interact with the box at all
	if( boxcountArray[0] == 0 || boxcountArray[1] == 0 )
	{
		return 0.0f;
	}

	// compute and return final value of box count dimension
	double countRatio = boxcountArray[0] / boxcountArray[1];
	double fracDim = log( countRatio ) / log ( 2.0f );
	
	#ifdef DEBUG_MODE
	fprintf( stderr, "Computed FD of the streamline within current box segment bunch: %g %g %g\n", boxcountArray[0], boxcountArray[1], fracDim );
	#endif

	// Return computed fractal dimension
	return fracDim;	
}

void
FDL_boxcountdim::computeBCD_SingleStreamline_windowed( FDL_streamline* streamline, int windowLength,
															  double** bcdArray 	)
{
	// Get the trace and trace length of the streamline
	int totalTraceLength = streamline->getLength();
	VECTOR3* totalTrace = new VECTOR3[totalTraceLength];
	streamline->getTotalTrace( totalTrace );

	// Allocate memory for storing box count for each segment
	(*bcdArray) = new double[totalTraceLength];

	// Slide window until end of streamline reached
	int startIndex = 0;
	VECTOR3* localTrace;
	while( startIndex < totalTraceLength )
	{
		// Determine the actual window length
		if( startIndex + windowLength > totalTraceLength )
			windowLength = totalTraceLength - startIndex;

		#ifdef DEBUG_MODE
		fprintf( stderr, "Moving to next window of length %d from trace point %d\n", windowLength, startIndex );
		#endif

		// Compute and store box count dimension for this window
		localTrace = totalTrace + startIndex;
		(*bcdArray)[startIndex] = computeBCD_SingleStreamline_SingleWindow( totalTrace, windowLength );

		// Move one step ahead along streamline
		startIndex ++;

	}// end while

	delete [] totalTrace;

}// end function

void
FDL_boxcountdim::computeBCD_SingleStreamline_DisjointSegments( FDL_streamline* streamline, int segmentLength, double *bcdArray )
{
	// Get the trace and trace length of the streamline
	int totalTraceLength = streamline->getLength();
	VECTOR3* totalTrace = new VECTOR3[totalTraceLength];
	streamline->getTotalTrace( totalTrace );

	// Compute number of windows possible
	int nSegment = (int) ceil( totalTraceLength / (float)segmentLength );
	#ifdef DEBUG_MODE	
	fprintf( stderr, "Stremline Length: %d\n", streamline->totalTraceLength );
	fprintf( stderr, "Number of segments: %d\n", nSegment );
	#endif
	
	// Clear up the fixed-length segment list of the current streamline
	streamline->fixedLengthSegmentList.clear();
	
	// Slide window until end of streamline reached
	int startIndex = 0;
	int bcdArrayIndex = 0;
	VECTOR3* localTrace;
	while( startIndex < totalTraceLength )
	{
		// Determine the actual window length
		if( startIndex + segmentLength > totalTraceLength )
			segmentLength = totalTraceLength - startIndex;

		#ifdef DEBUG_MODE
		fprintf( stderr, "Moving to next window of length %d from trace point %d\n", segmentLength, startIndex );
		#endif

		// Compute and store box count dimension for this segment
		localTrace = totalTrace + startIndex;
		bcdArray[bcdArrayIndex] = computeBCD_SingleStreamline_SingleWindow( localTrace, segmentLength );
		if( bcdArray[bcdArrayIndex] <0 || bcdArray[bcdArrayIndex] > 3 ) 
		{
			fprintf( stderr, "%d %d %d %g\n", totalTraceLength, nSegment, startIndex, bcdArray[bcdArrayIndex] );
			exit(0);
		}
		
		// Update the fixed-length segment list of this streamline		
		struct Segment newsegment = { 0, startIndex, segmentLength, bcdArray[bcdArrayIndex] }; 
		streamline->fixedLengthSegmentList.push_back( newsegment );

		// Move to the start point of next segment
		startIndex += segmentLength;
		//startIndex += (segmentLength-1);

		// Increment outpt array index
		bcdArrayIndex ++;
		

	}// end while

}// end function

void
FDL_boxcountdim::computeBCD_MultiStreamline_Full( FDL_streamline* streamlineList, int nStreamline, double *bcdArray )
{
	assert( bcdArray != NULL );

	for( int iS=0; iS<nStreamline; iS++ )
	{
		bcdArray[iS] = computeBCD_SingleStreamline_Full( streamlineList + iS );
		#ifdef DEBUG_MODE
		fprintf( stderr, "BCR of %d-th streamline: %g\n", iS, bcdArray[iS] );
		#endif
	}

}// end function

void
FDL_boxcountdim::computeBCD_MultiStreamline_DisjointSegments( FDL_streamline* streamlineList,
																	  int nStreamline, int segmentLength,
																	  double** segmentedFdArray )
{
	int nSegmentNextStreamline = 0;

	for( int iS=0; iS<nStreamline; iS++ )
	{
		// Compute segment wise FD of the next streamline
		computeBCD_SingleStreamline_DisjointSegments( streamlineList + iS, segmentLength, segmentedFdArray[iS] );
		//fprintf( stderr, "Segmentwise FD of %d-th streamline computed.\n", iS );
	}
	
}// end function

double
FDL_boxcountdim::computeBCD_MultiStreamline_SingleWindow( FDL_streamline* streamlineList,
																 int nStreamline,
																 int *windowIdList, int* windowLengthList )
{
	set<VECTOR3,Comparator2> *boxidset = NULL;

	// Allocate memory for storing box count for each resolution
	double boxcountArray[2];

	for( int iRes=0; iRes<2; iRes++ )
	{
		// Initialize set for current resolution
		boxidset = new set<VECTOR3,Comparator2>();

		// Add up box count for each contributing window
		for( int iS=0; iS<nStreamline; iS++ )
		{
			// Get the trace and trace length of the next streamline
			int totalTraceLength = streamlineList[iS].getLength();
			VECTOR3* nextTrace = new VECTOR3[totalTraceLength];
			streamlineList[iS].getTotalTrace( nextTrace );

			computeBC_SingleResolution2( &nextTrace[windowIdList[iS]], windowLengthList[iS], iRes, boxidset );

			delete [] nextTrace;
		}

		// Get total box count for this resolution
		boxcountArray[iRes] = boxidset->size();

		// Clear set
		delete boxidset;

	}// end for

	#ifdef DEBUG_MODE
	for( int i=0; i<this->nResolution; i++ )
		fprintf( stderr, "Box count for %d-th resolution %g\n", i, boxcountArray[i] );
	#endif

	// compute and return final value of box count dimension
	double countRatio = boxcountArray[0] / boxcountArray[1];
	double fracDim = log( countRatio ) / log ( 2.0f );

	#ifdef DEBUG_MODE
	fprintf( stderr, "Computed FD of trace segment bunch: %g %g %g\n", boxcountArray[0], boxcountArray[1], fracDim );
	#endif

	// Return computed fractal dimension
	return fracDim;

}// end function


double
FDL_boxcountdim::computeBC_SingleResolution( VECTOR3* trace, int traceLength, int iRes )
{
	int boxCount = 0;
	set<VECTOR3,Comparator2> boxidset;// = new set<VECTOR3,Comparator2>();

	// Data structute used to prevent counting the same box multiple times
	//bool *countcheckArray = new bool[this->nBox[iRes]];
	//for( int i=0; i<this->nBox[iRes]; i++ )
	//	countcheckArray[i] = false;

	// Resolve sigularity (Length 1 streamlines)
	if( traceLength == 1 )
		return 1.0f;

	// Traverse the streamline and record the boxes corresponding to each position
	for( int i=0; i<traceLength-1; i++ )
	{
		#ifdef DEBUG_MODE
		fprintf( stderr, "%d th point on trace: %g %g %g\n", i, trace[i].x(), trace[i].y(), trace[i].z() );
		fprintf( stderr, "%d th point on trace: %g %g %g\n", i+1, trace[i+1].x(), trace[i+1].y(), trace[i+1].z() );
		#endif

		// Run 3D DDA based algo to get box intersection count and check duplicates using array
		//boxCount += this->countBoxIntersectAlongRay_DDA_Array( *(trace+i), *(trace+i+1), iRes, countcheckArray );

		this->countBoxIntersectAlongRay_DDA_Set( *(trace+i), *(trace+i+1), iRes, &boxidset );

	}// end for

	//delete [] countcheckArray;
	//if( boxidset != NULL )
	//{
	boxCount = (int)boxidset.size();//->size();

	boxidset.clear();//->clear();

		//delete boxidset;
	//}
	
	return (double)boxCount;

}// end function

void
FDL_boxcountdim::computeBC_SingleResolution2( VECTOR3* trace, int traceLength,
													int iRes,
													set<VECTOR3,Comparator2> *boxidset )
{
	assert( boxidset != NULL );

	// Data structute used to prevent counting the same box multiple times
	//bool *countcheckArray = new bool[this->nBox[iRes]];
	//for( int i=0; i<this->nBox[iRes]; i++ )
	//	countcheckArray[i] = false;

	// Resolve sigularity ( Length 1 traces )
	if( traceLength == 1 )
	{
		VECTOR3 singlePointBoxid = this->positionToBox3D( *trace, iRes );

		// Insert point to box count set, if unique
		boxidset->insert( singlePointBoxid );

		return;
	}

	// Traverse the streamline and record the boxes corresponding to each position
	for( int i=0; i<traceLength-1; i++ )
	{
		#ifdef DEBUG_MODE
		fprintf( stderr, "%d th point on trace: %g %g %g\n", i, trace[i].x(), trace[i].y(), trace[i].z() );
		fprintf( stderr, "%d th point on trace: %g %g %g\n", i+1, trace[i+1].x(), trace[i+1].y(), trace[i+1].z() );
		#endif

		// Run 3D DDA based algo to get box intersection count and check duplicates using array
		//boxCount += this->countBoxIntersectAlongRay_DDA_Array( *(trace+i), *(trace+i+1), iRes, countcheckArray );
		this->countBoxIntersectAlongRay_DDA_Set( *(trace+i), *(trace+i+1), iRes, boxidset );

	}// end for

}// end function

int
FDL_boxcountdim::positionToBox( VECTOR3 position, int iRes )
{
	float xid = floor( ( position[0] - this->lowGrid[0] ) / this->boxLen[iRes][0] );
	float yid = floor( ( position[1] - this->lowGrid[1] ) / this->boxLen[iRes][1] );
	float zid = floor( ( position[2] - this->lowGrid[2] ) / this->boxLen[iRes][2] );

	#ifdef DEBUG_MODE
	fprintf( stderr, "Trace position: %g %g %g\n", position[0], position[1], position[2] );
	fprintf( stderr, "Grid lower bound: %g %g %g\n", lowGrid[0], lowGrid[1], lowGrid[2] );
	fprintf( stderr, "Box Length: %g %g %g\n", this->boxLen[iRes][0], this->boxLen[iRes][1], this->boxLen[iRes][2] );
	fprintf( stderr, "Box id: %g %g %g\n", xid, yid, zid );
	#endif

	return (int) ( zid * this->resolutions[iRes].x() * this->resolutions[iRes].y() + yid * this->resolutions[iRes].x() + xid );

}// end function

VECTOR3
FDL_boxcountdim::positionToBox3D( VECTOR3 position, int iRes )
{
	VECTOR3 boxid( floor( ( position[0] - this->lowGrid[0] ) / this->boxLen[iRes][0] ),
				   floor( ( position[1] - this->lowGrid[1] ) / this->boxLen[iRes][1] ),
				   floor( ( position[2] - this->lowGrid[2] ) / this->boxLen[iRes][2] ) );

	//boxid[0] = FDL_util<float>::clamp( boxid[0], 0, this->resolutions[iRes].x()-1.0 );
	//boxid[1] = FDL_util<float>::clamp( boxid[1], 0, this->resolutions[iRes].y()-1.0 );
	//boxid[2] = FDL_util<float>::clamp( boxid[2], 0, this->resolutions[iRes].z()-1.0 );

	// box id can be negative
	//for( int i=0; i<3; i++ )
	//	if( boxid[i] < 0 )	boxid[i] = 0;


	#ifdef DEBUG_MODE
	fprintf( stderr, "Trace position: %g %g %g\n", position[0], position[1], position[2] );
	fprintf( stderr, "Grid lower bound: %g %g %g\n", lowGrid[0], lowGrid[1], lowGrid[2] );
	fprintf( stderr, "Box Length: %g %g %g\n", this->boxLen[iRes][0], this->boxLen[iRes][1], this->boxLen[iRes][2] );
	fprintf( stderr, "Box id: %g %g %g\n", boxid[0], boxid[1], boxid[2] );
	#endif

	return boxid;

}// end function

int
FDL_boxcountdim::boxId3DTo1D( VECTOR3 id, int iRes )
{
	return (int)( id[2] * this->resolutions[iRes].x() * this->resolutions[iRes].y() + id[1] * this->resolutions[iRes].x() + id[0] );

}// end function

int
FDL_boxcountdim::computeNextWindowSize_VoxelBased( FDL_streamline* streamline,
														  int startIndex, int voxelLimit )
{
	int nVoxelsCovered = 0;
	int nVoxelsCoveredThisTrace = 0;
	int traceBasedWindowLength = 0;
	int iL = 0;

	// Get the trace and trace length of the streamline
	int totalTraceLength = streamline->getLength();
	VECTOR3* totalTrace = new VECTOR3[totalTraceLength];
	streamline->getTotalTrace( totalTrace );

	while(1)
	{
		// Check whether the other end of the window has reached the end of streamline
		if( startIndex + traceBasedWindowLength > totalTraceLength )
		{
			traceBasedWindowLength = totalTraceLength - startIndex;
			break;
		}

		nVoxelsCoveredThisTrace = this->countBoxIntersectAlongRay_MD( *(totalTrace + startIndex + iL),
																      *(totalTrace + startIndex + iL + 1),
																      0 );

		if( nVoxelsCovered + nVoxelsCoveredThisTrace + 1 > voxelLimit )
			break;
		else
		{
			nVoxelsCovered += nVoxelsCoveredThisTrace;
			traceBasedWindowLength ++;
			iL++;
		}
	}

	fprintf( stderr, "Window Length: %d Voxel Coverage: %d\n", traceBasedWindowLength, nVoxelsCovered );

	delete [] totalTrace;

	return traceBasedWindowLength;

}// end function

int
FDL_boxcountdim::computeNextWindowSize_ArclengthBased( FDL_streamline* streamline,
															  int startIndex, double arclengthLimit )
{
	float lengthTraversed = 0;
	float lengthTraversedThisTrace = 0;
	int traceBasedWindowLength = 0;
	int iL = 0;

	// Get the trace and trace length of the streamline
	int totalTraceLength = streamline->getLength();
	VECTOR3* totalTrace = new VECTOR3[totalTraceLength];
	streamline->getTotalTrace( totalTrace );

	while(1)
	{
		// Check whether the other end of the window has reached the end of streamline
		if( startIndex + traceBasedWindowLength > totalTraceLength )
		{
			traceBasedWindowLength = totalTraceLength - startIndex;
			break;
		}

		lengthTraversedThisTrace = this->computeLengthAlongRay( *(totalTrace + startIndex + iL),
															   *(totalTrace + startIndex + iL + 1) );

		if( lengthTraversed + lengthTraversedThisTrace >= arclengthLimit )
			break;
		else
		{
			lengthTraversed += lengthTraversedThisTrace;
			traceBasedWindowLength ++;
			iL++;
		}
	}

	fprintf( stderr, "Window Length: %d Arclength Coverage: %g\n", traceBasedWindowLength, lengthTraversed );

	delete [] totalTrace;

	return traceBasedWindowLength;

return 0;

}// end function

double
FDL_boxcountdim::computeLengthAlongRay( VECTOR3 startpoint, VECTOR3 endpoint )
{
	return sqrt( ( endpoint[0]-startpoint[0] )*( endpoint[0]-startpoint[0] ) +
					  ( endpoint[1]-startpoint[1] )*( endpoint[1]-startpoint[1] ) +
					  ( endpoint[2]-startpoint[2] )*( endpoint[2]-startpoint[2] ) );
}// end function

int
FDL_boxcountdim::countBoxIntersectAlongRay_MD( VECTOR3 startPoint, VECTOR3 endPoint, int iRes )
{
	// Get the box ids for the start and the end point
	VECTOR3 startPointBoxid = this->positionToBox3D( startPoint, iRes );
	VECTOR3 endPointBoxid = this->positionToBox3D( endPoint, iRes );

	return  (int) ( fabs(startPointBoxid[0]-endPointBoxid[0] ) +
					fabs(startPointBoxid[1]-endPointBoxid[1] ) +
					fabs(startPointBoxid[2]-endPointBoxid[2] ) );

}// end function

void
FDL_boxcountdim::countBoxIntersectAlongRay_DDA_Set( VECTOR3 startPoint, VECTOR3 endPoint,
															int iRes, set<VECTOR3, Comparator2> *boxidset )
{
	// Get the box ids for the start and the edn point
	VECTOR3 startPointBoxid = this->positionToBox3D( startPoint, iRes );
	VECTOR3 endPointBoxid = this->positionToBox3D( endPoint, iRes );

	// Insert points to box count set, if unique
	if( boxidset != NULL )
	{
		boxidset->insert( startPointBoxid );
		boxidset->insert( endPointBoxid );
	}

	#ifdef DEBUG_MODE
	fprintf( stderr, "Start point: %g %g %g\n", startPoint[0], startPoint[1], startPoint[2] );
	fprintf( stderr, "Start point box id: %g %g %g\n", startPointBoxid[0], startPointBoxid[1], startPointBoxid[2] );
	fprintf( stderr, "End point: %g %g %g\n", endPoint[0], endPoint[1], endPoint[2] );
	fprintf( stderr, "End point box id at current iteration: %g %g %g\n", endPointBoxid[0], endPointBoxid[1], endPointBoxid[2] );
	#endif

	// Return if both are in the same box
	if( startPointBoxid == endPointBoxid )	return;// 1;

	// Compute the direction vector between start and end point
	VECTOR3 dir( endPoint[0] - startPoint[0],
				 endPoint[1] - startPoint[1],
				 endPoint[2] - startPoint[2] );
	dir.Normalize();

	// Compute initial value of step
	VECTOR3 step;
	if( endPoint[0] >= startPoint[0] ) step[0] = 1;
	else							   step[0] = -1;
	if( endPoint[1] >= startPoint[1] ) step[1] = 1;
	else							   step[1] = -1;
	if( endPoint[2] >= startPoint[2] ) step[2] = 1;
	else							   step[2] = -1;

	// Compute t delta: t for moving one box length
	VECTOR3 tDelta( this->boxLen[iRes][0] / fabs( dir[0] ),
					this->boxLen[iRes][1] / fabs( dir[1] ),
					this->boxLen[iRes][2] / fabs( dir[2] ) );

	// Compute initial value of max vector
	float tmaxX, tmaxY, tmaxZ;
	if( endPoint[0] >= startPoint[0] )
		tmaxX = ( this->boxLen[iRes][0]*(startPointBoxid[0]+1) - startPoint[0] ) / dir[0];
	else
		tmaxX = ( this->boxLen[iRes][0]*startPointBoxid[0] - startPoint[0] ) / dir[0];
	if( endPoint[1] >= startPoint[1] )
		tmaxY = ( this->boxLen[iRes][1]*(startPointBoxid[1]+1) - startPoint[1] ) / dir[1];
	else
		tmaxY = ( this->boxLen[iRes][1]*startPointBoxid[1] - startPoint[1] ) / dir[1];
	if( endPoint[2] >= startPoint[2] )
		tmaxZ = ( this->boxLen[iRes][2]*(startPointBoxid[2]+1) - startPoint[2] ) / dir[2];
	else
		tmaxZ = ( this->boxLen[iRes][2]*startPointBoxid[2] - startPoint[2] ) / dir[2];
	VECTOR3 tMax( tmaxX, tmaxY, tmaxZ );

	// Variables for iteration
	VECTOR3 oldPoint = startPoint;
	VECTOR3 oldPointBoxid = startPointBoxid;
	VECTOR3 nextPoint;
	VECTOR3 nextPointBoxid;

	#ifdef DEBUG_MODE
	fprintf( stderr, "Box Length: %g %g %g\n", this->boxLen[iRes][0], this->boxLen[iRes][1], this->boxLen[iRes][2] );
	fprintf( stderr, "dir: %g %g %g\n", dir[0], dir[1], dir[2] );
	fprintf( stderr, "tDel: %g %g %g\n", tDelta[0], tDelta[1], tDelta[2] );
	fprintf( stderr, "tMax: %g %g %g\n", tMax[0], tMax[1], tMax[2] );
	#endif

	while( !( oldPointBoxid == endPointBoxid ) &&
			oldPointBoxid.x() >= 0 &&
			oldPointBoxid.y() >= 0 &&
			oldPointBoxid.z() >= 0 )
	{
		//fprintf( stderr, "tMax: %g %g %g\n", tMax[0], tMax[1], tMax[2] );
		if( tMax[0] < tMax[1] )
		{
			if( tMax[0] < tMax[2] )
			{
				nextPointBoxid.Set( oldPointBoxid[0] + step[0], oldPointBoxid[1], oldPointBoxid[2] );
				tMax.add( tDelta[0], 0, 0 );
			}
			else
			{
				nextPointBoxid.Set( oldPointBoxid[0], oldPointBoxid[1], oldPointBoxid[2] + step[2] );
				tMax.add( 0, 0, tDelta[2] );
			}
		}
		else
		{
			if( tMax[1] < tMax[2] )
			{
				nextPointBoxid.Set( oldPointBoxid[0], oldPointBoxid[1] + step[1], oldPointBoxid[2] );
				tMax.add( 0, tDelta[1], 0 );
			}
			else
			{
				nextPointBoxid.Set( oldPointBoxid[0], oldPointBoxid[1], oldPointBoxid[2] + step[2] );
				tMax.add( 0, 0, tDelta[2] );
			}
		}// end if-else

		#ifdef DEBUG_MODE
		fprintf( stderr, "Next point box id at current iteration: %g %g %g\n", nextPointBoxid[0], nextPointBoxid[1], nextPointBoxid[2] );
		#endif

		// Check if next point is within data grid
		//if( nextPointBoxid[0] < 0 || nextPointBoxid[0] >= this->resolutions[iRes].x()
		// || nextPointBoxid[1] < 0 || nextPointBoxid[1] >= this->resolutions[iRes].y()
		 //|| nextPointBoxid[2] < 0 || nextPointBoxid[2] >= this->resolutions[iRes].z() )
		//{
			//fprintf( stderr, "Start point: %g %g %g\n", startPoint[0], startPoint[1], startPoint[2] );
			//fprintf( stderr, "End point: %g %g %g\n", endPoint[0], endPoint[1], endPoint[2] );
			//fprintf( stderr, "Start point box id: %g %g %g\n", startPointBoxid[0], startPointBoxid[1], startPointBoxid[2] );
			//fprintf( stderr, "End point box id at current iteration: %g %g %g\n", endPointBoxid[0], endPointBoxid[1], endPointBoxid[2] );
		//	exit(0);
	//	}

		// Insert points to box count set, if unique
		if( boxidset != NULL ) boxidset->insert( nextPointBoxid );

		// Update information for next iteration
		//oldPoint = nextPoint;
		oldPointBoxid = nextPointBoxid;

	}// end while

	//return (int)boxidset->size();

}// end function

int
FDL_boxcountdim::countBoxIntersectAlongRay_DDA_Array( VECTOR3 startPoint, VECTOR3 endPoint, int iRes, bool *countCheckArray )
{
	int boxCount = 0;

	// Get the box ids for the start and the edn point
	VECTOR3 startPointBoxid = this->positionToBox3D( startPoint, iRes );
	VECTOR3 endPointBoxid = this->positionToBox3D( endPoint, iRes );

	// Maintain a 1-D array to check if boxed are unique
	if( countCheckArray[this->boxId3DTo1D( startPointBoxid, iRes )] == false )
	{
		countCheckArray[this->boxId3DTo1D( startPointBoxid, iRes )] = true;
		boxCount ++;
	}
	if( countCheckArray[this->boxId3DTo1D( endPointBoxid, iRes )] == false )
	{
		countCheckArray[this->boxId3DTo1D( endPointBoxid, iRes )] = true;
		boxCount ++;
	}

	#ifdef DEBUG_MODE
	fprintf( stderr, "Start point box id: %g %g %g\n", startPointBoxid[0], startPointBoxid[1], startPointBoxid[2] );
	fprintf( stderr, "End point box id at current iteration: %g %g %g\n", endPointBoxid[0], endPointBoxid[1], endPointBoxid[2] );
	#endif

	// Return if both are in the same box
	if( startPointBoxid == endPointBoxid )	return 1;

	// Compute the direction vector between start and end point
	VECTOR3 dir( endPoint[0] - startPoint[0],
				 endPoint[1] - startPoint[1],
				 endPoint[2] - startPoint[2] );
	dir.Normalize();

	// Compute initial value of step
	VECTOR3 step;
	if( endPoint[0] >= startPoint[0] ) step[0] = 1;
	else							   step[0] = -1;
	if( endPoint[1] >= startPoint[1] ) step[1] = 1;
	else							   step[1] = -1;
	if( endPoint[2] >= startPoint[2] ) step[2] = 1;
	else							   step[2] = -1;

	// Compute t delta: t for moving one box length
	VECTOR3 tDelta( this->boxLen[iRes][0] / fabs( dir[0] ),
					this->boxLen[iRes][1] / fabs( dir[1] ),
					this->boxLen[iRes][2] / fabs( dir[2] ) );

	// Compute initial value of max vector
	float tmaxX, tmaxY, tmaxZ;
	if( endPoint[0] >= startPoint[0] ) tmaxX = ( this->boxLen[iRes][0]*(startPointBoxid[0]+1) - startPoint[0] ) / dir[0];
	else							   tmaxX = ( this->boxLen[iRes][0]*startPointBoxid[0] - startPoint[0] ) / dir[0];
	if( endPoint[1] >= startPoint[1] ) tmaxY = ( this->boxLen[iRes][1]*(startPointBoxid[1]+1) - startPoint[1] ) / dir[1];
	else							   tmaxY = ( this->boxLen[iRes][1]*startPointBoxid[1] - startPoint[1] ) / dir[1];
	if( endPoint[2] >= startPoint[2] ) tmaxZ = ( this->boxLen[iRes][2]*(startPointBoxid[2]+1) - startPoint[2] ) / dir[2];
	else							   tmaxZ = ( this->boxLen[iRes][2]*startPointBoxid[2] - startPoint[2] ) / dir[2];
	VECTOR3 tMax( tmaxX, tmaxY, tmaxZ );

	// Variables for iteration
	VECTOR3 oldPoint = startPoint;
	VECTOR3 oldPointBoxid = startPointBoxid;
	VECTOR3 nextPoint;
	VECTOR3 nextPointBoxid;

	#ifdef DEBUG_MODE
	fprintf( stderr, "Box Length: %g %g %g\n", this->boxLen[iRes][0], this->boxLen[iRes][1], this->boxLen[iRes][2] );
	fprintf( stderr, "dir: %g %g %g\n", dir[0], dir[1], dir[2] );
	fprintf( stderr, "tDel: %g %g %g\n", tDelta[0], tDelta[1], tDelta[2] );
	fprintf( stderr, "tMax: %g %g %g\n", tMax[0], tMax[1], tMax[2] );
	#endif

	while( !( oldPointBoxid == endPointBoxid ) )
	{
		if( tMax[0] < tMax[1] )
		{
			if( tMax[0] < tMax[2] )
			{
				nextPointBoxid.Set( oldPointBoxid[0] + step[0], oldPointBoxid[1], oldPointBoxid[2] );
				tMax.add( tDelta[0], 0, 0 );
			}
			else
			{
				nextPointBoxid.Set( oldPointBoxid[0], oldPointBoxid[1], oldPointBoxid[2] + step[2] );
				tMax.add( 0, 0, tDelta[2] );
			}
		}
		else
		{
			if( tMax[1] < tMax[2] )
			{
				nextPointBoxid.Set( oldPointBoxid[0], oldPointBoxid[1] + step[1], oldPointBoxid[2] );
				tMax.add( 0, tDelta[1], 0 );
			}
			else
			{
				nextPointBoxid.Set( oldPointBoxid[0], oldPointBoxid[1], oldPointBoxid[2] + step[2] );
				tMax.add( 0, 0, tDelta[2] );
			}
		}// end if-else

		#ifdef DEBUG_MODE
		fprintf( stderr, "Next point box id at current iteration: %g %g %g\n", nextPointBoxid[0], nextPointBoxid[1], nextPointBoxid[2] );
		#endif

		// Check if next point is within data grid
		/*
		if( nextPointBoxid[0] < 0 || nextPointBoxid[0] >= this->resolutions[iRes].x()
		 || nextPointBoxid[1] < 0 || nextPointBoxid[1] >= this->resolutions[iRes].y()
		 || nextPointBoxid[2] < 0 || nextPointBoxid[2] >= this->resolutions[iRes].z() )
		{
			//fprintf( stderr, "Start point: %g %g %g\n", startPoint[0], startPoint[1], startPoint[2] );
			//fprintf( stderr, "End point: %g %g %g\n", endPoint[0], endPoint[1], endPoint[2] );
			//fprintf( stderr, "Start point box id: %g %g %g\n", startPointBoxid[0], startPointBoxid[1], startPointBoxid[2] );
			//fprintf( stderr, "End point box id at current iteration: %g %g %g\n", endPointBoxid[0], endPointBoxid[1], endPointBoxid[2] );
			exit(0);
		}
		*/

		// Insert points to box count set, if unique
		if( countCheckArray[this->boxId3DTo1D( nextPointBoxid, iRes )] == false )
		{
			countCheckArray[this->boxId3DTo1D( nextPointBoxid, iRes )] = true;
			boxCount ++;
		}


		// Update information for next iteration
		//oldPoint = nextPoint;
		oldPointBoxid = nextPointBoxid;

	}// end while

	return boxCount ++;

}// end function

void
FDL_boxcountdim::transformGeometry( VECTOR3* trace, int traceLength )
{
	float xmin = 10000, ymin = 10000, zmin = 10000;
	float xmax = -10000, ymax = -10000, zmax = -10000;
	
	// Apply hotelling transform
	list<VECTOR3*>* lpv3Trace = new list<VECTOR3*>();
	list<VECTOR3*>* lpv3NewTrace = new list<VECTOR3*>();
	
	// Copy data
	for(int iT = 0; iT < traceLength; iT++)
	{
		lpv3Trace->push_back( new VECTOR3(trace[iT].x(), trace[iT].y(), trace[iT].z()) );
		lpv3NewTrace->push_back( new VECTOR3(0.0f, 0.0f, 0.0f) );
	}

	// Transform
	//htApply( *lpv3Trace, *lpv3NewTrace );
	//htApply<VECTOR3, 3>( *lpv3Trace, *lpv3NewTrace );

	// Copy transformed data back
	list<VECTOR3*>::iterator iter = lpv3NewTrace->begin();
	for(int iT=0; iT < traceLength; iT++, iter++ )
	{
		VECTOR3 *ptr = *iter;
		trace[iT].Set( ptr->x(), ptr->y(), ptr->z() );	
	}

	// Free resources
	for( list<VECTOR3*>::iterator iter1 = lpv3Trace->begin();
		 iter1 != lpv3Trace->end(); 
		 iter1++ ) 
		delete *iter1;
	for( list<VECTOR3*>::iterator iter2 = lpv3NewTrace->begin();
		 iter2 != lpv3NewTrace->end(); 
		 iter2++ ) 
		delete *iter2;
	delete lpv3Trace;
	delete lpv3NewTrace;

	// Find bounding box that can fit the streamline
	for( int iT=0; iT<traceLength; iT++ )
	{
		if( trace[iT].x() < xmin ) xmin = trace[iT].x();
		if( trace[iT].y() < ymin ) ymin = trace[iT].y();
		if( trace[iT].z() < zmin ) zmin = trace[iT].z();

		//if( trace[iT].x() > xmax ) xmax = trace[iT].x();
		//if( trace[iT].y() > ymax ) ymax = trace[iT].y();
		//if( trace[iT].z() > zmax ) zmax = trace[iT].z();
	}

	// Apply translation
	for( int iT=0; iT<traceLength; iT++ )
		trace[iT].Set( trace[iT].x()-xmin, trace[iT].y()-ymin, trace[iT].z()-zmin ); 

}// end function

void
FDL_boxcountdim::transformHotelling( VECTOR3* trace, VECTOR3* newTrace, int traceLength )
{
	assert( trace != NULL );
	assert( newTrace != NULL );

	// Apply hotelling transform
	list<VECTOR3*>* lpv3Trace = new list<VECTOR3*>();
	list<VECTOR3*>* lpv3NewTrace = new list<VECTOR3*>();

	// Copy data
	for(int iT = 0; iT < traceLength; iT++)
	{
		#if DEBUG_MODE
		fprintf( stderr, "%d: [%g %g %g]\n", iT, trace[iT].x(), trace[iT].y(), trace[iT].z() );
		#endif

		lpv3Trace->push_back( new VECTOR3(trace[iT].x(), trace[iT].y(), trace[iT].z()) );
		lpv3NewTrace->push_back( new VECTOR3(0.0f, 0.0f, 0.0f) );
	}

	// Transform
	htApply<VECTOR3, 3>( *lpv3Trace, *lpv3NewTrace );

	// Copy transformed data back
	list<VECTOR3*>::iterator iter = lpv3NewTrace->begin();
	for(int iT=0; iT < traceLength; iT++, iter++ )
	{
		VECTOR3 *ptr = *iter;
		newTrace[iT].Set( ptr->x(), ptr->y(), ptr->z() );
	}

	// Free resources
	for( list<VECTOR3*>::iterator iter1 = lpv3Trace->begin();
		 iter1 != lpv3Trace->end();
		 iter1++ )
		delete *iter1;
	for( list<VECTOR3*>::iterator iter2 = lpv3NewTrace->begin();
		 iter2 != lpv3NewTrace->end();
		 iter2++ )
		delete *iter2;
	delete lpv3Trace;
	delete lpv3NewTrace;

}

void
FDL_boxcountdim::transformHotelling2( VECTOR3* cHullTrace,
										   VECTOR3* trace,
										   VECTOR3* newTrace,
										   int hullLength,
										   int traceLength )
{
	assert( cHullTrace != NULL );
	assert( trace != NULL );
	assert( newTrace != NULL );

	// Apply hotelling transform
	list<VECTOR3*>* lpv3Trace = new list<VECTOR3*>();
	list<VECTOR3*>* lpv3HullTrace = new list<VECTOR3*>();
	list<VECTOR3*>* lpv3NewTrace = new list<VECTOR3*>();
	list<VECTOR3*>* lpv3NewHullTrace = new list<VECTOR3*>();

	// Copy data
	for(int iT = 0; iT < traceLength; iT++)
	{
		lpv3Trace->push_back( new VECTOR3(trace[iT].x(), trace[iT].y(), trace[iT].z()) );
		lpv3NewTrace->push_back( new VECTOR3(0.0f, 0.0f, 0.0f) );
	}
	// Copy hull
	for(int iT = 0; iT < hullLength; iT++)
	{
		lpv3HullTrace->push_back( new VECTOR3( cHullTrace[iT].x(), cHullTrace[iT].y(), cHullTrace[iT].z()) );
		lpv3NewHullTrace->push_back( new VECTOR3(0.0f, 0.0f, 0.0f) );
	}

	// Transform
	htApply2<VECTOR3, 3>( *lpv3HullTrace, *lpv3Trace, *lpv3NewHullTrace, *lpv3NewTrace );

	// Copy transformed data back
	list<VECTOR3*>::iterator iter = lpv3NewTrace->begin();
	for(int iT=0; iT < traceLength; iT++, iter++ )
	{
		VECTOR3 *ptr = *iter;
		newTrace[iT].Set( ptr->x(), ptr->y(), ptr->z() );
	}

	// Free resources
	for( list<VECTOR3*>::iterator iter1 = lpv3Trace->begin();
		 iter1 != lpv3Trace->end();
		 iter1++ )
		delete *iter1;
	for( list<VECTOR3*>::iterator iter2 = lpv3NewTrace->begin();
		 iter2 != lpv3NewTrace->end();
		 iter2++ )
		delete *iter2;
	for( list<VECTOR3*>::iterator iter3 = lpv3HullTrace->begin();
		 iter3 != lpv3HullTrace->end();
		 iter3++ )
		delete *iter3;
	for( list<VECTOR3*>::iterator iter4 = lpv3NewHullTrace->begin();
		 iter4 != lpv3NewHullTrace->end();
		 iter4++ )
		delete *iter4;

	delete lpv3Trace;
	delete lpv3NewTrace;
	delete lpv3HullTrace;
	delete lpv3NewHullTrace;
}

void
FDL_boxcountdim::transformPCA( VECTOR3* trace, VECTOR3* newTrace, int traceLength )
{
	float* data = new float[traceLength*3];
	for(int iT = 0; iT < traceLength; iT++ )
	{
		data[iT*3] = trace[iT].x();
		data[iT*3+1] = trace[iT].y();
		data[iT*3+2] = trace[iT].z();
	}

	cv::Mat streamline( traceLength, 3, CV_32F, data );

	// Perform a PCA:
	PCA pca( streamline, Mat(), CV_PCA_DATA_AS_ROW, 3 );

	//OutputArray transformedLine( traceLength, 3, CV_32F );

	Mat coeffs;
	Mat reconStreamline;
	//for( int iT = 0; iT < traceLength; iT++ )
	//{
		//Mat vec = streamline.row(iT), coeffs = transformedLine.row(iT);
		// compress the vector, the result will be stored
		// in the i-th row of the output matrix
		pca.project( streamline, coeffs);
		// and then reconstruct it
		//pca.backProject( coeffs, reconStreamline );
		// and measure the error
		//printf("%d. diff = %g\n", i, norm(vec, reconstructed, NORM_L2));

	for(int iT = 0; iT < traceLength; iT++ )
	{
		//newTrace[iT].Set( reconStreamline.at<float>(iT,0),
		//				  reconStreamline.at<float>(iT,1),
		//				  reconStreamline.at<float>(iT,2) );
		newTrace[iT].Set( coeffs.at<float>(iT,0),
							coeffs.at<float>(iT,1),
							coeffs.at<float>(iT,2) );

		#if DEBUG_MODE
		fprintf( stderr, "[%g %g %g] - [%g %g %g]\n",
				trace[iT].x(), trace[iT].y(), trace[iT].z(),
				newTrace[iT].x(), newTrace[iT].y(), newTrace[iT].z() );
		#endif
	}

	streamline.release();
	reconStreamline.release();
	delete [] data;
}

void
FDL_boxcountdim::transformToLocal_FirstPointBased( VECTOR3* trace, VECTOR3* newTrace, int traceLength )
{
	// Apply translation
	for( int iT=0; iT<traceLength; iT++ )
		newTrace[iT].Set( trace[iT].x()-trace[0].x(),
					   	  trace[iT].y()-trace[0].y(),
					   	  trace[iT].z()-trace[0].z() );

}// end function

void
FDL_boxcountdim::transformToLocal_CenterBased( VECTOR3* trace, VECTOR3* newTrace, int traceLength )
{
	float xmin = 0.0, ymin = 0.0, zmin = 0.0;
	float xmax = 0.0, ymax = 0.0, zmax = 0.0;
	float  center[3];

	// Find bounding box that can fit the streamline
	for( int iT=0; iT<traceLength; iT++ )
	{
		if( iT == 0 )
		{
			xmin = trace[iT].x();			xmax = trace[iT].x();
			ymin = trace[iT].y();			ymax = trace[iT].y();
			zmin = trace[iT].z();			zmax = trace[iT].z();
		}

		if( trace[iT].x() < xmin ) xmin = trace[iT].x();
		if( trace[iT].y() < ymin ) ymin = trace[iT].y();
		if( trace[iT].z() < zmin ) zmin = trace[iT].z();

		if( trace[iT].x() > xmax ) xmax = trace[iT].x();
		if( trace[iT].y() > ymax ) ymax = trace[iT].y();
		if( trace[iT].z() > zmax ) zmax = trace[iT].z();
	}

	center[0] = ( xmin + xmax ) / 2.0f;
	center[1] = ( ymin + ymax ) / 2.0f;
	center[2] = ( zmin + zmax ) / 2.0f;

	// Apply translation
	for( int iT=0; iT<traceLength; iT++ )
		newTrace[iT].Set( trace[iT].x()-center[0],
					   	  trace[iT].y()-center[1],
					   	  trace[iT].z()-center[2] );


}// end function

void
FDL_boxcountdim::transformTranslate( VECTOR3* trace, VECTOR3* newTrace,
										  int traceLength, VECTOR3* delta )
{
	#ifdef DEBUG_MODE
	fprintf( stderr, "delta: %g %g %g ...\n", delta->x(), delta->y(), delta->z() );
	#endif

	// Find bounding box that can fit the streamline
	//for( int iT=0; iT<traceLength; iT++ )
	//{
	//	if( trace[iT].x() < xmin ) xmin = trace[iT].x();
	//	if( trace[iT].y() < ymin ) ymin = trace[iT].y();
	//	if( trace[iT].z() < zmin ) zmin = trace[iT].z();

		//if( trace[iT].x() > xmax ) xmax = trace[iT].x();
		//if( trace[iT].y() > ymax ) ymax = trace[iT].y();
		//if( trace[iT].z() > zmax ) zmax = trace[iT].z();
	//}

	// Apply translation
	for( int iT=0; iT<traceLength; iT++ )
		newTrace[iT].Set( trace[iT].x() - delta->x(),
					   	  trace[iT].y() - delta->y(),
					   	  trace[iT].z() - delta->z() );

}// end function

VECTOR3
FDL_boxcountdim::transformMeanShift( int traceLength, VECTOR3* trace, VECTOR3* newTrace )
{
	VECTOR3 mu( 0.0f, 0.0f, 0.0f );
	fprintf( stderr, "len: %d\n", traceLength );

	// Compute mean
	for( int iT=0; iT<traceLength; iT++ )
	{
		//fprintf( stderr, "trace: %g %g %g ...\n", trace[iT].x(), trace[iT].y(), trace[iT].z() );
		mu[0] += trace[iT].x();
		mu[1] += trace[iT].y();
		mu[2] += trace[iT].z();
	}
	fprintf( stderr, "mu: %g %g %g ...\n", mu.x(), mu.y(), mu.z() );
	mu.Set( mu.x() / traceLength, mu.y() / traceLength, mu.z() / traceLength );

	//#ifdef DEBUG_MODE
	fprintf( stderr, "mu: %g %g %g ...\n", mu.x(), mu.y(), mu.z() );
	//#endif

	// Apply translation
	for( int iT=0; iT<traceLength; iT++ )
		newTrace[iT].Set( trace[iT].x() - mu[0],
					   	  trace[iT].y() - mu[1],
					   	  trace[iT].z() - mu[2] );

	return mu;

}// end function


void
FDL_boxcountdim::transformRotate( VECTOR3* trace, VECTOR3* newTrace,
									  int traceLength, MATRIX3* rot )
{
	// Find bounding box that can fit the streamline
	//for( int iT=0; iT<traceLength; iT++ )
	//{
	//	if( trace[iT].x() < xmin ) xmin = trace[iT].x();
	//	if( trace[iT].y() < ymin ) ymin = trace[iT].y();
	//	if( trace[iT].z() < zmin ) zmin = trace[iT].z();

		//if( trace[iT].x() > xmax ) xmax = trace[iT].x();
		//if( trace[iT].y() > ymax ) ymax = trace[iT].y();
		//if( trace[iT].z() > zmax ) zmax = trace[iT].z();
	//}

	// Apply translation
	for( int iT=0; iT<traceLength; iT++ )
		newTrace[iT] = trace[iT] * (*rot);

}// end function


