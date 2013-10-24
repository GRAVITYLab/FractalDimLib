#include "FDL_curvature.h"

FDL_curvature::FDL_curvature()
{
}// end constructor

FDL_curvature::~FDL_curvature()
{
}// end destructor


double
FDL_curvature::computeCurvature_SingleStreamline_Point( VECTOR3* trace )
{
	VECTOR3 tempV;
	float tempF;

	VECTOR3 prevT = ( trace[1] - trace[0] );
	prevT.Normalize();

	VECTOR3 curT = ( trace[2] - trace[1] );
	curT.Normalize();
	
	tempV = ( curT - prevT );
	
	return tempV.GetMag();
}// end function

double
FDL_curvature::computeCurvature_SingleStreamline_SingleWindow( VECTOR3* trace, int traceLength )
{
	// Handle the degenerate case first
	if( traceLength < 3 )
		return 0.0f;

	float meanCurvature = 0;
	for( int i=0; i<traceLength-2; i++ )
		meanCurvature += computeCurvature_SingleStreamline_Point( trace+i ); 	
	
	return meanCurvature / (float)(traceLength-2);
}// end function

double
FDL_curvature::computeCurvature_MultiStreamline_SingleWindow( FDL_streamline* streamlineList,
																	  int nStreamline,
																	  int *windowIdList, int* windowLengthList )
{
	return 0.0f;
}// end function

double
FDL_curvature::computeCurvature_SingleStreamline_Full( FDL_streamline* streamline,
															double *localCurvatureArray,
															double *minScore, double *maxScore )
{
	// Get the trace and trace length of the streamline
	int totalTraceLength = streamline->getLength();
	VECTOR3* totalTrace = new VECTOR3[totalTraceLength];
	streamline->getTotalTrace( totalTrace );

	// Handle the degenerate case first
	if( totalTraceLength < 3 )
		return 0.0f;

	float meanCurvature = 0, nextCurvature = 0;
	for( int i=0; i<totalTraceLength-2; i++ )
	{
		nextCurvature = computeCurvature_SingleStreamline_Point( &totalTrace[i] );
		meanCurvature += nextCurvature;
		
		// Optional steps
		if( localCurvatureArray != NULL ) localCurvatureArray[i] = nextCurvature;
		if( minScore != NULL && nextCurvature < *minScore )	*minScore = nextCurvature;
		if( maxScore != NULL && nextCurvature > *maxScore )	*maxScore = nextCurvature;

	}// end for

	delete [] totalTrace;

	return meanCurvature / (totalTraceLength-2);

}// end function

double
FDL_curvature::computeCurvature_SingleStreamline_InAxisAlignedBox( VECTOR3* trace,
																			int traceLength,
																			double* boxMin, double* boxMax )
{
	return 0.0f;
}// end function

double*
FDL_curvature::computeCurvature_SingleStreamline_windowed( FDL_streamline* streamline, int windowLength )
{
	return NULL;
}// end function

void
FDL_curvature::computeCurvature_SingleStreamline_DisjointSegments( FDL_streamline* streamline,
																			int segmentLength, int nSegment,
																			double* curvatureArray )
{
	int startIndex = 0;
	VECTOR3* localTrace = NULL;
	
	// Get the trace and trace length of the streamline
	int totalTraceLength = streamline->getLength();
	VECTOR3* totalTrace = new VECTOR3[totalTraceLength];
	streamline->getTotalTrace( totalTrace );

	#ifdef DEBUG_MODE	
	printf( "Streamline Length: %d\n", streamline->totalTraceLength );
	printf( "Number of segments: %d\n", nSegment );
	#endif
	
	// Clear up the fixed-length segment list of the current streamline
	streamline->fixedLengthSegmentList.clear();
	
	// Slide window until end of streamline reached
	for( int i=0; i<nSegment; i++ )
	{
		// Determine the actual window length
		if( startIndex + segmentLength > totalTraceLength )
			segmentLength = totalTraceLength - startIndex;

		#ifdef DEBUG_MODE
		printf( "Moving to next window of length %d from trace point %d\n", segmentLength, startIndex );
		#endif

		// Compute and store box count dimension for this segment
		curvatureArray[i] = computeCurvature_SingleStreamline_SingleWindow( &totalTrace[startIndex], segmentLength );

		// Update the fixed-length segment list of this streamline		
		struct Segment newsegment = { 0, startIndex, segmentLength, curvatureArray[i] }; 
		streamline->fixedLengthSegmentList.push_back( newsegment );

		// Move to the start point of next segment
		startIndex += segmentLength;
		
	}// end while

	delete [] totalTrace;

}// end function

void
FDL_curvature::computeCurvature_MultiStreamline_Full( FDL_streamline* streamlineList,
															 int nStreamline, double* curvatureArray,
															 double* globalMinScore, double* globalMaxScore )
{
	assert( curvatureArray != NULL );
	for( int iS=0; iS<nStreamline; iS++ )
	{
		curvatureArray[iS] = computeCurvature_SingleStreamline_Full( streamlineList + iS );
		if( globalMinScore != NULL && curvatureArray[iS] < *globalMinScore ) *globalMinScore = curvatureArray[iS];
		if( globalMaxScore != NULL && curvatureArray[iS] > *globalMaxScore ) *globalMaxScore = curvatureArray[iS];		
	}

}// end function

void
FDL_curvature::computeCurvature_MultiStreamline_DisjointSegments( FDL_streamline* streamlineList,
																		   int nStreamline, int segmentLength,
																		   int *nSegmentList, double** segmentedScoreList )
{
	for( int iS=0; iS<nStreamline; iS++ )
	{
		// Compute segment wise mean curvatures of the next streamline
		computeCurvature_SingleStreamline_DisjointSegments( streamlineList + iS, segmentLength, nSegmentList[iS], segmentedScoreList[iS] );
		//printf( "Segmentwise FD of %d-th streamline computed.\n", iS );
		
	}// end for

}// end function
