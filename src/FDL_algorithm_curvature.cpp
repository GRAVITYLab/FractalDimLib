/*
 * FDL_algorithm_curvature.cpp
 *
 *  Created on: Feb 2, 2011
 *      Author: abon
 */

#include "FDL_algorithm_curvature.h"

FDL_algorithm_curvature::FDL_algorithm_curvature()
{
	curvatureComputer = NULL;
	
	boxBoundary = NULL;
	scoreList = NULL;
	sortedFeatureList = NULL;
	
}// end constructor

void
FDL_algorithm_curvature::setCurvatureComputer( FDL_curvature *cvcomputer )
{
	curvatureComputer = cvcomputer;
}// end function


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//													Functions related to computation of feature map
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
FDL_algorithm_curvature::compute_FeatureMap_SingleStreamline( FDL_streamline *streamline, int minLen )
{
	assert( streamline != NULL );
	assert( curvatureComputer != NULL );
	
	// Get the trace and trace length of the streamline
	int totalTraceLength = streamline->getLength();
	VECTOR3* totalTrace = new VECTOR3[totalTraceLength];
	streamline->getTotalTrace( totalTrace );

	// Based on the streamline's length and min length, determine how many levels to go
	int maxLevel = (int)( floor( log( (double)totalTraceLength ) / log( 2.0f ) ) - floor( log( (double)minLen ) / log( 2.0f ) ) );
	
	// Allocate memory for feature map if needed
	FDL_featuremap featureMap;

	// Initialize the recursive call 	
	compute_FeatureMap_SingleStreamline_Recursive( &featureMap, totalTrace, 0, 0, totalTraceLength, maxLevel );
	
	// Handle special case where streamline's total length is smaller than minimum threshold
	// No segment has been created by the  resursive function
	if( featureMap.hierarchicalSegmentList.size() == 0 )
	{
		// Create a segment containing entire streamline
		SLSegment segment = { 0, 0, totalTraceLength,
								curvatureComputer->computeCurvature_SingleStreamline_SingleWindow( totalTrace, totalTraceLength )  };
		
		// Store information about this segment to the feature map
		featureMap.add_Segment( segment );
	}
	
	// Sort the segments in the feature map from top level to bottom
	featureMap.sort_Segments_Level_Start();

	// Set number of levels and segments
	featureMap.setNumLevels( featureMap.hierarchicalSegmentList.back().level+1 );
	featureMap.setNumSegments( featureMap.hierarchicalSegmentList.size() );
	
	// Assign this feature map to the streamline
	streamline->setFeatureMap( featureMap  );
	
	#ifdef DEBUG_MODE
	fprintf( stderr, "Feature map of the current streamline has %d levels and %d segments.\n", nLevel, nSegment );
	#endif

	delete [] totalTrace;

}// end function

void
FDL_algorithm_curvature::compute_FeatureMap_SingleStreamline_Recursive( FDL_featuremap *featureMap, VECTOR3 *trace, int level, int start,  int traceLength, int maxLevel )
{
	float curSegmentScore = 0;
	
	// Return if trace length smaller than threshold, otherwise proceed to further subdivision
	if( level > maxLevel )	return;
		
	// Compute Score for this segment
	curSegmentScore = curvatureComputer->computeCurvature_SingleStreamline_SingleWindow( trace, traceLength );
	
	// Store info about this segment
	SLSegment curSegment = { level, start, traceLength, curSegmentScore };
	featureMap->add_Segment( curSegment );
	
	// Divide the trace into two halves and compute FD of each half
	int firstHalfLength = traceLength / 2;
	int secondHalfLength = traceLength - firstHalfLength + 1;
	
	// Recursively call the same function on the two halves	
	compute_FeatureMap_SingleStreamline_Recursive( featureMap, trace, level+1, start, firstHalfLength, maxLevel );
	compute_FeatureMap_SingleStreamline_Recursive( featureMap, trace+firstHalfLength-1, level+1, start+firstHalfLength-1, secondHalfLength, maxLevel );
	
}// end function

void
FDL_algorithm_curvature::compute_FeatureMap_MultiStreamline( FDL_streamline *streamlineList, int nStreamline, int minLen )
{
	assert( streamlineList != NULL );
	
	for( int iS=0; iS<nStreamline; iS++ )
	{
		fprintf( stderr, "Computing feature map of streamline %d\n", iS );
		compute_FeatureMap_SingleStreamline( streamlineList+iS, minLen );
	}

}// end function

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//											Feature based segmentation of streamlines 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
FDL_algorithm_curvature::segment_FeatureMapBased_Recursive_SingleStreamline( FDL_streamline *streamline  )
{
	assert( streamline != NULL );
	
	// Get the trace and trace length of the streamline
	int totalTraceLength = streamline->getLength();
	VECTOR3* totalTrace = new VECTOR3[totalTraceLength];
	streamline->getTotalTrace( totalTrace );

	// Get pointer to the featuremap of the streamline
	FDL_featuremap featureMap = streamline->getFeatureMap();
	
	// This function assumes that the streamline already has a valid computed feature map
	//assert( featureMap != NULL );
	assert( featureMap.getNumSegments() > 0 );
	
	// Clear segment list
	streamline->featureBasedSegmentList.clear();

	// Sort all segments in feature map in increasing order of start point (and then decreasing order of score)
	list<SLSegment> tempList = featureMap.hierarchicalSegmentList;
	tempList.sort( Segment::compare_Segment_Start_Score2 ); 
	
	// Create a shorter list containing ONLY the highest scoring segment from each streamline
	list<SLSegment> highScoreSegmentsList;
	int laststart = -1;
	for( list<SLSegment>::const_iterator pcIter1 = tempList.begin(); pcIter1 != tempList.end(); pcIter1++ )
	{
		if( (*pcIter1).start != laststart ) 
		{
			highScoreSegmentsList.push_back( *pcIter1 );
			laststart = (*pcIter1).start;
		}
	}	
	
	#ifdef DEBUG_MODE
	fprintf( stderr, "printing highest scoring segment at each starting point of the feature map ...\n" );
	Segment::printSegmentList( highScoreSegmentsList );
	#endif	

	// Recursive call to include segments from start to end
	SLSegment matchFullStreamline = featureMap.filter_Segment_Start_Length( 0, totalTraceLength );
	getNextSegment( streamline, highScoreSegmentsList, matchFullStreamline, totalTraceLength );

	// Sort all segments of the output segment list based on start point 
	streamline->featureBasedSegmentList.sort( Segment::compare_Segment_Start_Score2 ); 

	#ifdef DEBUG_MODE
	fprintf( stderr, "printing results of segmentation based on feature map ...\n" );
	Segment::printSegmentList( streamline->featureBasedSegmentList );
	#endif

}// end function

void
FDL_algorithm_curvature::segment_FeatureMapBased_Iterative_SingleStreamline( FDL_streamline *streamline  )
{
	assert( streamline != NULL );
	
	// Get pointer to the featuremap of the streamline
	FDL_featuremap featureMap = streamline->getFeatureMap();
	
	// This function assumes that the streamline already has a valid computed feature map
	//assert( featureMap != NULL );
	
	// Clear segment list
	streamline->featureBasedSegmentList.clear();

	// Sort all segments of the chosen resolution pair in increasing order of start point (and then decreasing order of score)
	list<SLSegment> tempList = featureMap.hierarchicalSegmentList;
	tempList.sort( SLSegment::compare_Segment_Start_Score2 ); 
	
	// Create a shorter list containing ONLY the highest scoring segment from each streamline
	list<SLSegment> highScoreSegmentsList;
	int laststart = -1;
	for( list<SLSegment>::const_iterator pcIter1 = tempList.begin(); pcIter1 != tempList.end(); pcIter1++ )
	{
		if( (*pcIter1).start != laststart ) 
		{
			highScoreSegmentsList.push_back( *pcIter1 );
			laststart = (*pcIter1).start;
		}
	}	
	
	#ifdef DEBUG_MODE
	fprintf( stderr, "printing highest scoring segment at each starting point of the feature map ...\n" );
	Segment::printSegmentList( highScoreSegmentsList );
	#endif

	// Create a shorter list containing ONLY the segments at from each streamline (Use just one grid resolution)
	list<struct Segment> highestLevelSegmentsList;
	int maxLevel = featureMap.getNumLevels() - 1;
	for( list<SLSegment>::const_iterator pcIter1 = tempList.begin(); pcIter1 != tempList.end(); pcIter1++ )
	{
		if( (*pcIter1).level == maxLevel ) 
			highestLevelSegmentsList.push_back( *pcIter1 );
	}	
	
	// Iterate through the high scoring segments
	for( list<SLSegment>::const_iterator pcIter2 = highestLevelSegmentsList.begin(); pcIter2 != highestLevelSegmentsList.end(); pcIter2++ )
	{
		// Find the maximum value at that start point
		float max = (*pcIter2).score;
		for( list<SLSegment>::const_iterator pcIter1 = highScoreSegmentsList.begin(); pcIter1 != highScoreSegmentsList.end(); pcIter1++ )
		{
			if( SLSegment::hasOverlap( *pcIter1,  *pcIter2 ) &&  (*pcIter1).score > max )
				max = (*pcIter1).score;			
		}
		
		SLSegment newsegment = { maxLevel, (*pcIter2).start, (*pcIter2).length, max }; 
		streamline->featureBasedSegmentList.push_back( newsegment );
			
	}// end for
	
	// Sort all segments of the output segment list based on start point 
	streamline->featureBasedSegmentList.sort( SLSegment::compare_Segment_Start_Score2 ); 

	#ifdef DEBUG_MODE
	fprintf( stderr, "printing results of segmentation based on feature map ...\n" );
	Segment::printSegmentList( streamline->featureBasedSegmentList );
	#endif

}// end function

void
FDL_algorithm_curvature::segment_FeatureMapBased_MultiStreamline( FDL_streamline *streamlineList, int nStreamline, int minLen )
{
	assert( streamlineList != NULL );
	
	for( int iS=0; iS<nStreamline; iS++ )
	{
		// Clear any stored segmentation


		//fprintf( stderr, "Segmenting streamline %d ... \n", iS );
		segment_FeatureMapBased_Iterative_SingleStreamline( streamlineList + iS );
		joinSegments_SingleStreamline( streamlineList + iS );
	}

}// end function


void 
FDL_algorithm_curvature::joinSegments_SingleStreamline( FDL_streamline *streamline )
{
	// Get the trace and trace length of the streamline
	int totalTraceLength = streamline->getLength();
	VECTOR3* totalTrace = new VECTOR3[totalTraceLength];
	streamline->getTotalTrace( totalTrace );

	// Get segmentation list
	list<SLSegment> tempList = streamline->featureBasedSegmentList;
	
	// Clear compressed segment list of streamline
	streamline->joinedFeatureBasedSegmentList.clear();

	// Handle special case: only one segment in the streamline 
	if( tempList.size() == 1 )
	{
		streamline->joinedFeatureBasedSegmentList = streamline->featureBasedSegmentList;
		return;
	}

	// Iterate through the copied list of segments and join them
	list<SLSegment>::const_iterator pcIter = tempList.begin();
	SLSegment lastSegment = (*pcIter);
	SLSegment nextSegment;
	int stretchStart = 0;
	int stretchLen = 1;
	float stretchScore = 0;
	float tmpScore = 0;
	float joinThreshold = 0.005f;
	for(; pcIter != tempList.end(); pcIter++ )
	{
		// Compare with last segment
		// If the values are similar, continue to be in the strech
		// Otherwise break and create the next segment
		if( fabs( lastSegment.score - (*pcIter).score ) >= 0.1 * lastSegment.score  )
		{
			// Create a segment upto this point
			stretchScore = curvatureComputer->computeCurvature_SingleStreamline_SingleWindow( &totalTrace[stretchStart], stretchLen );

#if defined( _WIN32 ) || defined( _WIN64 )
			nextSegment.level = 0;
			nextSegment.start = stretchStart;
			nextSegment.length = stretchLen;
			nextSegment.score = stretchScore;
#endif
#if defined (_LINUX)
			nextSegment = { boxdimComputer->resType, stretchStart, stretchLen, stretchScore }; 
#endif

			// Add segment to the streamline
			streamline->joinedFeatureBasedSegmentList.push_back( nextSegment );
			
			// Start a new segment
			stretchStart = (*pcIter).start;
			stretchLen = (*pcIter).length;
		}
		else
		{
			// Continue this stretch
			stretchLen += (*pcIter).length-1;
		}
		
		// Remember last segment
		lastSegment = (*pcIter);
			
	}// end for
	
	// Add the last stretch to the streamline
	stretchScore = curvatureComputer->computeCurvature_SingleStreamline_SingleWindow( &totalTrace[stretchStart], stretchLen );

#if defined( _WIN32 ) || defined( _WIN64 )
	nextSegment.level = 0;
	nextSegment.start = stretchStart;
	nextSegment.length = stretchLen;
	nextSegment.score = stretchScore;
#endif
#if defined (_LINUX)
	nextSegment = { boxdimComputer->resType, stretchStart, stretchLen, stretchScore }; 
#endif

	streamline->joinedFeatureBasedSegmentList.push_back( nextSegment );

	#ifdef DEBUG_MODE
	fprintf( stderr, "printing results of joining segments ...\n" );
	Segment::printSegmentList( streamline->compressedFeatureBasedSegmentList );
	#endif	

	// Clear
	delete [] totalTrace;

}// end function

int
FDL_algorithm_curvature::getNextSegment( FDL_streamline *streamline,
											   list<SLSegment> tempList, SLSegment segment,
											   int totalLength )
{
	// Get the trace and trace length of the streamline
	int totalTraceLength = streamline->getLength();
	VECTOR3* totalTrace = new VECTOR3[totalTraceLength];
	streamline->getTotalTrace( totalTrace );

	// Pointers to feature map and destination list
	FDL_featuremap featureMap = streamline->getFeatureMap();
	list<SLSegment> *destList = &streamline->featureBasedSegmentList;
	
	//fprintf( stderr, "Recursive call with <start: %d, end: %d, score:%f>\n", segment.start, segment.start+segment.length-1, segment.bcd );
	if( segment.start >= totalLength-1 )
		return ( segment.start + segment.length - 1 );
	
	list<SLSegment>::const_iterator pcIter = tempList.begin();
	int nextstart = 0;
	for(; pcIter != tempList.end(); pcIter++ )
	{
		//fprintf( stderr, "In for loop with <start: %d, end: %d, score:%f>\n", (*pcIter).start, (*pcIter).start+(*pcIter).length-1, (*pcIter).bcd );
		//std::cin.get();
		
		// Skip if the segment in the list has its start point before the segment being examined 
		if( (*pcIter).start < segment.start )
			continue;
		
		// Skip if the segment in the list is same as the segment being examined
		if( (*pcIter).start == segment.start || (*pcIter).length == segment.length  )
			continue;
		
		// If the segment in the list has its start point after the end point of the segment being examined,
		// No more segment in the list needs to be compared with. Return.
		if( (*pcIter).start >= segment.start+segment.length-1 )
		{
			// No sub-segment with higher score was found
			// Add this full segment to the final list
			destList->push_back( segment );
	
			// Return the last point of this segment
			return ( segment.start + segment.length - 1 );
		}
				
		// If a sub-segment has a higher score than the current one, exit the loop and make a recursive call		
		if( (*pcIter).score > segment.score ) 
			break;
	}
		
	if( pcIter == tempList.end() )
	{
		destList->push_back( segment );

		// Return the last point of this segment
		return ( segment.start + segment.length - 1 );		
	}
		
	nextstart = (*pcIter).start;
	
	// Subdivide the segment under study into two 
	// Add the first part which has already been traversed to the final list
	if( nextstart != segment.start )
	{
		SLSegment firstPart = { 0,  segment.start, (*pcIter).start - segment.start + 1,
							curvatureComputer->computeCurvature_SingleStreamline_SingleWindow( &totalTrace[segment.start], (*pcIter).start - segment.start + 1 )  							};
	
		destList->push_back( firstPart );
	}
	
	while( nextstart < segment.start + segment.length - 1 )
	{		
		//fprintf( stderr, "Looking for the highest score segment with <start: %d>\n", nextstart );
		list<SLSegment>::const_iterator pcIter2;
		struct  Segment match;
		for( pcIter2 = tempList.begin(); pcIter2 != tempList.end(); pcIter2++ )
		{
			if( (*pcIter2).start == nextstart )
			{
				match = (*pcIter2);
				break;
			}
			
		}
							
		// Recursive call with remaining portion of the segment
		// Returns the last point upto which the recursion covers
		nextstart = getNextSegment( streamline, tempList, match, totalLength );
		//fprintf( stderr, "returned new start point: %d\n", nextstart );
	}

	// Clear
	delete [] totalTrace;

	// Return the last point of this segment
	return nextstart;
	
}// end function

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//															File I/O functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FDL_algorithm_curvature::save_SegmentationResult( FDL_streamline *streamlineList,
														 int nstreamline, const char *featureFileName  )
{
	// Open file for reading
	FILE* featureFile = fopen( featureFileName, "wb" );
	assert( featureFile != NULL );
	
	int id = 0;
	int nSegment = 0;
	list<SLSegment>::const_iterator pIter;
	for( int iS = 0; iS<nstreamline; iS++ )
	{
		nSegment = streamlineList[iS].joinedFeatureBasedSegmentList.size();
		fwrite( &nSegment, 1, sizeof(int), featureFile );
		for( pIter = streamlineList[iS].joinedFeatureBasedSegmentList.begin();
		 	 pIter != streamlineList[iS].joinedFeatureBasedSegmentList.end();
		 	 pIter++ )
		{
			int start = (*pIter).start;
			int length = (*pIter).length;
			float score = (*pIter).score;
			
			// Write start point
			fwrite( &start, 1, sizeof(int), featureFile );
			fwrite( &length, 1, sizeof(int), featureFile );
			fwrite( &score, 1, sizeof(float), featureFile );
		}
	}
	
	// Close file
	fclose( featureFile );
	
}// end function

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//						Functions related to flow field subdivision and computation of feature vectors from there (NOT needed now)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FDL_algorithm_curvature::subdivide_FlowField( VECTOR3 *min, VECTOR3 *max, int npartx, int nparty, int npartz  )
{
	// Compute the total number of boxes to be created
	nPartition[0] = npartx;
	nPartition[1] = nparty;
	nPartition[2] = npartz;
	nTotalBox = nPartition[0] * nPartition[1] * nPartition[2];
	cout << nTotalBox << endl;
	
	// Compute the dimension of  a box 
	float boxDim[3];
	boxDim[0] = ( max->x() - min->x() ) / (float)nPartition[0];
	boxDim[1] = ( max->y() - min->y() ) / (float)nPartition[1];
	boxDim[2] = ( max->z() - min->z() ) / (float)nPartition[2];
	cout << boxDim[0] << " " << boxDim[1] << " " << boxDim[2] << endl;	

	// Allocate memory for storing box boundaries
	if( boxBoundary != NULL )	delete [] boxBoundary;
	boxBoundary = new double[nTotalBox*6];
	
	// Subdivide
	int index = 0;
	for( int k=0; k<nPartition[2]; k++ )
	{
		for( int j=0; j<nPartition[1]; j++ )
		{
			for( int i=0; i<nPartition[0]; i++ )
			{
				// Compute boundary of the current box
				boxBoundary[index] = i*boxDim[0];
				boxBoundary[index+1] = j*boxDim[1];
				boxBoundary[index+2] = k*boxDim[2];
				boxBoundary[index+3] = (i+1)*boxDim[0];
				boxBoundary[index+4] = (j+1)*boxDim[1];
				boxBoundary[index+5] = (k+1)*boxDim[2];
				
				#ifdef DEBUG_MODE
				fprintf( stderr, "[%d %d %d]-th box: (%f %f %f), (%f %f %f)\n", i, j, k, boxBoundary[index], boxBoundary[index+1], boxBoundary[index+2], boxBoundary[index+3], boxBoundary[index+4], boxBoundary[index+5] );
				#endif

				// Increment index
				index += 6;
			}
		}
	}
	
}// end function

void
FDL_algorithm_curvature::computeFeatureList_AllBox( FDL_streamline *streamlinelist, int nstreamline )
{
	for( int iBox=0; iBox<nTotalBox; iBox++ )
		computeFeatureList_SingleBox( streamlinelist, nstreamline, iBox );
	fprintf( stderr, "List contains %d features.\n", partitionwiseFeatureList.size() );
	
}// end function

void
FDL_algorithm_curvature::computeFeatureList_SingleBox( FDL_streamline *streamlineList,
															  int nstreamline, int boxid )
{
	int totalTraceLength = 0;
	float score = 0;
	VECTOR3* totalTrace = NULL;
	//fprintf( stderr, "Box id: %d\n", boxid );
	
	// Iterate through all the streamlines
	for( int iS = 0; iS<nstreamline; iS++ )
	{
		// Get the trace and trace length of the streamline
		totalTraceLength = streamlineList[iS].getLength();
		totalTrace = new VECTOR3[totalTraceLength];
		streamlineList[iS].getTotalTrace( totalTrace );

		// Compute the fractal dimensnion of this streamline with respect to this axis aligned box
		score = curvatureComputer->computeCurvature_SingleStreamline_InAxisAlignedBox( totalTrace, totalTraceLength,
																					   boxBoundary + boxid*6, boxBoundary + boxid*6+3 );
		//fprintf( stderr, "FD of Streamline %d in %d-th box: %f\n", iS, boxid, score );
		
		// Store FD as a feature if this streamline is present 
		if( score > 1.0 )
		{
			struct feature_in_box nextFeature = { boxid, iS, score };
			partitionwiseFeatureList.push_back( nextFeature ); 
		}

		// Clear
		delete [] totalTrace;

	}// end fir

}// end function

void
FDL_algorithm_curvature::computeFeatureList_AllBox2( FDL_streamline *streamlineList,
															int nstreamline,
															double* scoreList )
{
	if( scoreList == NULL )	scoreList = new double[nTotalBox];

	for( int iBox=0; iBox<nTotalBox; iBox++ )
	{
		scoreList[iBox] = computeFeatureList_SingleBox2( streamlineList, nstreamline, iBox );
		//cout << iBox << " " << fdList[iBox] << endl;
	}	

}// end function

double
FDL_algorithm_curvature::computeFeatureList_SingleBox2( FDL_streamline *streamlineList,
															   int nstreamline, int boxid )
{
	int totalTraceLength = 0;
	double score = 0;
	VECTOR3* totalTrace = NULL;
	//fprintf( stderr, "Box id: %d\n", boxid );
	//fprintf( stderr, "Number of streamlines: %d\n", nstreamline );
	
	// Iterate through all the streamlines
	float maxScore = -1;	
	for( int iS = 0; iS<nstreamline; iS++ )
	{
		//fprintf( stderr, "Streamline ID: %d\n", iS );
		// Get the trace and trace length of the streamline
		totalTraceLength = streamlineList[iS].getLength();
		totalTrace = new VECTOR3[totalTraceLength];
		streamlineList[iS].getTotalTrace( totalTrace );
		
		// Compute the fractal dimensnion of this streamline with respect to this axis aligned box
		score = curvatureComputer->computeCurvature_SingleStreamline_InAxisAlignedBox( totalTrace, totalTraceLength,
																					   boxBoundary + boxid*6, boxBoundary + boxid*6+3 );
		//fprintf( stderr, "FD of Streamline %d in %d-th box: %f\n", iS, boxid, score );
		
		// Store FD as a feature if this streamline is present 
		if( score > maxScore )
		{
			maxScore = score;
		}

		// Clear
		delete [] totalTrace;

	}// end for

	return maxScore;

}// end function

/*
void
FDL_algorithm_curvature::computeFeatureVectors_AllBox_AllStreamline( FDL_streamline *streamlinelist, int nstreamline, float *featureVector )
{
	//fprintf( stderr, "%d\n", nTotalBox );
	for( int iS=0; iS<nstreamline; iS++ )
	{
		computeFeatureVectors_AllBox_SingleStreamline( streamlinelist+iS,  featureVector + iS*nTotalBox );
		fprintf( stderr, "Computed feature vector for %d-th streamline\n", iS );
	}
	
	for( int iS=0; iS<nstreamline; iS++ )
	{
		list<float> featureList( featureVector + iS*nTotalBox, featureVector + iS*nTotalBox + nTotalBox );
		featureList.sort( FDL_algorithm_curvature::comparator_feature2 );	

		// Copy back
		for( int iB=0; iB<nTotalBox; iB++ )
		{
			*(featureVector + iS*nTotalBox + iB) = featureList.front();
			featureList.pop_front();
		}
	}
	
	
}// end function

void
FDL_algorithm_curvature::computeFeatureVectors_AllBox_SingleStreamline( FDL_streamline *streamline, float *featureVector  )
{
	float score = 0;
	for( int iBox=0; iBox<nTotalBox; iBox++ )
	{
		score = boxdimComputer->computeBCD_SingleStreamline_InAxisAlignedBox( streamline->totalTrace, streamline->totalTraceLength,
																			       boxBoundary + iBox*6, boxBoundary + iBox*6+3 );
		//fprintf( stderr, "FD of Streamline in %d-th box: %f\n", iBox, score );
		
		// Store FD as a feature if this streamline is present 
		if( score > 1.0 ) *(featureVector+iBox) = score;
		else			 *(featureVector+iBox) = 0;
	}
	
}// end function

void
FDL_algorithm_curvature::saveFeatureVectors( float *featureVector, int nstreamline, const char *featureFileName  )
{
	// Open file for reading
	FILE* featureFile = fopen( featureFileName, "w" );
	assert( featureFile != NULL );
	
	int id = 0;
	for( int iS = 0; iS<nstreamline; iS++ )
	{
		for( int iB=0; iB<nTotalBox; iB++ )
		{
			if( iB == nTotalBox-1 )	ffprintf( stderr, featureFile, "%f\n", featureVector[id]  );
			else					ffprintf( stderr, featureFile, "%f, ", featureVector[id]  );
			id ++;
		}
	}
	
	// Close file
	fclose( featureFile );
	
}// end function

void
FDL_algorithm_curvature::projectFeatures_AllStreamlines( int nstreamline, int K, float *featureSpan, float *featureScore )
{
	
	list<struct feature_in_box> featureList_curStreamline;
	list<struct feature_in_box> featureList_topK_currentStreamline;
	
	for( int iS=0; iS<nstreamline; iS++ )
	{
		// Identify features for the next streamline
		featureList_curStreamline = listFeatures_SingleStreamline( iS );
		
		// Select top-K features for the next streamline
		featureList_topK_currentStreamline =  listTopKFeatures_SingleStreamline( featureList_curStreamline, K );
		
		// Project Top-k features into a span and a score
		featureScore[iS] = compute_TopKFeaturesProjection_SingleStreamline(  featureList_topK_currentStreamline, K );
		featureSpan[iS] = compute_FeatureSpan_SingleStreamline(  featureList_topK_currentStreamline, K );

		// Clear
		featureList_curStreamline.clear();
		featureList_topK_currentStreamline.clear();
		
		
		fprintf( stderr, "Feature selection and projection for %d-th streamline completed\n", iS );
 	}
	
}// end function


list<struct feature_in_box>
FDL_algorithm_curvature::listFeatures_SingleStreamline( int streamlineID )
{
	// Iterator to all features list
	list<struct feature_in_box>::const_iterator pIter = partitionwiseFeatureList.begin();
	
	// Create a new list for features from this streamline
	list<struct feature_in_box> featureList_currentStreamline;
	
	// Iterate over list of all features to select the ones with given streamline ID
	while( pIter != partitionwiseFeatureList.end() )
	{
		if( (*pIter).streamline_id == streamlineID ) 	
			featureList_currentStreamline.push_back( (*pIter) );
		
		// Increment to the next feature
		pIter++;
		
	}// end while
	
	// Return the list of features
	return featureList_currentStreamline;
	
}// end function

bool
FDL_algorithm_curvature::comparator_feature( struct feature_in_box first, struct feature_in_box second )
{	
	if( first.feature_score >= second.feature_score )	return true;
	return false;
	
}// end function

bool
FDL_algorithm_curvature::comparator_feature2( float first, float second )
{	
	if( first >= second )	return true;
	return false;
	
}// end function


list<struct feature_in_box> 
FDL_algorithm_curvature::listTopKFeatures_SingleStreamline( list<struct feature_in_box> featureList_currentStreamline, int K )
{
	// Get total number of features 
	int nFeature = featureList_currentStreamline.size();
	
	// Sort features in decending order
	featureList_currentStreamline.sort( FDL_algorithm_curvature::comparator_feature );
	
	// Allocate memory for top-K features and corresponding
	list<struct feature_in_box> featureList_topK_currentStreamline;
	
	// Extract first K from the list of all features
	// If the feature list contains less than K elements, put 0 for the rest
	struct feature_in_box dummy = { 0, 0, 0 };
	for( int i=0; i<K; i++ )
	{
		if( featureList_currentStreamline.empty() )
			featureList_topK_currentStreamline.push_back( dummy );	
		else
		{
			featureList_topK_currentStreamline.push_back( featureList_currentStreamline.front() );
			featureList_currentStreamline.pop_front();
		}		
	}
	
	#ifdef DEBUG_MODE
	fprintf( stderr, "Printing top-%d features: ", K );
	list<struct feature_in_box>::const_iterator pIter;
	for( pIter = featureList_topK_currentStreamline.begin(); pIter != featureList_topK_currentStreamline.end(); pIter++ ) 
		fprintf( stderr, "%f ", (*pIter).feature_score  );
	fprintf( stderr, "\n" );
	#endif
	
	return featureList_topK_currentStreamline;
	
}// end function

float
FDL_algorithm_curvature::compute_TopKFeaturesProjection_SingleStreamline(  list<struct feature_in_box>  featureList_topK_currentStreamline, int K )
{
	float score = 0;
	list<struct feature_in_box>::const_iterator pIter;
	int i=0;
	for( pIter = featureList_topK_currentStreamline.begin(); pIter != featureList_topK_currentStreamline.end(); pIter++ ) 
	{
		//fprintf( stderr, "%f\t", (*pIter).feature_score );
		score += (*pIter).feature_score * (1 - (float)i / (float)K );
		i++;
	}
	//fprintf( stderr, "\n" );
	
	return score;	
}// end function

float
FDL_algorithm_curvature::compute_FeatureSpan_SingleStreamline(  list<struct feature_in_box>  featureList_topK_currentStreamline, int K )
{
	float minBoxId[3];
	float maxBoxId[3];
	
	list<struct feature_in_box>::const_iterator pIter;
	int i=0;
	for( pIter = featureList_topK_currentStreamline.begin(); pIter != featureList_topK_currentStreamline.end(); pIter++ ) 
	{
		if( i>0 && (*pIter).feature_score == 0 )
			continue;
		
		// Convert 1D boxID to 3D
		int zid = (int)floor( (*pIter).box_id / (float)( nPartition[0] * nPartition[1] ) );
		int xy = (*pIter).box_id % ( nPartition[0] * nPartition[1] );
		int yid = (int)floor( xy / (float)nPartition[1] );
		int xid = xy % nPartition[1];
		//fprintf( stderr, "%d %d %d %d\n", (*pIter).box_id, xid, yid, zid );
		
		// Compute min max box ID
		if( i == 0 )
		{
			minBoxId[0] = maxBoxId[0] = xid;
			minBoxId[1] = maxBoxId[1] = yid;
			minBoxId[2] = maxBoxId[2] = zid;
		}
		else
		{
			if( xid < minBoxId[0] )	minBoxId[0] = xid;  if( xid > maxBoxId[0] )	maxBoxId[0] = xid;  
			if( yid < minBoxId[1] )	minBoxId[1] = yid;  if( yid > maxBoxId[1] )	maxBoxId[1] = yid;
			if( zid < minBoxId[2] )	minBoxId[2] = xid;  if( zid > maxBoxId[2] )	maxBoxId[2] = zid;  
		}	
		
		i++;
	}
	return abs( maxBoxId[0] - minBoxId[0] ) + abs( maxBoxId[1] - minBoxId[1] ) + abs( maxBoxId[2] - minBoxId[2] );	

}// end function

*/

