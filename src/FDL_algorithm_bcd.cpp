/*
 * FDL_algorithm_bcd.cpp
 *
 *  Created on: Feb 2, 2011
 *      Author: abon
 */

#include "FDL_algorithm_bcd.h"

FDL_algorithm_bcd::FDL_algorithm_bcd()
{
	boxdimComputer = NULL;
	boxBoundary = NULL;
	scoreList = NULL;
	sortedFeatureList = NULL;

}// end constructor

FDL_algorithm_bcd::FDL_algorithm_bcd( FDL_boxcountdim *boxcounter )
{
	boxdimComputer = NULL;
	boxBoundary = NULL;
	scoreList = NULL;
	sortedFeatureList = NULL;
	
	boxdimComputer = boxcounter;

}// end constructor

FDL_algorithm_bcd::FDL_algorithm_bcd( const FDL_algorithm_bcd& that )
{
	this->boxdimComputer = NULL;
	this->boxBoundary = NULL;
	this->scoreList = NULL;
	this->sortedFeatureList = NULL;

	if( that.boxdimComputer != NULL )
		this->boxdimComputer = that.boxdimComputer;

}

FDL_algorithm_bcd&
FDL_algorithm_bcd::operator= ( const FDL_algorithm_bcd& that )
{
	if (this != &that ) // protect against invalid self-assignment
	{
		this->boxdimComputer = NULL;
		this->boxBoundary = NULL;
		this->scoreList = NULL;
		this->sortedFeatureList = NULL;

		if( that.boxdimComputer != NULL )
			this->boxdimComputer = that.boxdimComputer;
	}
	// by convention, always return *this
	return *this;

}// end constructor

FDL_algorithm_bcd::~FDL_algorithm_bcd()
{
	if( boxBoundary != NULL )	delete [] boxBoundary;
	if( scoreList != NULL )	delete [] scoreList;

}// end destructor


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//													Functions related to computation of feature map
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FDL_algorithm_bcd::compute_FeatureMap_SingleStreamline( FDL_streamline *streamline, int minLen, int slId )
{
	int maxLevel = 0, traceLength = 0;
	VECTOR3* totalTrace = NULL;
	FDL_featuremap featureMap;

	assert( streamline != NULL );
	assert( boxdimComputer != NULL );
	
	// Get the trace and trace length of the streamline
	traceLength = streamline->getLength();
	totalTrace = new VECTOR3[traceLength];
	streamline->getTotalTrace( totalTrace );

	// Based on the streamline's length and min length, determine how many levels to go
	if( traceLength <= minLen )
		maxLevel = 0;
	else
		maxLevel = (int)( floor( log( (float)traceLength ) / log( 2.0f ) ) - floor( log( (float)minLen ) / log( 2.0f ) ) );

	#ifdef DEBUG_MODE
	fprintf( stderr, "Streamline Length: %d \n", streamline->getLength() );
	fprintf( stderr, "Max level id of feature map: %d\n", maxLevel );
	#endif

	// Initialize the recursive call
	// for subdivision and segment creation
	compute_FeatureMap_SingleStreamline_Recursive( slId,
												   &featureMap,
												   totalTrace,
												   0, 0,
												   traceLength, maxLevel );
	
	// Handle special case where streamline's total length is smaller than minimum threshold
	// No segment has been created by the  recursive function
	if( featureMap.hierarchicalSegmentList.size() == 0 )
	{
		// Create a segment containing entire streamline
		SLSegment segment = { slId,
							  0, 0,
							  traceLength,
							  boxdimComputer->computeBCD_SingleStreamline_SingleWindow( totalTrace,
									  	  	  	  	  	  	  	  	  	  	  	  	    traceLength )
							};
		
		// Store information about this segment to the feature map
		featureMap.add_Segment( segment );
	}
	
	// Sort the segments in the feature map from top level to bottom
	featureMap.sort_Segments_Level_Start();


	// Set number of levels and segments
	//int hIndex = FDL_featuremap::resolutionToArrayIndex( boxdimComputer->resType );

	//featureMap.setNumLevels( featureMap.hierarchicalSegmentList[hIndex].back().level+1 );
	featureMap.setNumLevels( maxLevel+1 );
	featureMap.setNumSegments( featureMap.hierarchicalSegmentList.size() );
	
	// Clear any feature map the streamline may already have
	streamline->clearFeatureMap();

	// Assign this feature map to the streamline
	streamline->setFeatureMap( featureMap  );
	
	#ifdef DEBUG_MODE
	fprintf( stderr, "Feature map of the current streamline has %d levels and %d segments.\n",
			featureMap.hierarchicalSegmentList.back().level+1,
			featureMap.hierarchicalSegmentList.size() );
	#endif
	
	// Clear
	delete [] totalTrace;

}// end function

void
FDL_algorithm_bcd::compute_FeatureMap_SingleStreamline_Recursive( int slId,
																  FDL_featuremap *featureMap,
																  VECTOR3 *trace,
																  int level, int start,
																  int traceLength, int maxLevel )
{
	double curSegmentScore = 0;
	//fprintf( stderr, "%d: %d\n", level, maxLevel );

	// Return if trace length smaller than threshold, otherwise proceed to further subdivision
	if( level > maxLevel )	return;

	#ifdef DEBUG_MODE
	for( int i=0; i<traceLength; i++ )
		fprintf( stderr, "%d: %g %g %g\n", i, trace[i][0], trace[i][1], trace[i][2] );
	#endif

	// Compute Score for this segment
	curSegmentScore = boxdimComputer->computeBCD_SingleStreamline_SingleWindow( trace, traceLength );
	//fprintf( stderr, "%g\n", curSegmentScore );
	
	// Store info about this segment
	SLSegment curSegment = { slId, level, start, traceLength, curSegmentScore };
	featureMap->add_Segment( curSegment );
	//fprintf( stderr, "%d, %d, %d\n", level, start, traceLength );
	
	// Divide the trace into two halves and compute FD of each half
	int firstHalfLength = traceLength / 2;
	int secondHalfLength = traceLength - firstHalfLength + 1;
	//fprintf( stderr, "%d, %d\n", firstHalfLength, secondHalfLength );
	
	// Recursively call the same function on the two halves	
	compute_FeatureMap_SingleStreamline_Recursive( slId, featureMap, trace, level+1, start, firstHalfLength, maxLevel );
	compute_FeatureMap_SingleStreamline_Recursive( slId, featureMap, trace+firstHalfLength-1, level+1, start+firstHalfLength-1, secondHalfLength, maxLevel );
	
}// end function

void
FDL_algorithm_bcd::compute_FeatureMap_MultiStreamline( FDL_streamline *streamlineList,
															  int nStreamline, int minLen )
{
	assert( streamlineList != NULL );
	
	for( int iS=0; iS<nStreamline; iS++ )
	{
		compute_FeatureMap_SingleStreamline( streamlineList+iS, minLen, iS );

		#ifdef DEBUG_MODE
		fprintf( stderr, "Computed feature map of %d-th streamline\n", iS );
		//fprintf( stderr, "Printing feature map of %d-th streamline\n", iS );
		//streamlineList[iS].getFeatureMap().printFeatureMap();
		#endif

	}

}// end function

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//											Feature based segmentation of streamlines 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
FDL_algorithm_bcd::segment_FeatureMapBased_Recursive_SingleStreamline( FDL_streamline *streamline, int slId )
{
	assert( streamline != NULL );
	
	// Get the trace and trace length of the streamline
	int totalTraceLength = streamline->getLength();
	VECTOR3* totalTrace = new VECTOR3[totalTraceLength];
	streamline->getTotalTrace( totalTrace );

	// Get pointer to the feature map of the streamline
	FDL_featuremap featureMap = streamline->getFeatureMap();
	
	// This function assumes that the streamline already has a valid computed feature map
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

	delete [] totalTrace;

}// end function

void
FDL_algorithm_bcd::segment_FeatureMapBased_Iterative_SingleStreamline( FDL_streamline *streamline, int slId  )
{
	assert( streamline != NULL );
	
	// Get pointer to the featuremap of the streamline
	FDL_featuremap featureMap = streamline->getFeatureMap();
	
	// Clear segment list
	streamline->featureBasedSegmentList.clear();

	// Sort all segments of the chosen resolution pair in increasing order of start point (and then decreasing order of score)
	list<SLSegment> tempList = featureMap.hierarchicalSegmentList;
	tempList.sort( SLSegment::compare_Segment_Start_Score2 ); 
	
	// Create a shorter list containing ONLY the highest scoring segment from each start point
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
			if( SLSegment::hasOverlap( *pcIter1, *pcIter2 ) &&  (*pcIter1).score > max )
				max = (*pcIter1).score;
		}
		
		SLSegment newsegment = { slId, maxLevel, (*pcIter2).start, (*pcIter2).length, max };
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
FDL_algorithm_bcd::segment_FeatureMapBased_Iterative_SingleStreamline( FDL_streamline **streamline, int slId  )
{
	assert( (*streamline) != NULL );

	// Get pointer to the featuremap of the streamline
	FDL_featuremap featureMap = (*streamline)->getFeatureMap();

	// Clear segment list
	(*streamline)->featureBasedSegmentList.clear();

	// Get all segments of the chosen resolution pair
	list<SLSegment> tempList = featureMap.hierarchicalSegmentList;
	#ifdef DEBUG_MODE
	fprintf( stderr, "Number of segments in the feature map: %d\n", (int)featureMap.hierarchicalSegmentList.size() );
	#endif

	// Get all highest level segments
	int maxLevel = featureMap.getNumLevels() - 1;
	list<SLSegment> highestLevelSegmentsList = featureMap.filter_Segments_Level( maxLevel );
	#ifdef DEBUG_MODE
	fprintf( stderr, "Number of levels in the feature map: %d\n", featureMap.getNumLevels() );
	fprintf( stderr, "Number of segments in the highest level of the feature map: %d\n", (int)highestLevelSegmentsList.size() );
	#endif

	// Sort the highest level segments based on start point
	highestLevelSegmentsList.sort( SLSegment::compare_Segment_Start_Score2 );

	#ifdef DEBUG_MODE
	fprintf( stderr, "printing highest scoring segment at each starting point of the feature map ...\n" );
	Segment::printSegmentList( highScoreSegmentsList );
	#endif

	//cout << "6" << endl;
	// For each highest level segment, find the optimal score
	for( list<SLSegment>::const_iterator pcIter1 = highestLevelSegmentsList.begin();
		 pcIter1 != highestLevelSegmentsList.end();
		 pcIter1++ )
	{
		// Find the maximum value at that start point
		float max = (*pcIter1).score;
		for( list<SLSegment>::const_iterator pcIter2 = tempList.begin();
			 pcIter2 != tempList.end();
			 pcIter2++ )
		{
			if( SLSegment::hasOverlap( *pcIter1, *pcIter2 ) &&  (*pcIter2).score > max )
				max = (*pcIter2).score;
		}

		SLSegment newsegment = { slId, maxLevel, (*pcIter1).start, (*pcIter1).length, max };
		(*streamline)->featureBasedSegmentList.push_back( newsegment );

	}// end for

	#ifdef DEBUG_MODE
	fprintf( stderr, "printing results of segmentation based on feature map ...\n" );
	Segment::printSegmentList( (*streamline)->featureBasedSegmentList );
	#endif

}// end function

void
FDL_algorithm_bcd::segment_FeatureMapBased_MultiStreamline( FDL_streamline *streamlineList,
																	int nStreamline,
																	int minLen, double joinThreshold )
{
	assert( streamlineList != NULL );

	for( int iS=0; iS<nStreamline; iS++ )
	{
		//??
		//if( iS == 272 || iS == 732 )
		//	continue;
		//??
		//if( streamlineList[iS].getFeatureMap().nLevel == 0 )
		//	this->compute_FeatureMap_SingleStreamline( streamlineList+iS, minLen, iS );
		
		// Clear any segmentation the streamline may already have
		streamlineList[iS].clearSegmentation();

		FDL_streamline* ptr = (streamlineList + iS);
		segment_FeatureMapBased_Iterative_SingleStreamline( &ptr, iS );

		joinSegments_SingleStreamline( (streamlineList + iS), iS, joinThreshold );

		#ifdef DEBUG_MODE
		fprintf( stderr, "Segmented %d-th streamline.\n", iS );
		#endif
	}

}// end function

void 
FDL_algorithm_bcd::joinSegments_SingleStreamline( FDL_streamline *streamline, int slId,
														 double joinThreshold )
{
	// Get the trace and trace length of the streamline
	int traceLength = streamline->getLength();
	VECTOR3* totalTrace = new VECTOR3[traceLength];
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
	int stretchStart = 0, stretchLen = 1;
	double stretchScore = 0, tmpScore = 0;

	for(; pcIter != tempList.end(); pcIter++ )
	{
		// Compare with last segment
		// If the values are similar, continue to be in the stretch
		// Otherwise break and create the next segment
		if( abs( lastSegment.score - (*pcIter).score ) >= joinThreshold )
		{
			// Create a segment upto this point
			stretchScore = boxdimComputer->computeBCD_SingleStreamline_SingleWindow( &totalTrace[stretchStart], stretchLen );

#if defined( _WIN32 ) || defined( _WIN64 )
			nextSegment.slId = slId;
			nextSegment.level = boxdimComputer->resType;
			nextSegment.start = stretchStart;
			nextSegment.length = stretchLen;
			nextSegment.score = stretchScore;
#else
			nextSegment = { slId, 0, stretchStart, stretchLen, stretchScore };
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
	stretchScore = boxdimComputer->computeBCD_SingleStreamline_SingleWindow( &totalTrace[stretchStart], stretchLen );

#if defined( _WIN32 ) || defined( _WIN64 )
	nextSegment.slId = slId;
	nextSegment.level = boxdimComputer->resType;
	nextSegment.start = stretchStart;
	nextSegment.length = stretchLen;
	nextSegment.score = stretchScore;
#else
	nextSegment = { slId, 0, stretchStart, stretchLen, stretchScore };
#endif

	streamline->joinedFeatureBasedSegmentList.push_back( nextSegment );

	#ifdef DEBUG_MODE
	fprintf( stderr, "printing results of joining segments ...\n" );
	fprintf( stderr, "Finally joined segments: %d\n\n", streamline->joinedFeatureBasedSegmentList.size() );
	Segment::printSegmentList( streamline->compressedFeatureBasedSegmentList );
	#endif

	delete [] totalTrace;

}// end function

int
FDL_algorithm_bcd::getNextSegment( FDL_streamline *streamline,
								   list<SLSegment> tempList,
								   SLSegment segment, int totalLength )
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
		SLSegment firstPart = { segment.slId,
								segment.start, (*pcIter).start - segment.start + 1,
								boxdimComputer->computeBCD_SingleStreamline_SingleWindow( &totalTrace[segment.start],
																						  (*pcIter).start - segment.start + 1 )
							  };
	
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

	delete [] totalTrace;

	// Return the last point of this segment
	return nextstart;
	
}// end function

list<SLSegment>
FDL_algorithm_bcd::sortSegmentsBySize( FDL_streamline *streamlineList, int nstreamline,
									   	    double** bbox, double** unitbox,
									   	    double** scoreList )
{
	list<SLSegment> allFeatureSegmentList;

	// Collect all feature segments from different streamlines into a list
	for( int iS=0; iS<nstreamline; iS++ )
	{
		allFeatureSegmentList.insert( allFeatureSegmentList.end(),
									  streamlineList[iS].joinedFeatureBasedSegmentList.begin(),
									  streamlineList[iS].joinedFeatureBasedSegmentList.end() );
	}// end for
	fprintf( stderr, "Total number of segments identified from all streamlines: %d\n", allFeatureSegmentList.size() );

	// Sort segments by score
	allFeatureSegmentList.sort( SLSegment::compare_Segment_Score2_Length2 );

	// Display first few high scoring segments
	list<SLSegment>::const_iterator pIter = allFeatureSegmentList.begin();
	int nFeaturesToBeDisplayed = min( 50, (int)allFeatureSegmentList.size() );

	// Allocate memory for bounding box and unit box of each segment, also for
	// other properties to be computed per segment
	(*bbox) = new double[nFeaturesToBeDisplayed*6];
	(*unitbox) = new double[nFeaturesToBeDisplayed*6];
	(*scoreList) = new double[nFeaturesToBeDisplayed*7];

	float featureBbox[6];
	float unitBbox[6];
	float featureVolume, arcLength, volumeovern, volumeoverlength, sparsity;
	float avgFeatureVolume = 0;
	float avgFeatureLength = 0;
	float avgFeatureVolOverN = 0;
	float avgFeatureVolOverLength = 0;
	float avgSparsity = 0;
	float minFeatureScore, maxFeatureScore;

	// Analyze the segments
	for( int i=0; i<nFeaturesToBeDisplayed; i++ )
	{
		// Get next feature segment
		SLSegment nextSeg = *pIter;

		if( i == 0 )
			maxFeatureScore = nextSeg.score;
		if( i == nFeaturesToBeDisplayed-1 )
			minFeatureScore = nextSeg.score;

		// Compute bounding box
		int totalTraceLength = streamlineList[nextSeg.slId].getLength();
		VECTOR3* totalTrace = new VECTOR3[totalTraceLength];
		streamlineList[nextSeg.slId].getTotalTrace( totalTrace );
		Streamline::computeBoundingBox( &totalTrace[nextSeg.start],
											nextSeg.length,
											featureBbox,
											&featureVolume );

		// Compute unit box tied to the bounding box
		unitBbox[0] = featureBbox[0];
		unitBbox[1] = featureBbox[1];
		unitBbox[2] = featureBbox[2];
		unitBbox[3] = featureBbox[0] + boxdimComputer->boxLen->x();
		unitBbox[4] = featureBbox[1] + boxdimComputer->boxLen->y();
		unitBbox[5] = featureBbox[2] + boxdimComputer->boxLen->z();

		// Store bounding box and unit box for rendering
		memcpy( (*bbox)+6*i, featureBbox, sizeof(float)*6 );
		memcpy( (*unitbox)+6*i, unitBbox, sizeof(float)*6 );

		#ifdef DEBUG_MODE
		fprintf( stderr, "Streamline ID: %d start:%d end: %d score:%g feature size: %g\n", nextSeg.slId, nextSeg.start,
																				  nextSeg.start + nextSeg.length - 1,
																				  nextSeg.score, featureSize );
		#endif

		// Compute other size related metrics
		FDL_streamline::computeArcLength( &totalTrace[nextSeg.start],
										  nextSeg.length,
										  &arcLength );
		volumeovern = featureVolume / (float)nextSeg.length;
		volumeoverlength = featureVolume / arcLength;
		FDL_streamline::computeSparsity( &totalTrace[nextSeg.start],
									     nextSeg.length,
									     &sparsity );


		avgFeatureVolume += featureVolume;
		avgFeatureLength += arcLength;
		avgFeatureVolOverN += volumeovern ;
		avgFeatureVolOverLength += volumeoverlength;
		avgSparsity += sparsity;

		// Clear
		delete [] totalTrace;

		// Go to next feature segment
		pIter++;
	}

	avgFeatureVolume /= (float)nFeaturesToBeDisplayed;
	avgFeatureLength /= (float)nFeaturesToBeDisplayed;
	avgFeatureVolOverN /= (float)nFeaturesToBeDisplayed;
	avgFeatureVolOverLength /= (float)nFeaturesToBeDisplayed;
	avgSparsity /= (float)nFeaturesToBeDisplayed;

	fprintf( stderr, "Feature score range: %g, %g\n", minFeatureScore, maxFeatureScore );
	fprintf( stderr, "Average feature volume, arc length, (volume/N), (volume/arc length), sparsity: %g, %g, %g, %g, %g\n",
			avgFeatureVolume, avgFeatureLength, avgFeatureVolOverN, avgFeatureVolOverLength, avgSparsity );


	//int resT[5] = { -2, -1, 0, 1, 2 };
	int resT[7] = { 0, 1, 2, 3, 4, 5, 6 };
	int nr = 7;

	pIter = allFeatureSegmentList.begin();
	int id=0;
	for( int i=0; i<nFeaturesToBeDisplayed; i++ )
	{
		SLSegment nextSeg = *pIter;

		// Get the trace and trace length of the streamline
		int totalTraceLength = streamlineList[nextSeg.slId].getLength();
		VECTOR3* totalTrace = new VECTOR3[totalTraceLength];
		streamlineList[nextSeg.slId].getTotalTrace( totalTrace );

		for( int j=0; j<nr; j++ )
		{
			//boxdimComputer->setResolutions( resT[j] );
			*((*scoreList)+id) = boxdimComputer->computeBCD_SingleStreamline_SingleWindow( &totalTrace[nextSeg.start],
				 																	  nextSeg.length );
			id++;
		}

		#ifdef DEBUG_MODE
		fprintf( stderr, "%d, %g, %g, %g, %g, %g\n", nextSeg.length,
											scoreList[i][0], scoreList[i][1],
											scoreList[i][2], scoreList[i][3],
											scoreList[i][4], scoreList[i][5],
											scoreList[i][6] );
		#endif

		// Clear
		delete [] totalTrace;

		pIter++;
	}


	// Return first few segments (the high scoring ones)
	list<SLSegment> highScoringSegmentList;
	pIter = allFeatureSegmentList.begin();
	for( int i=0; i<nFeaturesToBeDisplayed; i++ )
	{
		SLSegment nextSeg = *pIter;
		highScoringSegmentList.push_back( nextSeg );
		pIter++;

	}

	// Clear temporary resources
	allFeatureSegmentList.clear();

	return highScoringSegmentList;

}// end function

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//															File I/O functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
FDL_algorithm_bcd::save_SegmentationResult( FDL_streamline *streamlineList, int nstreamline, const char *featureFileName  )
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
			float bcd = (*pIter).score;
			
			// Write start point
			fwrite( &start, 1, sizeof(int), featureFile );
			fwrite( &length, 1, sizeof(int), featureFile );
			fwrite( &bcd, 1, sizeof(float), featureFile );
		}
	}
	
	// Close file
	fclose( featureFile );
	
}// end function



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//						Functions related to flow field subdivision and computation of feature vectors from there (NOT needed now)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
FDL_algorithm_bcd::subdivide_FlowField( VECTOR3 *min, VECTOR3 *max, int npartx, int nparty, int npartz  )
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
FDL_algorithm_bcd::computeFeatureList_AllBox( FDL_streamline *streamlinelist, int nstreamline )
{
	for( int iBox=0; iBox<nTotalBox; iBox++ )
		computeFeatureList_SingleBox( streamlinelist, nstreamline, iBox );
	fprintf( stderr, "List contains %d features.\n", partitionwiseFeatureList.size() );
	
}// end function

void
FDL_algorithm_bcd::computeFeatureList_SingleBox( FDL_streamline *streamlineList,
														int nstreamline, int boxid )
{
	int totalTraceLength = 0;
	double score = 0;
	VECTOR3* totalTrace = NULL;
	
	// Iterate through all the streamlines
	for( int iS = 0; iS<nstreamline; iS++ )
	{
		// Get the trace and trace length of the streamline
		totalTraceLength = streamlineList[iS].getLength();
		totalTrace = new VECTOR3[totalTraceLength];
		streamlineList[iS].getTotalTrace( totalTrace );

		// Compute the fractal dimensnion of this streamline with respect to this axis aligned box
		//score = boxdimComputer->computeBCD_SingleStreamline_InAxisAlignedBox( totalTrace,
		//																	  totalTraceLength,
		//																	  boxBoundary + boxid*6, boxBoundary + boxid*6+3 );
		//fprintf( stderr, "FD of Streamline %d in %d-th box: %f\n", iS, boxid, score );
		
		// Store FD as a feature if this streamline is present 
		if( score > 1.0 )
		{
			struct feature_in_box nextFeature = { boxid, iS, score };
			partitionwiseFeatureList.push_back( nextFeature ); 
		}			

		// Clear
		delete [] totalTrace;

	}// end for

}// end function

double*
FDL_algorithm_bcd::computeFeatureList_AllBox2( FDL_streamline *streamlinelist, int nstreamline )
{
	if( scoreList == NULL )	scoreList = new double[nTotalBox];

	for( int iBox=0; iBox<nTotalBox; iBox++ )
	{
		scoreList[iBox] = computeFeatureList_SingleBox2( streamlinelist, nstreamline, iBox );
		//cout << iBox << " " << scoreList[iBox] << endl;
	}	

	return scoreList;

}// end function

double
FDL_algorithm_bcd::computeFeatureList_SingleBox2( FDL_streamline *streamlineList, int nstreamline, int boxid )
{
	int totalTraceLength = 0;
	double score = 0;
	VECTOR3* totalTrace = NULL;
	//fprintf( stderr, "Box id: %d\n", boxid );
	//fprintf( stderr, "Number of streamlines: %d\n", nstreamline );
	
	// Iterate through all the streamlines
	double maxScore = -1;
	for( int iS = 0; iS<nstreamline; iS++ )
	{
		//fprintf( stderr, "Streamline ID: %d\n", iS );
		
		// Get the trace and trace length of the streamline
		totalTraceLength = streamlineList[iS].getLength();
		totalTrace = new VECTOR3[totalTraceLength];
		streamlineList[iS].getTotalTrace( totalTrace );

		// Compute the fractal dimensnion of this streamline with respect to this axis aligned box
		//score = boxdimComputer->computeBCD_SingleStreamline_InAxisAlignedBox( totalTrace,
		//																	  totalTraceLength,
		//																	  boxBoundary + boxid*6, boxBoundary + boxid*6+3 );
		//fprintf( stderr, "FD of Streamline %d in %d-th box: %f\n", iS, boxid, score );
		
		// Store FD as a feature if this streamline is present 
		if( score > maxScore )
			maxScore = score;

		// Clear
		delete [] totalTrace;

	}// end for

	return maxScore;

}// end function

/*
void
FDL_algorithm_bcd::computeFeatureVectors_AllBox_AllStreamline( FDL_streamline *streamlinelist, int nstreamline, float *featureVector )
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
		featureList.sort( FDL_algorithm_bcd::comparator_feature2 );	

		// Copy back
		for( int iB=0; iB<nTotalBox; iB++ )
		{
			*(featureVector + iS*nTotalBox + iB) = featureList.front();
			featureList.pop_front();
		}
	}
	
	
}// end function

void
FDL_algorithm_bcd::computeFeatureVectors_AllBox_SingleStreamline( FDL_streamline *streamline, float *featureVector  )
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
FDL_algorithm_bcd::saveFeatureVectors( float *featureVector, int nstreamline, const char *featureFileName  )
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
FDL_algorithm_bcd::projectFeatures_AllStreamlines( int nstreamline, int K, float *featureSpan, float *featureScore )
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
FDL_algorithm_bcd::listFeatures_SingleStreamline( int streamlineID )
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
FDL_algorithm_bcd::comparator_feature( struct feature_in_box first, struct feature_in_box second )
{	
	if( first.feature_score >= second.feature_score )	return true;
	return false;
	
}// end function

bool
FDL_algorithm_bcd::comparator_feature2( float first, float second )
{	
	if( first >= second )	return true;
	return false;
	
}// end function


list<struct feature_in_box> 
FDL_algorithm_bcd::listTopKFeatures_SingleStreamline( list<struct feature_in_box> featureList_currentStreamline, int K )
{
	// Get total number of features 
	int nFeature = featureList_currentStreamline.size();
	
	// Sort features in decending order
	featureList_currentStreamline.sort( FDL_algorithm_bcd::comparator_feature );
	
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
FDL_algorithm_bcd::compute_TopKFeaturesProjection_SingleStreamline(  list<struct feature_in_box>  featureList_topK_currentStreamline, int K )
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
FDL_algorithm_bcd::compute_FeatureSpan_SingleStreamline(  list<struct feature_in_box>  featureList_topK_currentStreamline, int K )
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

