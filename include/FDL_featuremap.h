/*
 * FDL_featuremap.h
 *
 *  Created on: Feb 19, 2011
 *      Author: abon
 */

#ifndef FDL_FEATUREMAP_H_
#define FDL_FEATUREMAP_H_

#include <list>
#include "FDL_header.h"
#include "FDL_util.h"

class FDL_featuremap
{
public:

	int scoreType;
	int nLevel;												// Number of levels (should be same for any resolution pair)
	int nSegment;											// Total number of segments in all segments (should be same for any resolution pair)
	list<struct Segment> hierarchicalSegmentList;			// List of segments

public:

	// Constructor
	FDL_featuremap();
	FDL_featuremap( const FDL_featuremap& that );
	FDL_featuremap& operator= ( const FDL_featuremap& that );

	// Destructor
	~FDL_featuremap();

	// Get-set functions
	void setNumLevels( int nL );
	void setNumSegments( int nS );
	int getNumLevels();
	int getNumSegments();

	// Add segments 
	void add_Segment( struct Segment nextSegment );

	// Sort segments based on different criteria
	void sort_Segments_Level_Start();
	void sort_Segments_Start_Score2();
	
	// Filter segments based on different criteria
	list<struct Segment> filter_Segments_Level( int level );
	list<struct Segment> filter_Segments_Start( int start );
	struct Segment filter_Segment_Start_Length( int startPoint, int length );
	struct Segment filter_TopSegment_Level( int level );
	struct Segment filter_TopSegment_Start( int startPoint );

	void printFeatureMap();

};

#endif /* FDL_featuremap_H_ */
