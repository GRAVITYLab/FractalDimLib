/*
 * FDL_featuremap.cpp
 *
 *  Created on: Feb 19, 2011
 *      Author: abon
 */

#include "FDL_featuremap.h"

FDL_featuremap::FDL_featuremap()
{
	scoreType = 0;
	nLevel = 0;
	nSegment = 0;

}// constructor

FDL_featuremap::FDL_featuremap( const FDL_featuremap& that )
{
	this->scoreType = that.scoreType;
	this->nLevel = that.nLevel;
	this->nSegment = that.nSegment;

	if( that.hierarchicalSegmentList.size() > 0 )
		this->hierarchicalSegmentList.insert( this->hierarchicalSegmentList.end(),
											  that.hierarchicalSegmentList.begin(),
											  that.hierarchicalSegmentList.end() );

}// constructor

FDL_featuremap&
FDL_featuremap::operator= ( const FDL_featuremap& that )
{
	if (this != &that ) // protect against invalid self-assignment
	{
		this->scoreType = that.scoreType;
		this->nLevel = that.nLevel;
		this->nSegment = that.nSegment;

		if( that.hierarchicalSegmentList.size() > 0 )
			this->hierarchicalSegmentList.insert( this->hierarchicalSegmentList.end(),
												  that.hierarchicalSegmentList.begin(),
												  that.hierarchicalSegmentList.end() );


	}
	// by convention, always return *this
	return *this;

}// constructor

FDL_featuremap::~FDL_featuremap()
{
	hierarchicalSegmentList.clear();

}// destructor

void
FDL_featuremap::setNumLevels( int nL )
{
	nLevel = nL;
}// end function

void
FDL_featuremap::setNumSegments( int nS )
{
	nSegment = nS;
}// end function

int
FDL_featuremap::getNumLevels()
{
	return nLevel;
}// end function

int
FDL_featuremap::getNumSegments()
{
	return nSegment;
}// end function

void
FDL_featuremap::add_Segment( struct Segment nextSegment )
{
	hierarchicalSegmentList.push_back( nextSegment );
}// end function

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//													Sorting related functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FDL_featuremap::sort_Segments_Level_Start()
{
	hierarchicalSegmentList.sort( Segment::compare_Segment_Level_Start );
}// end function

void
FDL_featuremap::sort_Segments_Start_Score2()
{
	hierarchicalSegmentList.sort( Segment::compare_Segment_Start_Score2 );
}// end function

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//													Filtering related functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
list<struct Segment>
FDL_featuremap::filter_Segments_Level( int level )
{
	list<struct Segment> filteredList;
	list<struct Segment>::const_iterator pIter;

	for( pIter=hierarchicalSegmentList.begin();
		 pIter != hierarchicalSegmentList.end();
		 pIter++ )
	{
		if( (*pIter).level > level)
			break;
		if( (*pIter).level == level)
			filteredList.push_back( *pIter );
	}
	
	return filteredList;
	
}// end function

list<struct Segment>
FDL_featuremap::filter_Segments_Start( int start )
{
	list<struct Segment> filteredList;
	list<struct Segment>::const_iterator pIter;
	
	for( pIter = hierarchicalSegmentList.begin();
		 pIter != hierarchicalSegmentList.end();
		 pIter++ )
	{
		if( (*pIter).start !=  start )
			continue;
		filteredList.push_back( *pIter );
	}
	
	return filteredList;
	
}// end function

struct Segment
FDL_featuremap::filter_Segment_Start_Length( int start, int length )
{
	list<struct Segment>::const_iterator pIter;

	for( pIter = hierarchicalSegmentList.begin();
		 pIter != hierarchicalSegmentList.end();
		 pIter++ )
	{
		if( (*pIter).start ==  start &&  (*pIter).length ==  length )
			return (*pIter);
	}
	
	struct Segment dummy = { -1, -1, -1, -1 };
	return dummy;
	
}// end function

/*
 * Get all segments of the selected level
 */
struct Segment
FDL_featuremap::filter_TopSegment_Level( int level )
{
	list<struct Segment> filteredList = filter_Segments_Level( level );
	
	// Sort all segments of this level in decreasing order of score
	filteredList.sort( Segment::compare_Segment_Score2_Length2 );
	
	return filteredList.front();
	
}// end function

/*
 * Get all segments having a given start point
 */
struct Segment
FDL_featuremap::filter_TopSegment_Start( int start )
{
	list<struct Segment> filteredList = filter_Segments_Start( start );
	
	// Sort all segment of the filtered segments based on score
	filteredList.sort( Segment::compare_Segment_Start_Score2 );
	
	return filteredList.front();
	
}// end function

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//													Printing/Saving related functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
FDL_featuremap::printFeatureMap()
{
	list<struct Segment>::const_iterator pIter;

	fprintf( stderr, "Feature map:\n" );
	for( pIter = hierarchicalSegmentList.begin();
		 pIter != hierarchicalSegmentList.end();
		 pIter++ )
		fprintf( stderr, "<level: %d start: %d end: %d score: %f>\n", (*pIter).level,
																	   (*pIter).start,
																	   (*pIter).start + (*pIter).length-1,
																	   (*pIter).score );
	
}// end function
