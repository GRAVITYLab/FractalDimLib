/*
 * A Collection of all includes and definitions.
 * Adopted from OSU Flow Vis Library
 * Created: Han-Wei Shen, Liya Li
 * The Ohio State University
 * Date: 06/2005
 *
 */

#ifndef _FDL_HEADER_H_
#define _FDL_HEADER_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
//#include <cmath>
#include <cassert>
//#include <float.h>
#include <ctime>
#include <vector>
#include <string>
#include <list>
#include <set>
#include <queue>
#include <algorithm>
#include <iterator>

#include "OSUFlow.h"

using namespace std;

// Definitions of some constant parameters
#define CLOCKS_PER_MS CLOCKS_PER_SEC/1000.0f

//#define NRES 5
#define NRES 25

//#define DEBUG_MODE

struct Segment
{
	int slId;
	int level;
	int start;
	int length;
	double score;

	static bool
	hasOverlap( struct Segment seg1, struct Segment seg2  )
	{
		if( seg1.start < seg2.start && ( seg1.start + seg1.length-1) <= seg2.start ) return false;
		if( seg2.start < seg1.start && ( seg2.start + seg2.length-1) <= seg1.start ) return false;
		return true;
	}// end function

	static bool
	hasOverlap( struct Segment seg, list<struct Segment> segList )
	{
		bool ret = false;

		for( list<struct Segment>::const_iterator pIter = segList.begin(); pIter != segList.end(); pIter++ )
		{
			if( Segment::hasOverlap( seg, (*pIter) ) )
				return true;
		}

		return ret;
	}// end function

	static bool
	isEqual( struct Segment seg1, struct Segment seg2 )
	{
		return ( seg1.level == seg2.level ) && ( seg1.start == seg2.start ) && ( seg1.length == seg2.length ) && ( seg1.score == seg2.score );
	}

	static bool
	compare_Segment_Level_Start( struct Segment seg1, struct Segment seg2 )
	{
		if( seg1.level < seg2.level ) return true;
		if( seg1.level > seg2.level ) return false;
		if( seg1.level == seg2.level ) return ( seg1.start < seg2.start );
		//return false;
	}// end function

	static bool
	compare_Segment_Score2_Length2( struct Segment seg1, struct Segment seg2 )
	{
		if( seg1.score > seg2.score ) return true;
		if( seg1.score == seg2.score ) return ( seg1.length >= seg2.length );
		if( seg1.score < seg2.score ) return false;
		//return false;

	}// end function

	static bool
	compare_Segment_Start_Score2( struct Segment seg1, struct Segment seg2 )
	{
		assert( seg1.score >= 0 && seg2.score >= 0 ); 
		
		if( seg1.start < seg2.start ) return true;
		if( seg1.start == seg2.start ) return ( seg1.score > seg2.score );
		if( seg1.start > seg2.start ) return false;
		//return false;
	}// end function

	static bool
	compare_Segment_Start_Length2( struct Segment seg1, struct Segment seg2 )
	{
		if( seg1.start < seg2.start ) return true;
		if( seg1.start == seg2.start ) return ( seg1.length >= seg2.length );
		if( seg1.start > seg2.start ) return false;

	}// end function

	static void
	printSegmentList( list<struct Segment> segList )
	{
		for( list<struct Segment>::const_iterator pcIter = segList.begin(); pcIter != segList.end(); pcIter++ )
			fprintf( stderr, "start:%d end: %d score:%g \n",
					 	 	  (*pcIter).start,
					 	 	  (*pcIter).start + (*pcIter).length - 1,
					 	 	  (*pcIter).score );
	}// end function

};

typedef struct Segment SLSegment;

struct feature_in_box
{
	int box_id;
	int streamline_id;
	float feature_score;
};

#endif
