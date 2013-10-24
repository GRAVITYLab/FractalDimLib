/*
 * FDL_streamline.h
 *
 *  Created on: Dec 15, 2010
 *      Author: abon
 */

#ifndef FDL_STREAMLINE_H_
#define FDL_STREAMLINE_H_

#include "Streamline.h"
#include "FDL_header.h"
#include "FDL_featuremap.h"

class FDL_streamline : public Streamline
{
private:

	double boxCountDim;
	int *inoutList;

public:

	list<VECTOR4> color_InOutBased();

	// Hierarchical fractal dimension associated to the streamline
	FDL_featuremap featureMap;
	
	// Different possible srgmentations of the streamline
	list<struct Segment> fixedLengthSegmentList;
	list<struct Segment> featureBasedSegmentList;
	list<struct Segment> joinedFeatureBasedSegmentList;
	float minSegmentScore[3];
	float maxSegmentScore[3];

	// Highlighted features from this streamline
	list<struct Segment> highlightedSegmentList;

	FDL_streamline();
	FDL_streamline( const FDL_streamline& that );
	FDL_streamline& operator= ( const FDL_streamline& that );
	~FDL_streamline();

	// Get functions
	VECTOR3 getForwardTraceAt( int t );
	VECTOR3 getBackwardTraceAt( int t );
	int getLength();
	float getBoxCountDimension();
	FDL_featuremap getFeatureMap();

	void getTotalTraceData( VECTOR3* totaltrace, int totaltracelen );

	// Set functions	
	void setSeed( VECTOR3 seed );
	void setTraceData( VECTOR3* curfwdtrace, VECTOR3* curbacktrace );
	void setTraceLength( int fwdtracelen, int backtracelen );
	void setTotalTraceData( VECTOR3* totaltrace  );
	void setTotalTraceLength( int totaltracelen );
	void setFeatureMap( FDL_featuremap fMap );

	void clearFeatureMap();
	void clearSegmentation();

	// Create a smaller streamline from its part
	void getSegment( FDL_streamline *streamlinePart, int startPoint, int length );

	// Geometric operations
	void computeBoundingBox();
	static VECTOR3* joinTwoTraceDirections( FDL_streamline *streamline );

	/*
	void computeBoxStreamlineIntersection( VECTOR3 *corners, MATRIX4 inv );
	int isPointInsideBox( VECTOR3 point, VECTOR3 *corners, MATRIX4 inv );
	*/
	
	VECTOR3 transformPoint2( VECTOR3 point, MATRIX4 transformMat );
};


#endif /* FDL_STREAMLINE_H_ */
