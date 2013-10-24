/**
 * FDL_streamline.cpp
 *
 *  Created on: Dec 15, 2010
 *      Author: abon
 */

#include "FDL_streamline.h"

FDL_streamline::FDL_streamline()
{
	this->fwdTraceLength = 0;
	this->backTraceLength = 0;
	this->totalTraceLength = 0;
	this->boxCountDim = 0;

	// NULL initialize trace
	fwdTrace = NULL;
	backTrace = NULL;
	totalTrace = NULL;
	inoutList = NULL;

}// end constructor

FDL_streamline::FDL_streamline( const FDL_streamline& that )
{
	this->fwdTraceLength = that.fwdTraceLength;
	this->backTraceLength = that.backTraceLength;
	this->totalTraceLength = that.totalTraceLength;
	this->boxCountDim = that.boxCountDim;
	this->featureMap = that.featureMap;
	memcpy( this->boundingBox, that.boundingBox, sizeof(float)*6 );

	fwdTrace = NULL;
	backTrace = NULL;
	totalTrace = NULL;
	inoutList = NULL;

	if( that.fwdTrace != NULL )
	{
		this->fwdTrace = new VECTOR3[fwdTraceLength];
		memcpy( this->fwdTrace, that.fwdTrace, sizeof(VECTOR3)*fwdTraceLength );
	}
	if( that.backTrace != NULL )
	{
		this->backTrace = new VECTOR3[backTraceLength];
		memcpy( this->backTrace, that.backTrace, sizeof(VECTOR3)*backTraceLength );
	}
	if( that.totalTrace != NULL )
	{
		this->totalTrace = new VECTOR3[totalTraceLength];
		memcpy( this->totalTrace, that.totalTrace, sizeof(VECTOR3)*totalTraceLength );
	}
	if( that.inoutList != NULL )
	{
		this->inoutList = new int[totalTraceLength];
		memcpy( this->inoutList, that.inoutList, sizeof(int)*totalTraceLength );
	}

}// end constructor

FDL_streamline&
FDL_streamline::operator= ( const FDL_streamline& that )
{
	if (this != &that ) // protect against invalid self-assignment
	{
		this->fwdTraceLength = that.fwdTraceLength;
		this->backTraceLength = that.backTraceLength;
		this->totalTraceLength = that.totalTraceLength;
		this->boxCountDim = that.boxCountDim;
		this->featureMap = that.featureMap;
		memcpy( this->boundingBox, that.boundingBox, sizeof(float)*6 );

		fwdTrace = NULL;
		backTrace = NULL;
		totalTrace = NULL;
		inoutList = NULL;

		if( that.fwdTrace != NULL )
		{
			this->fwdTrace = new VECTOR3[fwdTraceLength];
			memcpy( this->fwdTrace, that.fwdTrace, sizeof(VECTOR3)*fwdTraceLength );
		}
		if( that.backTrace != NULL )
		{
			this->backTrace = new VECTOR3[backTraceLength];
			memcpy( this->backTrace, that.backTrace, sizeof(VECTOR3)*backTraceLength );
		}
		if( that.totalTrace != NULL )
		{
			this->totalTrace = new VECTOR3[totalTraceLength];
			memcpy( this->totalTrace, that.totalTrace, sizeof(VECTOR3)*totalTraceLength );
		}
		if( that.inoutList != NULL )
		{
			this->inoutList = new int[totalTraceLength];
			memcpy( this->inoutList, that.inoutList, sizeof(int)*totalTraceLength );
		}
	}
	// by convention, always return *this
	return *this;
}


FDL_streamline::~FDL_streamline()
{
	if( fwdTrace != NULL ) delete [] fwdTrace;
	if( backTrace != NULL ) delete [] backTrace;
	if( totalTrace != NULL ) delete [] totalTrace;
	if( inoutList != NULL ) delete [] inoutList;

}// end destructor

void
FDL_streamline::setSeed( VECTOR3 seed )
{
	this->seed.Set( seed[0], seed[1], seed[2] );
}// end function

void
FDL_streamline::setTraceData( VECTOR3* fwdtrace, VECTOR3* backtrace )
{
	this->fwdTrace = new VECTOR3[this->fwdTraceLength];
	for( int i=0; i<this->fwdTraceLength; i++ )
		this->fwdTrace[i] = fwdtrace[i];

	this->backTrace = new VECTOR3[this->backTraceLength];
	for( int i=0; i<this->backTraceLength; i++ )
		this->backTrace[i] = backtrace[i];

}// end function

void
FDL_streamline::setTraceLength( int fwdtracelen, int backtracelen )
{
	this->fwdTraceLength = fwdtracelen;
	this->backTraceLength = backtracelen;
	this->totalTraceLength = this->fwdTraceLength + this->backTraceLength -	1;
}// end function

void
FDL_streamline::setTotalTraceData( VECTOR3* totaltrace  )
{
	this->totalTrace = new VECTOR3[this->totalTraceLength];
	for( int i=0; i<this->totalTraceLength; i++ )
		this->totalTrace[i] = totaltrace[i];
}// end function

void
FDL_streamline::setTotalTraceLength( int totaltracelen )
{
	this->totalTraceLength = totaltracelen;
}// end function

void
FDL_streamline::setFeatureMap( FDL_featuremap fMap )
{
	// Copy the segment information for all resolutions
	for( int iR=0; iR<NRES; iR++ )
	{
		featureMap.hierarchicalSegmentList.clear();
		featureMap.hierarchicalSegmentList = fMap.hierarchicalSegmentList;
	}
	
	// Copy total number of levels and segments
	featureMap.setNumLevels( fMap.getNumLevels() );
	featureMap.setNumSegments( fMap.getNumSegments() );

}// end function

void
FDL_streamline::clearFeatureMap()
{
	featureMap.hierarchicalSegmentList.clear();
	featureMap.nLevel = 0;
	featureMap.nSegment = 0;

}// end function

void
FDL_streamline::clearSegmentation()
{
	fixedLengthSegmentList.clear();
	featureBasedSegmentList.clear();
	joinedFeatureBasedSegmentList.clear();
}


void
FDL_streamline::computeBoundingBox()
{
	assert( totalTrace != NULL );
	assert( totalTraceLength > 0 );
	
	// Scan through all trace points	
	for( int iT=0; iT<totalTraceLength; iT++ )
	{
		if( iT == 0 )
		{
			boundingBox[0] = boundingBox[3] = totalTrace[iT].x();
			boundingBox[1] = boundingBox[4] = totalTrace[iT].y();
			boundingBox[2] = boundingBox[5] = totalTrace[iT].z();
		}
		else
		{
			if( totalTrace[iT].x() < boundingBox[0] ) boundingBox[0] = totalTrace[iT].x();
			if( totalTrace[iT].x() > boundingBox[3] ) boundingBox[3] = totalTrace[iT].x();
			if( totalTrace[iT].y() < boundingBox[1] ) boundingBox[1] = totalTrace[iT].y();
			if( totalTrace[iT].y() > boundingBox[4] ) boundingBox[4] = totalTrace[iT].y();
			if( totalTrace[iT].z() < boundingBox[2] ) boundingBox[2] = totalTrace[iT].z();
			if( totalTrace[iT].z() > boundingBox[5] ) boundingBox[5] = totalTrace[iT].z();
			
		}
	}	

}// end function

VECTOR3
FDL_streamline::getForwardTraceAt( int t )
{
	assert( t < fwdTraceLength );
	return this->fwdTrace[t];
}// end function

VECTOR3
FDL_streamline::getBackwardTraceAt( int t )
{
	assert( t < backTraceLength );
	return this->backTrace[t];
}// end function

int
FDL_streamline::getLength()
{
	return this->totalTraceLength;
}// end function

float
FDL_streamline::getBoxCountDimension()
{
	return boxCountDim;
}// end function

FDL_featuremap
FDL_streamline::getFeatureMap()
{
	return featureMap;
}// end function

void
FDL_streamline::getSegment( FDL_streamline *streamlinePart, int startPoint, int length )
{
	assert( streamlinePart != NULL );
	streamlinePart->setTotalTraceLength( length );
	streamlinePart->setTotalTraceData( this->totalTrace+startPoint );
}// end function

void
FDL_streamline::getTotalTraceData( VECTOR3* totaltrace, int totaltracelen )
{
	for( int i=0; i<this->totalTraceLength; i++ )
		 totaltrace[i] = this->totalTrace[i];

}

VECTOR3*
FDL_streamline::joinTwoTraceDirections( FDL_streamline *streamline )
{
	VECTOR3* totalTrace = new VECTOR3[streamline->totalTraceLength];

	#ifdef DEBUG_MODE
	fprintf( stderr, "Next streamline\n" );
	fprintf( stderr, "Printing back trace\n" );
	for(int i=0; i<streamline->backTraceLength; i++ )
		fprintf( stderr, "%g %g %g\n", streamline->backTrace[i][0], streamline->backTrace[i][1], streamline->backTrace[i][2] );
	fprintf( stderr, "Printing forward trace\n" );
	for(int i=0; i<streamline->fwdTraceLength; i++ )
		fprintf( stderr, "%g %g %g\n", streamline->fwdTrace[i][0], streamline->fwdTrace[i][1], streamline->fwdTrace[i][2] );
	#endif

	// Copy location information from streamline
	// Goes from -1 to -bakT and then 0 to fwdT
	for( int i=0; i<streamline->backTraceLength; i++ )
		totalTrace[i] = streamline->backTrace[streamline->backTraceLength-1-i];
	for( int i=0; i<streamline->fwdTraceLength-1; i++ )
		totalTrace[i+streamline->backTraceLength] = streamline->fwdTrace[i+1];

	return totalTrace;

}// end function

/*
void
FDL_streamline::computeBoxStreamlineIntersection( VECTOR3 *corners, MATRIX4 inv )
{
	// Initialize in-out list
	if(  this->inoutList == NULL )	this->inoutList = new int[this->totalTraceLength];
	for( int i=0; i<this->totalTraceLength; i++ )	this->inoutList[i] = 0;
		
	// Compute the inverse transformation of the box
	VECTOR3 *invTransformedCorners = new VECTOR3[8];
	for( int i=0; i<8; i++ )
	{
		invTransformedCorners[i] = this->transformPoint2( corners[i], inv );
		//cout << invTransformedCorners[i].x() << " " << invTransformedCorners[i].y() << " " << invTransformedCorners[i].z() << endl;
	}	
	
	// Go through each point along the trace and determine if it is inside the box outside
	for( int i=0; i<this->totalTraceLength; i++ )
	{
		this->inoutList[i] = this->isPointInsideBox( this->totalTrace[i], invTransformedCorners, inv );
		//if( this->inoutList[i] == 1 )	fprintf( stderr, "%d: in\n", i );
	}
	
	delete [] invTransformedCorners;
	
}// end function

int
FDL_streamline::isPointInsideBox( VECTOR3 point, VECTOR3 *corners, MATRIX4 inv )
{
	// Compute inverse transformation of the point
	VECTOR3 invPoint = this->transformPoint2( point, inv );
	//cout << invPoint[0] << " " << invPoint[1] << " " << invPoint[2] << endl;
	
	// Check if transformed point is inside transformed box or not
	if( invPoint[0] < corners[0].x() ||  invPoint[0] > corners[1].x() )	return 0;
	if( invPoint[1] < corners[0].y() ||  invPoint[1] > corners[3].y() )	return 0;
	if( invPoint[2] < corners[0].z() ||  invPoint[2] > corners[4].z() )	return 0;
	
	return 1;
	
}// end function
*/

VECTOR3
FDL_streamline::transformPoint2( VECTOR3 point, MATRIX4 transformMat )
{
	VECTOR4 temp( point[0], point[1], point[2], 1.0f );

	temp =  transformMat * temp;
	
	return VECTOR3( temp[0], temp[1], temp[2] );		
	
}// end function

list<VECTOR4>
FDL_streamline::color_InOutBased()
{
	assert( this->inoutList != NULL );
	
	// Initialize list of colors
	list<VECTOR4> streamlineColorList;
	VECTOR4 inColor( ( rand() % 1000 ) / (float)1000.0f, ( rand() % 1000 ) / (float)1000.0f, ( rand() % 1000 ) / (float)1000.0f, 1.0f );
	VECTOR4 outColor( 0.5f, 0.5f, 0.5f, 1.0f );
	
	// Go through each point along the trace and determine color
	for( int i=0; i<this->totalTraceLength; i++ )
	{
		if( this->inoutList[i] == 1 )	streamlineColorList.push_back(  inColor );
		else						streamlineColorList.push_back( outColor );		
	}	
	
	return streamlineColorList;
	
}// end function


