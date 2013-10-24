/**
 * FDL_osuflow.cpp
 *
 *  Created on: Jan 26, 2010
 *  Author: Abon
 */

#include "FDL_osuflow.h"

FDL_osuflow::FDL_osuflow( OSUFlow *osuflow )
{
	// Get pointer to osuflow
	this->osuflow = osuflow;

	// Get global bounds of the vector field
	osuflow->GetGlobalBounds( this->dataLow, this->dataHigh );
	for( int i=0; i<3; i++ )
		dataDim[i] = (int) ( dataHigh[i] - dataLow[i] + 1 );

	// NULL initialize rest of variables
	this->seeds = NULL;
	this->streamlineLengthList = NULL;
	this->seedFile = NULL;
	this->streamlineFile = NULL;

}// end constructor

void FDL_osuflow::generateSeedAtEachVoxel()
{
	#ifdef DEBUG_MODE
	printf( "Generating seed at each vertex of the vector field\n" );
	#endif

	// Compute and store the total number of seeds generated
	this->nSeed = dataDim[0] * dataDim[1] * dataDim[2];
	#ifdef DEBUG_MODE
	printf( "Total number of seeds to be generated: %d\n", this->nSeed );
	#endif

	// Allocate memory for seeds
	if( this->seeds != NULL )	delete [] this->seeds;
	this->seeds = new VECTOR3[this->nSeed];

	// Generate the seed points at each grid vertex
	int seedIndex = 0;
	for( float z=dataLow[2]; z<=this->dataHigh[2]; z++ )
	{
		for( float y=dataLow[1]; y<this->dataHigh[1]; y++ )
		{
			for( float x=dataLow[0]; x<this->dataHigh[0]; x++ )
			{
				// Set Seed
				this->seeds[seedIndex].Set( x, y, z );

				// Move to the next seed
				seedIndex ++;
			}
		}
	}

	// Print the seed points (optional)
	#ifdef DEBUG_MODE
	printf( "Seeds generated\n" );
	printf( "Printing seeds ...\n" );
	for (int i=0; i<this->nSeed; i++)
	    printf("Seed %d: [%f %f %f]\n", i, this->seeds[i][0], this->seeds[i][1], this->seeds[i][2]);
	#endif

}// end function

void FDL_osuflow::generateSeedAtEachVoxel( float* low, float* high )
{
	#ifdef DEBUG_MODE
	printf( "Generating seed at each grid vertex within a given bounding box\n" );
	#endif

	// Compute and store the total number of seeds generated in this part
	float *tmp = FDL_util<float>::subtractArrays( high, low, 3 );
	FDL_util<float>::addArrayScalar( tmp, 1, 3 );
	this->nSeed = (int) FDL_util<float>::prod( tmp, 3 );

	#ifdef DEBUG_MODE
	printf( "Total number of seeds to be generated: %d\n", this->nSeed );
	#endif

	// Allocate memory for seeds
	if( this->seeds != NULL )	delete [] this->seeds;
	this->seeds = new VECTOR3[this->nSeed];

	// Generate the seed points at each grid vertex
	int seedIndex = 0;
	for( float z=low[2]; z<=high[2]; z++ )
	{
		for( float y=low[1]; y<=high[1]; y++ )
		{
			for( float x=low[0]; x<=high[0]; x++ )
			{
				this->seeds[seedIndex].Set( x, y, z );

				// Move to the next seed
				seedIndex ++;
			}
		}
	}

	// Print the seed points (optional)
	#ifdef DEBUG_MODE
	printf( "Seeds generated within bounding box\n" );
	printf( "Printing seeds ...\n" );
	for (int i=0; i<this->nSeed; i++)
		printf("Seed %d: [%f %f %f]\n", i, this->seeds[i][0], this->seeds[i][1], this->seeds[i][2]);
	#endif

	delete [] tmp;

}// end function

void FDL_osuflow::generateSeedAlongLine( float* start, float* end, int nseed )
{
	#ifdef DEBUG_MODE
	printf( "Generating seed along a given line\n" );
	#endif

	// Check if both start and end points are within the data domain
	assert( start[0] >= dataLow[0] );
	assert( start[1] >= dataLow[1] );
	assert( start[2] >= dataLow[2] );
	assert( end[0] <= dataHigh[0] );
	assert( end[1] <= dataHigh[1] );
	assert( end[2] <= dataHigh[2] );

	// Store the total number of seeds to be generated along the line
	this->nSeed = nseed;

	#ifdef DEBUG_MODE
	printf( "Total number of seeds to be generated: %d\n", this->nSeed );
	#endif

	// Allocate memory for seeds
	if( this->seeds != NULL ) delete [] this->seeds;
	this->seeds = new VECTOR3[this->nSeed];

	// Direction
	VECTOR3 dir( end[0] - start[0], end[1] - start[1], end[2] - start[2] );

	// Increment or step size along the line
	float delt = 1.0f / (float)( this->nSeed-1 );
	#ifdef DEBUG_MODE
	printf( "Step size: %f\n", delt );
	#endif

	float t = 0;
	for( int i=0; i<this->nSeed; i++ )
	{
		this->seeds[i].Set( start[0] + t*dir[0], start[1] + t*dir[1], start[2] + t*dir[2] );

		// Increment t
		t += delt;
	}

	// Print the seed points (optional)
	#ifdef DEBUG_MODE
	printf( "Seeds generated along line\n" );
	printf( "Printing seeds ...\n" );
	for (int i=0; i<this->nSeed; i++)
		printf("Seed %d: [%f %f %f]\n", i, this->seeds[i][0], this->seeds[i][1], this->seeds[i][2]);
	#endif

}// end function

void
FDL_osuflow::generateSeedOnPlane( VECTOR3 *corners )
{
	#ifdef DEBUG_MODE
	printf( "Generating seed on a given plane\n" );
	#endif
	
	//cout <<  "0: " << corners[0].x() << " " << corners[0].y() << " " << corners[0].z() << endl;
	//cout <<  "1: " << corners[1].x() << " " << corners[1].y() << " " << corners[1].z() << endl;
	//cout <<  "2: " << corners[2].x() << " " << corners[2].y() << " " << corners[2].z() << endl;
	//cout <<  "3: " << corners[3].x() << " " << corners[3].y() << " " << corners[3].z() << endl;
	
	// Compute the span of the box
	float spanX = sqrt( ( corners[1].x() - corners[0].x() ) * ( corners[1].x() - corners[0].x() ) +   ( corners[1].y() - corners[0].y() ) * ( corners[1].y() - corners[0].y() ) + ( corners[1].z() - corners[0].z() ) * ( corners[1].z() - corners[0].z() ) );
	float spanY = sqrt( ( corners[3].x() - corners[0].x() ) * ( corners[3].x() - corners[0].x() ) +   ( corners[3].y() - corners[0].y() ) * ( corners[3].y() - corners[0].y() ) + ( corners[3].z() - corners[0].z() ) * ( corners[3].z() - corners[0].z() ) );
	
	//cout << spanX << " " << spanY << endl;
	
	// Compute the directions of movement
	VECTOR3 dir1(   corners[1].x() - corners[0].x(),  corners[1].y() - corners[0].y(),  corners[1].z() - corners[0].z() );
	VECTOR3 dir2(   corners[3].x() - corners[0].x(),  corners[3].y() - corners[0].y(),  corners[3].z() - corners[0].z() );
	
	// Compute the total number of seeds to be generated on the plane
	this->nSeed = (int) ( floor( spanX ) * floor( spanY ) ); 
	int nseed1 = (int) floor( spanX );
	int nseed2 = (int) floor( spanY );
	
	//cout << "dir1: " << dir1[0] << " " << dir1[1] << " " << dir1[2] << endl;
	//cout << "dir2: " << dir2[0] << " " << dir2[1] << " " << dir2[2] <<  endl;
	//cout << "nseeds: " << nseed1 << " " << nseed2 << endl;
		
	#ifdef DEBUG_MODE
	printf( "Total number of seeds to be generated: %d\n", this->nSeed );
	#endif

	// Allocate memory for seeds
	if( this->seeds != NULL ) delete [] this->seeds;
	this->seeds = new VECTOR3[this->nSeed];

	// Increment or step size along the plane
	float delt1 = 1.0f / (float)( nseed1-1 );
	float delt2 = 1.0f / (float)( nseed2-1 );
	#ifdef DEBUG_MODE
	//printf( "Step size: %f %f\n", delt1, delt2 );
	#endif

	int seedIndex = 0;
	float t1 = 0;
	float t2 = 0;
	VECTOR3 base;
	for( int i=0; i<nseed1; i++ )
	{
		// Move one step along first direction
		base.Set( corners[0].x() + t1 * dir1[0], corners[0].y() + t1 * dir1[1], corners[0].z() + t1 * dir1[2] );
		
		t2 = 0;
		for( int j=0; j<nseed2; j++ )
		{
			// Move one step along second direction
			this->seeds[seedIndex].Set(  base[0] + t2 * dir2[0], base[1] + t2 * dir2[1], base[2] + t2 * dir2[2] );
			
			// Move to next seed
			seedIndex++;
			
			// Increment t2
			t2 += delt2;
		}
		
		// Increment t1
		t1 += delt1;
	}

	// Print the seed points (optional)
	#ifdef DEBUG_MODE
	printf( "Seeds generated along line\n" );
	printf( "Printing seeds ...\n" );
	for (int i=0; i<this->nSeed; i++)
		printf("Seed %d: [%f %f %f]\n", i, this->seeds[i][0], this->seeds[i][1], this->seeds[i][2]);
	#endif

}// end function

void
FDL_osuflow::generateSeedInBox( VECTOR3 *corners )
{
	#ifdef DEBUG_MODE
	printf( "Generating seed on a given plane\n" );
	#endif
	
	// Compute the span of the box
	float spanX = sqrt( ( corners[1].x() - corners[0].x() ) * ( corners[1].x() - corners[0].x() ) +   ( corners[1].y() - corners[0].y() ) * ( corners[1].y() - corners[0].y() ) + ( corners[1].z() - corners[0].z() ) * ( corners[1].z() - corners[0].z() ) );
	float spanY = sqrt( ( corners[3].x() - corners[0].x() ) * ( corners[3].x() - corners[0].x() ) +   ( corners[3].y() - corners[0].y() ) * ( corners[3].y() - corners[0].y() ) + ( corners[3].z() - corners[0].z() ) * ( corners[3].z() - corners[0].z() ) );
	float spanZ = sqrt( ( corners[4].x() - corners[0].x() ) * ( corners[4].x() - corners[0].x() ) +   ( corners[4].y() - corners[0].y() ) * ( corners[4].y() - corners[0].y() ) + ( corners[4].z() - corners[0].z() ) * ( corners[4].z() - corners[0].z() ) );
	
	//cout << spanX << " " << spanY << endl;
	
	// Compute the directions of movement
	VECTOR3 dir1(   corners[1].x() - corners[0].x(),  corners[1].y() - corners[0].y(),  corners[1].z() - corners[0].z() );
	VECTOR3 dir2(   corners[3].x() - corners[0].x(),  corners[3].y() - corners[0].y(),  corners[3].z() - corners[0].z() );
	VECTOR3 dir3(   corners[4].x() - corners[0].x(),  corners[4].y() - corners[0].y(),  corners[4].z() - corners[0].z() );
	
	// Compute the total number of seeds to be generated on the plane
	this->nSeed = 125; //(int) ( floor( spanX ) * floor( spanY ) * floor( spanZ ) ); 
	int nseed1 = 5; //(int) floor( spanX );
	int nseed2 = 5; //(int) floor( spanY );
	int nseed3 = 5; //(int) floor( spanZ );
		
	//cout << "dir1: " << dir1[0] << " " << dir1[1] << " " << dir1[2] << endl;
	//cout << "dir2: " << dir2[0] << " " << dir2[1] << " " << dir2[2] <<  endl;
	//cout << "nseeds: " << nseed1 << " " << nseed2 << endl;
		
	#ifdef DEBUG_MODE
	printf( "Total number of seeds to be generated: %d\n", this->nSeed );
	#endif

	// Allocate memory for seeds
	if( this->seeds != NULL ) delete [] this->seeds;
	this->seeds = new VECTOR3[this->nSeed];

	// Increment or step size along the plane
	float delt1 = 1.0f / (float)( nseed1-1 );
	float delt2 = 1.0f / (float)( nseed2-1 );
	float delt3 = 1.0f / (float)( nseed3-1 );
	#ifdef DEBUG_MODE
	//printf( "Step size: %f %f\n", delt1, delt2 );
	#endif

	int seedIndex = 0;
	float t1 = 0;
	float t2 = 0;
	float t3 = 0;
	VECTOR3 base;
	VECTOR3 base2;
	for( int i=0; i<nseed1; i++ )
	{
		// Move one step along first direction
		base.Set( corners[0].x() + t1 * dir1[0], corners[0].y() + t1 * dir1[1], corners[0].z() + t1 * dir1[2] );
		
		t2 = 0;
		for( int j=0; j<nseed2; j++ )
		{
			// Move one step along second direction
			base2.Set(  base[0] + t2 * dir2[0], base[1] + t2 * dir2[1], base[2] + t2 * dir2[2] );
			
			t3 = 0;
			for( int j=0; j<nseed2; j++ )
			{
				// Move one step along third direction
				this->seeds[seedIndex].Set(  base2[0] + t3 * dir3[0], base2[1] + t3 * dir3[1], base2[2] + t3 * dir3[2] );
		
				// Increment t3
				t3 += delt3;
				
				// Move to next seed
				seedIndex++;			
			
			}
			
			// Increment t2
			t2 += delt2;
		}
		
		// Increment t1
		t1 += delt1;
	}

	// Print the seed points (optional)
	#ifdef DEBUG_MODE
	printf( "Seeds generated along line\n" );
	printf( "Printing seeds ...\n" );
	for (int i=0; i<this->nSeed; i++)
		printf("Seed %d: [%f %f %f]\n", i, this->seeds[i][0], this->seeds[i][1], this->seeds[i][2]);
	#endif
	
}

void
FDL_osuflow::generateSeedInBox( float* low, float* high, int nseed )
{
	printf( "Generating random seeds within a subvolume...\n");
	this->nSeed = nseed;

	// Allocate memory for seeds
	if( this->seeds != NULL ) delete [] this->seeds;
	this->seeds = new VECTOR3[this->nSeed];

	// Generate seeds
	srand ( time(NULL) );
	for( int i=0; i<this->nSeed; i++ )
	{
		seeds[i].Set( low[0] + (high[0]-low[0])*( rand() % 1000 ) / (float)1000.0f,
					  low[1] + (high[1]-low[1])*( rand() % 1000 ) / (float)1000.0f,
					  low[2] + (high[2]-low[2])*( rand() % 1000 ) / (float)1000.0f );
	}

	// Print the seed points (optional)
	#ifdef DEBUG_MODE
	printf( "Seeds generated along line\n" );
	printf( "Printing seeds ...\n" );
	for (int i=0; i<this->nSeed; i++ )
		printf("Seed %d: [%f %f %f]\n", i, this->seeds[i][0], this->seeds[i][1], this->seeds[i][2]);
	#endif

}// end function

void
FDL_osuflow::generateSeedsRandom( int nSeed )
{
	printf( "Generating random seeds...\n");
	this->nSeed = nSeed;

	float *low = new float[3];
	float *high = new float[3];
	low[0] = dataLow[0]; low[1] = dataLow[1]; low[2] = dataLow[2];
	high[0] = dataHigh[0]; high[1] = dataHigh[1]; high[2] = dataHigh[2];

	osuflow->SetRandomSeedPoints( low, high, this->nSeed );

	if( this->seeds != NULL )	delete [] this->seeds;
	this->seeds = osuflow->GetSeeds( this->nSeed );

	delete [] low;
	delete [] high;
}// end function

void FDL_osuflow::readSeeds( int nSeed, const char *seedFileName )
{
	// Open file
	this->seedFile = fopen( seedFileName, "r" );
	assert( seedFile != NULL );

	// Allocate memory for seeds
	this->nSeed = nSeed;
	this->seeds = new VECTOR3[this->nSeed];
	#ifdef DEBUG_MODE
	printf( "Memory allocated for %d seeds\n", this->nSeed );
	#endif

	// Read seeds (one per line)
	float x, y, z;
	for( int i=0; i<this->nSeed; i++ )
	{
		fscanf( seedFile, "%f, %f, %f", &x, &y, &z );
		//printf( "%f, %f, %f\n", x, y, z );
		this->seeds[i].Set( x, y, z );
	}// end for

	// Close file
	fclose( this->seedFile );

}// end function

void FDL_osuflow::writeSeeds( const char *seedFileName )
{
	assert( this->nSeed > 0 );
	assert( this->seeds != NULL );

	// Open file for writing seeds
	this->seedFile = fopen( seedFileName, "w" );

	// Write total number of seeds in the first line
	// Optional:: may create problem while reading as a csv file
	// fprintf( this->seedFile, "%d\n", this->nSeed );

	// Write seeds
	for( int i=0; i<this->nSeed; i++ )
		fprintf( this->seedFile, "%f, %f, %f\n", this->seeds[i].x(),
											     this->seeds[i].y(),
											     this->seeds[i].z() );

	// Close output file
	fclose( this->seedFile );

}// end function

void FDL_osuflow::advect( int nStep, float p1, float p2 )
{
	VECTOR3* fwdTraceData = NULL;
	VECTOR3* backTraceData = NULL;
	VECTOR3* totalTrace = NULL;
	list<vtListSeedTrace*> sl_list;
	int iP = 0;
	int iT = 0;
	int fwdLen, backLen;

	assert( this->nSeed > 0 );
	assert( this->seeds != NULL );

	// Allocate memory for storing length of each streamline
	this->streamlineLengthList = new int[this->nSeed];

	#ifdef DEBUG_MODE
	printf( "Computing streamlines..\n" );
	#endif

	// Advect
	this->osuflow->SetIntegrationParams( p1, p2 );
	this->osuflow->GenStreamLines( this->seeds, BACKWARD_AND_FORWARD, this->nSeed, nStep, sl_list );

	#ifdef DEBUG_MODE
	printf( "Integration done.\n" );
	#endif

	// Create list of streamlines from advection output
	for( list<vtListSeedTrace*>::const_iterator pIter = sl_list.begin();
		 pIter!=sl_list.end();
		 pIter++, iT++ )
	{
		// Get forward trace data and copy it to the streamline
		vtListSeedTrace *fwdTrace = *pIter;
		if( fwdTrace->size() == 0 )
		{
			//printf( "f: %f %f %f\n", this->seeds[iT][0], this->seeds[iT][1], this->seeds[iT][2] );
			fwdTraceData = new VECTOR3[1];
			fwdTraceData->Set( this->seeds[iT][0], this->seeds[iT][1], this->seeds[iT][2] );
			fwdLen = 1;
		}
		else
		{
			fwdTraceData = new VECTOR3[(int)(fwdTrace->size()-1)];
			iP = 0; fwdLen = 0;
			for( list<VECTOR3*>::const_iterator pnIter = fwdTrace->begin();
				 iP < (fwdTrace->size()-1);
			     pnIter++, iP++ )
			{
				fwdTraceData[iP] = **pnIter;
				fwdLen++;
			}
		}// end else: forward trace

		// Advance pointer
		pIter++;

		// Get backward trace data and copy it to the streamline
		vtListSeedTrace *backTrace = *pIter;
		if( backTrace->size() == 0 )
		{
			//printf( "b: %f %f %f\n", this->seeds[iT][0], this->seeds[iT][1], this->seeds[iT][2] );
			backTraceData = new VECTOR3[1];
			backTraceData->Set( this->seeds[iT][0], this->seeds[iT][1], this->seeds[iT][2] );
			backLen = 1;
		}
		else
		{
			backTraceData = new VECTOR3[(int)(backTrace->size()-1)];
			iP = 0; backLen = 0;
			for( list<VECTOR3*>::const_iterator pnIter = backTrace->begin();
				iP < (backTrace->size()-1);
				pnIter++, iP++ )
			{
				backTraceData[iP] = **pnIter;
				backLen++;
			}
		}// end else: back trace

		// Join the two traces
		totalTrace = this->joinTwoTraceDirections( fwdTraceData, fwdLen, backTraceData, backLen );

		// Store the length of the current trace
		this->streamlineLengthList[iT] = fwdLen + backLen - 1;
		//printf( "%d\n", this->streamlineLengthList[iT] );

		// Add the trace to the list of all traces
		this->streamlineList.push_back( totalTrace );

		// Free
		// delete [] fwdTraceData;
		// delete [] backTraceData;

	}// end for

	// OSUFlow clean up operation **** VVI ****
	for( list<vtListSeedTrace*>::const_iterator pIter = sl_list.begin();
		 pIter!=sl_list.end();
		 pIter++ )
	{
		vtListSeedTrace *temp = *pIter;
		// Clear each vector in the nested list
		for( list<VECTOR3*>::const_iterator pnIter = temp->begin();
			 pnIter != temp->end();
			 pnIter++ )
			delete *pnIter;

		// Clear the nested list
		(*temp).clear();
	}
	sl_list.clear();

}// end function

void FDL_osuflow::advect_unidirection( int nStep, float p1, float p2 )
{
	VECTOR3* fwdTraceData = NULL;
	list<vtListSeedTrace*> sl_list;
	int iP = 0;
	int iT = 0;
	int fwdLen;

	assert( this->nSeed > 0 );
	assert( this->seeds != NULL );

	// Clear up memory if already allocated
	if( this->streamlineLengthList != NULL )	delete [] this->streamlineLengthList;
	while( !this->streamlineList.empty() )
	{
		delete this->streamlineList.back();
		this->streamlineList.pop_back();
	}

	// Allocate memory for storing length of each streamline
	this->streamlineLengthList = new int[this->nSeed];

	#ifdef DEBUG_MODE
	printf( "Computing streamlines..\n" );
	#endif

	// Advect
	this->osuflow->SetIntegrationParams( p1, p2 );
	this->osuflow->GenStreamLines( this->seeds, FORWARD_DIR, this->nSeed, nStep, sl_list );

	#ifdef DEBUG_MODE
	printf( "Integration done.\n" );
	#endif

	// Create list of streamlines from advection output
	for( list<vtListSeedTrace*>::const_iterator pIter = sl_list.begin();
		 pIter!=sl_list.end();
		 pIter++, iT++ )
	{
	    // Get forward trace data and copy it to the streamline
		vtListSeedTrace *fwdTrace = *pIter;
		if( fwdTrace->size() == 0 )
		{
			//printf( "f: %f %f %f\n", this->seeds[iT][0], this->seeds[iT][1], this->seeds[iT][2] );
			fwdTraceData = new VECTOR3[1];
			fwdTraceData->Set( this->seeds[iT][0], this->seeds[iT][1], this->seeds[iT][2] );
			fwdLen = 1;
		}
		else
		{
			fwdTraceData = new VECTOR3[(int)(fwdTrace->size()-1)];
			iP = 0; fwdLen = 0;
			for( list<VECTOR3*>::const_iterator pnIter = fwdTrace->begin();
				 iP < (fwdTrace->size()-1);
			     pnIter++, iP++ )
			{
				fwdTraceData[iP] = **pnIter;
				//printf( "advecting: %f %f %f %f\n", fwdTraceData[iP][0], fwdTraceData[iP][1],
				//						 fwdTraceData[iP][2], fwdTraceData[iP][3] );
				fwdLen++;
			}
		}// end else: forward trace

	    // Store the length of the current trace
	    this->streamlineLengthList[iT] = fwdLen;
	    //printf( "%d\n", this->streamlineLengthList[iT] );

	    // Add the trace to the list of all traces
	    this->streamlineList.push_back( fwdTraceData );

	}// end for

	// OSUFlow clean up operation **** VVI ****
	for( list<vtListSeedTrace*>::const_iterator pIter = sl_list.begin();
		 pIter!=sl_list.end();
		 pIter++ )
	{
		vtListSeedTrace *temp = *pIter;
		// Clear each vector in the nested list
		for( list<VECTOR3*>::const_iterator pnIter = temp->begin();
			 pnIter != temp->end();
			 pnIter++ )
			delete *pnIter;

		// Clear the nested list
		(*temp).clear();
	}
	sl_list.clear();

}// end function

VECTOR3* FDL_osuflow::joinTwoTraceDirections( VECTOR3 *fwdTraceData, int fwdTraceLength,
											  VECTOR3 *backTraceData, int backTraceLength )
{
	// Compute total length of a trace
	int totalTraceLength = fwdTraceLength + backTraceLength - 1;

	// Allocate memory for total trace
	VECTOR3* totalTrace = new VECTOR3[totalTraceLength];

	#ifdef DEBUG_MODE
	printf( "Next streamline\n" );
	printf( "Printing back trace\n" );
	for( int i=0; i<backTraceLength; i++ )
		printf( "%f %f %f\n", backTraceData[i][0], backTraceData[i][1], backTraceData[i][2] );
	printf( "Printing forward trace\n" );
	for( int i=0; i<fwdTraceLength; i++ )
		printf( "%f %f %f\n", fwdTraceData[i][0], fwdTraceData[i][1], fwdTraceData[i][2] );
	#endif

	// Copy location information from streamline
	// Goes from -1 to -bakT and then 0 to fwdT
	for( int i=0; i<backTraceLength; i++ )
		totalTrace[i] = backTraceData[backTraceLength-1-i];
	for( int i=0; i<fwdTraceLength-1; i++ )
		totalTrace[i+backTraceLength] = fwdTraceData[i+1];

	return totalTrace;

}// end function

void FDL_osuflow::writeStreamlines( const char *streamlineFileName )
{
	// Open file for saving entropy field
	this->streamlineFile = fopen( streamlineFileName, "wb" );
	assert( this->streamlineFile != NULL );

	// Write number of streamlines
	fwrite( &this->nSeed, sizeof(int), 1, this->streamlineFile );

	// Write length of each streamline
	fwrite( this->streamlineLengthList, sizeof(int), this->nSeed, this->streamlineFile );

	// Write the actual traces
	int i=0;
	for( list<VECTOR3*>::const_iterator pIter = this->streamlineList.begin();
		 pIter != this->streamlineList.end();
		 pIter++, i++ )
	{
		fwrite( *pIter, sizeof(VECTOR3), this->streamlineLengthList[i], this->streamlineFile );
		VECTOR3* tmp = *pIter;
		//for( int j=0; j<this->streamlineLengthList[i]; j++ )
		//	printf( "writing: %f %f %f %f\n", tmp[j][0], tmp[j][1],
		//								 tmp[j][2], tmp[j][3] );
	}


	// Close file
	fclose( this->streamlineFile );

}// end function

FDL_osuflow::~FDL_osuflow()
{
	if( this->osuflow != NULL ) 				delete this->osuflow;
	if( this->seeds != NULL )					delete [] this->seeds;
	if( this->streamlineLengthList != NULL )	delete [] this->streamlineLengthList;
	for( list<VECTOR3*>::const_iterator pIter = this->streamlineList.begin();
		 pIter != this->streamlineList.end();
		 pIter++ )
		delete [] *pIter;
	this->streamlineList.clear();

}// end desctructor

