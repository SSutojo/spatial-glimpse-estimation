#include <math.h>
#include "mex.h"

/* Helper functions */
#define max(x, y)   ((x) > (y) ? (x) : (y))
#define	min(A, B)	((A) < (B) ? (A) : (B))
#define swap(A,B)   temp = (A); (A)=(B); (B) = temp;


/* Input Arguments */
#define   INPUT                prhs[0]  // Input vector or matrix
#define   MODE                 prhs[1]  // Mode

/* Output Arguments */
#define   LOCALIDX        	   plhs[0]  // all local maximum 
#define   MAXIDX        	   plhs[1]  // global maximum per block
#define   LOCALITERIDX         plhs[2]  // 
#define   MAXITERDX        	   plhs[3]  // 

/*  Set resolution of floating point numbers */
const double EPS = pow(2.00,-52.00);

void usage()
{
	mexPrintf("\n findLocalPeaks   MEX-File for local and global maximum search.");		
	mexPrintf("\n");
	mexPrintf("\n\t A value is considered as local/global peak if the adjacent amplitudes ");
	mexPrintf("\n\t decrease.");
	mexPrintf("\n");
	mexPrintf("\n USAGE");
	mexPrintf("\n 	[LOCMAXIDX,MAXIDX,LOCMAXITER,MAXITER] = findLocalPeaks(IN,MODE);");
	mexPrintf("\n");
	mexPrintf("\n INPUT ARGUMENTS");
	mexPrintf("\n\t    IN : input matrix [nSamples x nChannels]");
	mexPrintf("\n\t  MODE : mode for maxima search");
	mexPrintf("\n\t         1 = endpoints are not accepted as local or global maxima. So the");
	mexPrintf("\n\t         number of global maxima does not necessarily corresponds to the ");
	mexPrintf("\n\t         number of channels as there might be channels without any ");
	mexPrintf("\n\t         maximum.");
	mexPrintf("\n");
	mexPrintf("\n\t         2 = also endpoints are accepted as local or global maxima. The"); 
	mexPrintf("\n\t         number of global maxima always corresponds to the number of ");
	mexPrintf("\n\t         channels. (default, MODE = 2)");
	mexPrintf("\n");
	mexPrintf("\n OUTPUT ARGUMENTS");
	mexPrintf("\n\t  LOCMAXIDX : indices of all local maxima. If IN is a matrix the");
	mexPrintf("\n\t              can be used to directly access the maxima.");
	mexPrintf("\n\t     MAXIDX : index of global maximum per channel");
	mexPrintf("\n\t LOCMAXITER : indices identifying the channel of local maxima.");
	mexPrintf("\n\t    MAXITER : indices identifying the channel of global maxima.");
	mexPrintf("\n");
	mexPrintf("\n EXAMPLE");
	mexPrintf("\n\t %s Create input data","%");
	mexPrintf("\n\t blockSize = 100;");             
	mexPrintf("\n\t in        = rand(blockSize,10);    ");
	mexPrintf("\n\t in        = filter([0.15 0.15],[1 -0.6],in);");
	mexPrintf("\n");
	mexPrintf("\n\t %s Find local and global maxima","%");
	mexPrintf("\n\t [a,b,c,d] = findLocalPeaks(in);");
	mexPrintf("\n\t %s Use local and global indices for plotting on a 2D surface","%");
	mexPrintf("\n\t a = mod(a-1,blockSize)+1;");
	mexPrintf("\n\t b = mod(b-1,blockSize)+1;");
	mexPrintf("\n");
	mexPrintf("\n\t %s Plot 2D surface","%");
	mexPrintf("\n\t figure;hold on;");
	mexPrintf("\n\t imagesc(in);axis xy");
	mexPrintf("\n\t plot(c,a,'xw','MarkerSize',4)");
	mexPrintf("\n\t plot(d,b,'*w','MarkerSize',8)");
	mexPrintf("\n\t axis tight; colorbar; ");
	mexPrintf("\n");
	mexPrintf("\n\tImplementation by Tobias May ? 2008");
	mexPrintf("\n\tPlease send bug reports to: tobias.may@philips.com");
	mexPrintf("\n"); 
}    


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *input, *maxIdx, *maxIterIdx, *localIdx, *localIterIdx;
	double currMax, currMaxIdx;
	int    ii, jj, nSamples, nChannels, mode;
	int    localMaxCtr, globalMaxCtr, maxCtr, nLocPeaks;
 	int    exit = false;

	// Check for proper number of arguments
  	if ((nrhs < 1) || (nrhs > 2) || (nlhs > 4)){
		usage();
    	exit = true;   
  	} 

 	if(!exit){
  		// Check the dimensions of the feature space
	  	nSamples     = (int)mxGetM(INPUT);
  		nChannels    = (int)mxGetN(INPUT);
 
		if (nSamples <= 1)
			mexErrMsgTxt("ERROR@findLocalPeaks.dll: IN must contain at least two samples per column.");
		else			
			input = mxGetPr(INPUT);

		if (nrhs<2)
			mode = 2;
		else
			mode = (int) mxGetScalar(MODE);

		if ((mode != 1) && (mode != 2)){
			usage();
			mexErrMsgTxt("ERROR@findLocalPeaks.dll: MODE must be either one or two.");
		}


		// Create a matrix for the return argument 
		MAXIDX       = mxCreateDoubleMatrix(nChannels, 1, mxREAL);
		MAXITERDX    = mxCreateDoubleMatrix(nChannels, 1, mxREAL);
		LOCALIDX     = mxCreateDoubleMatrix(nSamples * nChannels, 1, mxREAL);
		LOCALITERIDX = mxCreateDoubleMatrix(nSamples * nChannels, 1, mxREAL);

		// Asign pointers
		// NOTE: The size of all these output pointers needs to be re-allocated
		//       as the number of peaks is unknown (using -> mxRealloc)
  		maxIdx       = mxGetPr(MAXIDX);
		maxIterIdx   = mxGetPr(MAXITERDX);
		localIdx     = mxGetPr(LOCALIDX);
		localIterIdx = mxGetPr(LOCALITERIDX);

		// Initialize counters
		localMaxCtr  = 0;
		globalMaxCtr = 0;
		maxCtr       = 0;

		if (mode == 1){
			
		// Loop over the number of channels
		for (ii = 0; ii < nChannels; ii++){
			
			// (Re-)Initialize global maximum parameters 
			currMax    = -99999999;
			currMaxIdx = -99999999;

			nLocPeaks  = 0;

			// Loop over number of samples		
			for (jj = 1; jj < nSamples-1; jj++){

				// Check if new local maximum is found ...
				if (input[jj + ii*nSamples] > input[jj + 1 + ii*nSamples] && (input[jj + ii*nSamples] > input[jj - 1 + ii*nSamples]))
				{
					// Store index of local maximum
					localIdx[localMaxCtr] = jj + 1 + ii * nSamples;
						
					// Store corresponding row index (useful for 2D plotting)
					localIterIdx[localMaxCtr] = ii + 1;

					// Increase local maximum counter (global)
					localMaxCtr = localMaxCtr + 1;

					// Increase channel-based local maximum counter
					nLocPeaks = nLocPeaks + 1;

					// Check for maximum peak of current column
					if (input[jj + ii*nSamples] > currMax){

						// Store current maximum
						currMax = input[jj + ii*nSamples];
							
						// Store maximum location 
						// Small reminder : C/C++ -> start index 0, that's why real index in Matlab is jj+1
						currMaxIdx = jj + 1;
					}
				}
			}
				
			// Note that this implementation will always find a maximum peak per column.
			// In contrast to simply comparing the amplitude height, the local peak has the 
			// additional contraint that the adjacent amplitude values must decrease ....
	
			/*
			// Check if first amplitude value is larger then the current maximum
			if (input[0 + ii*nSamples] > currMax){
  				// Store maximum location 
				currMaxIdx = 1;
			}
			*/

			if (currMax > -99999999)
			{
				// Store overall maximum for current column
				maxIdx[globalMaxCtr] = currMaxIdx + ii * nSamples;
				// Small reminder : C/C++ -> start index 0, that's why real index in Matlab is ii+1
				maxIterIdx[globalMaxCtr] = ii + 1;

			    // Increase global maximum counter
				globalMaxCtr++;
			}

			// If no local peak was found due to the constraints ... insert the channel maximum
			// so that the number of local peaks is always >= the number of the channel-dependent
			// maximum

			/*
			if (nLocPeaks==0){
				// Store index of local maximum
				localIdx[localMaxCtr] = maxIdx[ii];
					
				// Store corresponding row index (useful for 2D plotting)
				localIterIdx[localMaxCtr] = ii + 1;

				// Increase local maximum counter (global)
				localMaxCtr++;
			}
			*/
		}

		// Change array size to fit the number of found maxima
		mxSetM(LOCALIDX,     localMaxCtr);
		mxSetM(LOCALITERIDX, localMaxCtr);

		mxSetM(MAXIDX,       globalMaxCtr);
		mxSetM(MAXITERDX,    globalMaxCtr);
	    } 
		else{

			// Loop over the number of channels
			for (ii = 0; ii < nChannels; ii++){

				// (Re-)Initialize global maximum parameters 
				currMax    = -99999999;
			    currMaxIdx = -99999999;

				nLocPeaks  = 0;

				// Note that this implementation will always find a maximum peak per column.

				// ****************************
				// Include endpoints ...
				// ****************************

				// Check starting point
				if (input[0 + ii*nSamples] > input[1 + ii*nSamples]){
					// Store index of local maximum
					localIdx[localMaxCtr] = 1 + ii * nSamples;
					
					// Store corresponding row index (useful for 2D plotting)
  				    localIterIdx[localMaxCtr] = ii + 1;

					// Increase local maximum counter (global)
                    localMaxCtr = localMaxCtr + 1;
						
					// Increase channel-based local maximum counter
					nLocPeaks = nLocPeaks + 1;

					// Store current maximum
					currMax = input[0 + ii*nSamples];
							
					// Store maximum location 
					currMaxIdx = 1;
				}					

				// Loop over number of samples		
				for (jj = 1; jj < nSamples-1; jj++){
					
					// Check if new local maximum is found ...
					if (input[jj + ii*nSamples] > input[jj + 1 + ii*nSamples] && (input[jj + ii*nSamples] > input[jj - 1 + ii*nSamples]))
					{
						// Store index of local maximum
						localIdx[localMaxCtr] = jj + 1 + ii * nSamples;
						
						// Store corresponding row index (useful for 2D plotting)
						localIterIdx[localMaxCtr] = ii + 1;

						// Increase local maximum counter (global)
                        localMaxCtr = localMaxCtr + 1;
						
						// Increase channel-based local maximum counter
						nLocPeaks = nLocPeaks + 1;

						// Check for maximum peak of current column
						if (input[jj + ii*nSamples] > currMax){
							
							// Store current maximum
							currMax = input[jj + ii*nSamples];
							
							// Store maximum location 
							// Small reminder : C/C++ -> start index 0, that's why real index in Matlab is jj+1
							currMaxIdx = jj + 1;
						}
					}
				}

				// ****************************
				// Include endpoints ...
				// ****************************
				
				// Include endpoint
				if (input[nSamples - 1 + ii*nSamples] > input[nSamples - 2 + ii*nSamples]){
					// Store index of local maximum
					localIdx[localMaxCtr] = nSamples + ii * nSamples;
					
					// Store corresponding row index (useful for 2D plotting)
  				    localIterIdx[localMaxCtr] = ii + 1;

					// Increase local maximum counter (global)
                    localMaxCtr = localMaxCtr + 1;
						
					// Increase channel-based local maximum counter
					nLocPeaks = nLocPeaks + 1;

					// Check for maximum peak of current column
					if (input[nSamples - 1 + ii*nSamples] > currMax){
						// Store current maximum
						currMax = input[nSamples - 1 + ii*nSamples];
						
						// Store maximum location 
						currMaxIdx = nSamples;
					}
				}				

				if (nLocPeaks==0){
					localIdx[localMaxCtr]     = 1 + ii*nSamples;
					localIterIdx[localMaxCtr] = 1;
					localMaxCtr               = localMaxCtr + 1;
                    nLocPeaks                 = nLocPeaks + 1;
					currMax                   = input[0 + ii*nSamples];
					currMaxIdx                = 1;
				}


				/*
				// Check if first amplitude value is larger then the current maximum
				if (input[0 + ii*nSamples] > currMax){

					currMax = input[0 + ii*nSamples];
					// Store maximum location 
					currMaxIdx = 1;
				}
				
				// Check if last amplitude value is larger then the current maximum
				if (input[jj + ii*nSamples] > currMax){
					// Store maximum location 
					currMaxIdx = jj + 1;
				}
				*/

				// Store overall maximum for current column
				maxIdx[ii] = currMaxIdx + ii * nSamples;
				// Small reminder : C/C++ -> start index 0, that's why real index in Matlab is ii+1
				maxIterIdx[ii] = ii + 1;

				// If no local peak was found due to the constraints ... insert the channel maximum
				// so that the number of local peaks is always >= the number of the channel-dependent
				// maximum
				
				if (nLocPeaks==0){
					// Store index of local maximum
					localIdx[localMaxCtr] = maxIdx[ii];
			
					// Store corresponding row index (useful for 2D plotting)
					localIterIdx[localMaxCtr] = ii + 1;
					
					// Increase local maximum counter (global)
					localMaxCtr++;
				}
			}
			
			// Change array size to fit the number of found maxima
			mxSetM(LOCALIDX,     localMaxCtr);
			mxSetM(LOCALITERIDX, localMaxCtr);
		}
		}
					
} // end of mex file