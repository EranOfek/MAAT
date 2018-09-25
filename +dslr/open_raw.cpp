#include "libraw.h"
#include "mex.h"
#include "matrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] ){

	
	if(nrhs<1){ printf("USAGE: I = open_raw(filename)\n"); return; }
  
	if(mxIsChar(prhs[0])==0) mexErrMsgIdAndTxt("MATLAB:nikon:scan_nef:filename_not_string","Must input a filename string!");
  
	char *filename=mxArrayToString(prhs[0]);

	// Let us create an image processor
	LibRaw iProcessor;
  //Linearization Step
    iProcessor.imgdata.params.output_bps = 16;
    
  //SET BRIGHTNESS VALUE (default = 1.0, divides white value)   
    iProcessor.imgdata.params.bright = 1.0;
    
  //SET WHITE BALANCE MULTIPLIERS
  	iProcessor.imgdata.params.user_mul[0] = 300;
    iProcessor.imgdata.params.user_mul[1] = 300;
    iProcessor.imgdata.params.user_mul[2] = 300;
    iProcessor.imgdata.params.user_mul[3] = 300;
    
    //set gamma (set both to 1 for linear)
    iProcessor.imgdata.params.gamm[0] = 1;
    iProcessor.imgdata.params.gamm[1] = 1;
       
  //set saturation value (integers < 1, default = -1)
    iProcessor.imgdata.params.user_sat = -1;    
        
    
	// Open the file and read the metadata
	iProcessor.open_file(filename);

      
	// Let us unpack the image
	iProcessor.unpack();

    //process unpacked image with parameters set above
    iProcessor.dcraw_process(); 
    
    
    
          /* printf("WB Coeff 0: %f\n",iProcessor.imgdata.color.cam_mul[0]);
            printf("WB Coeff 1: %f\n",iProcessor.imgdata.color.cam_mul[1]);
            printf("WB Coeff 2: %f\n",iProcessor.imgdata.color.cam_mul[2]);
            printf("WB Coeff 3: %f\n",iProcessor.imgdata.color.cam_mul[3]);
            printf("Shutter Speed: %f\n",iProcessor.imgdata.other.shutter);
            printf("Aperture: %f\n",iProcessor.imgdata.other.aperture);
            printf("ISO: %f\n",iProcessor.imgdata.other.iso_speed);
            printf("sat: %i\n",iProcessor.imgdata.params.user_sat); */
            
           
            
	// Finally, let us free the image processor for work with the next image

	size_t dims[3] = {0};
	
	dims[0]=4;
	dims[1]=iProcessor.imgdata.sizes.iwidth;
	dims[2]=iProcessor.imgdata.sizes.iheight;
	
	plhs[0]=mxCreateNumericArray(3, (const size_t*)dims, mxUINT16_CLASS, mxREAL);
  
	// printf("dims[0]= %d | dims[1]= %d | dims[2]= %d\n", dims[0], dims[1], dims[2]);
	memcpy(mxGetData(plhs[0]), iProcessor.imgdata.image, dims[0]*dims[1]*dims[2]*sizeof(ushort));

  	iProcessor.recycle();

	// ushort *output=(ushort*) mxGetData(plhs[0]);
	// for(int i=0; i<dims[0]*dims[1]*dims[2]; i++) {
	  
		// printf("i= %d\n", i);
		// output[i]=data[i];
	// }
  
  
}
