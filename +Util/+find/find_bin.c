#include "mex.h"
#include <math.h>

/*
 * find_bin.c
 * Description: given a sorted vector and a vector, find the
 *              the index of the element in the sorted vector which
 *              is nearest to each one of the elements in the vector.
 *
 * Example: I=find_bin(R,0.5);
 *
 * compile using: mex find_bin.c
 *
 * Eran O. Ofek           March 2011
 *
 *
 */


void find_bin_c(double *X, double *Val, double *IndVal, mwSize N, mwSize Nv)
{
    
  int Low, High, Middle;
  int Counter;
  int Found;
  int iV;
  
  for (iV=0; iV<Nv; iV++) {
     High      = N-1;
     Low       = 1-1;
   
     Found  = 0;
     Counter = 0;
     while (Low<=High & Found==0) {
         
          Counter = Counter + 1;
          Middle = floor(Low + (High-Low)*0.5);
          
          if (Val[iV]<X[Middle]) {
              High = Middle - 1;
          } else if (Val[iV]>X[Middle]) {
              Low = Middle + 1;
          } else {
              Found  = 1;
          }
     }
   IndVal[iV] = Middle+1;
 }          



  
  /*  
  
  int Ind1, Ind2, IndM;
  double Y1, Y2, Ym;
  int Found;
    
   if (N==1) {
      *IndVal = 1;
   }
   else {
      Ind1   = 1;
      Ind2   = N;
      IndM   = floor(0.5*N);
      Y1     = X[Ind1-1];
      Y2     = X[Ind2-1];
      Ym     = X[IndM-1];
   
      Found  = 0;
      while (Found==0) {
         if (Val>Ym) {
            Ind1 = IndM;
            Y1   = X[Ind1-1];
         
            if ((Ind2-Ind1)>=2) {
               IndM = floor(floor(0.5*(Ind2+Ind1)));
            } else {
               Found = 1;
               if (abs(Val-Y1)<abs(Val-Y2)) {
                  *IndVal = Ind1;
               } else {
                  *IndVal = Ind2;
               }
            }
         
            Ym   = X[IndM-1];
         } else if (Val<Ym) {
            Ind2 = IndM;
            Y2   = X[Ind2-1];

            if ((Ind2-Ind1)>=2) {
               IndM = floor(floor(0.5*(Ind2+Ind1)));
            } else {
               Found = 1;
               if (abs(Val-Y1)<abs(Val-Y2)) {
                  *IndVal = Ind1;
               } else {
                  *IndVal = Ind2;
               }
            }
         
            Ym   = X[IndM-1];

         } else {
            Found  = 1;
            *IndVal = IndM;
         }
      }
   }
  */
   
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  double *Vec;
  double *Ind;
  double *Val;
  mwSize N;
  mwSize Nv;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=2) 
    mexErrMsgTxt("Two inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");
  
  /* check to make sure the second input argument is a scalar */
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) {
    /*      mxGetN(prhs[1])*mxGetM(prhs[1])!=1 ) { */
    mexErrMsgTxt("Input x must be a scalar.");
  }
  
  /*  get the scalar input x */
  Val = mxGetPr(prhs[1]);
  /*  get the dimensions of the matrix input y */
  Nv = mxGetM(prhs[1]);
  
  
  /*  create a pointer to the input vector Vec */
  Vec = mxGetPr(prhs[0]);
  
  /*  get the dimensions of the matrix input y */
  N = mxGetM(prhs[0]);

  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(Nv,1, mxREAL);
  
  /*  create a C pointer to a copy of the output matrix */
  Ind = mxGetPr(plhs[0]);
  
  /*  call the C subroutine */
  find_bin_c(Vec,Val,Ind,N,Nv);
  
}

