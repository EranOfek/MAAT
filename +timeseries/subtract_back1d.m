function VecSub=subtract_back1d(Vec,Algo,varargin)
%------------------------------------------------------------------------------
% subtract_back1d function                                             General
% Description: Subtract background level from a 1-D vector.
% Input  : - Vector from which to subtract background.
%            If two columns matrix then the first column indicate position,
%            while the second column is for the value from which to subtract
%            the background.
%          - Algorithm:
%            'none'     - do not subtract background.
%            'mean'     - subtract the mean.
%            'median'   - subtract the median.
%            'medfilt'  - subtract a median filter, where the block size
%                         is provided by the next input argument.
%          * Additional parameters to pass to the subtraction algorithm.
% Output : - A background subtracted vector.
%            If the input is a two column matrix than the output will
%            be a two column matrix.
% Tested : Matlab R2012A
%     By : Eran O. Ofek                    Oct 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: VecSub=subtract_back1d(rand(100,1),'median');
% Reliable: 2
%------------------------------------------------------------------------------

[Ni,Nj] = size(Vec);
if (Nj==1),
   VecPix = [1:1:Ni].';
   VecVal = Vec;
elseif (Nj>1),
   VecPix = Vec(:,1);
   VecVal = Vec(:,2);
else
   error('Spec input should be a column vector or two columns matrix');
end

switch lower(Algo)
 case 'none'
    VecSub = VecVal;
 case 'median'
    VecSub  = VecVal(:,1) - median(VecVal(:,1));
 case 'mean'
    VecSub  = VecVal(:,1) - mean(VecVal(:,1));
 case 'medfilt'
    VecSub  = VecVal(:,1) - medfilt1(VecVal(:,1),varargin{:});
 otherwise
    error('Unknown Back option');
end

if (Nj>1),
   VecSub = [VecPix, VecSub];
end
