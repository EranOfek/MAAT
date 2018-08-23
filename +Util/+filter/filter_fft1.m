function [CS]=filter_fft1(S,Freq,CutLowFreq,Smooth,Dim)
% Filter a 1D equally spaced series using FFT.
% Package: Util.cell
% Description: Filter a 1D equally spaced series using FFT.
% Input  : - Series.
%          - Cutoff frequency in fraction of total length.
%            Default is 0.1 of series length.
%          - Low or high frequencies {true|false}. true for cutting low
%            frequencies. false for removing high frequencies.
%            Default is true.
%          - Tapering smoothness (using a Fermi function) in units
%            of series length. Default is 0.0001.
%          - If series contains more than one dimension, than this is
%            the dimension at which to operate. Default is 1.
% Output : - Cleaned series.
% Tested : Matlab R2013a
%     By : Eran O. Ofek                    Aug 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [CS]=Util.cell.filter_fft1(A);
% Reliable: 2
%--------------------------------------------------------------------------


Def.Freq       = 0.1;
Def.CutLowFreq = true;
Def.Smooth     = 0.0001;
Def.Dim        = 1;
if (nargin==1)
    Freq       = Def.Freq;
    CutLowFreq = Def.CutLowFreq;
    Smooth     = Def.Smooth;
    Dim        = Def.Dim;
elseif (nargin==2)
    CutLowFreq = Def.CutLowFreq;
    Smooth     = Def.Smooth;
    Dim        = Def.Dim;
elseif (nargin==3)
    Smooth     = Def.Smooth;
    Dim        = Def.Dim;
elseif (nargin==4)
    Dim        = Def.Dim;
elseif (nargin==5)
    % do nothing
else
     errId = 'filter_fft1.m:TooManyInputArguments';
     errMsg = 'S,[Freq,CutLowFreq,Smooth,Dim]';
     error(errId, errMsg);
end

Len = size(S,Dim);
SX = (1:1:Len);
if (Dim==1)
    SX = SX.';
end

Taper = fermi_fun(SX,Freq.*Len,Smooth.*Len,~CutLowFreq);
CS    = abs(ifft(bsxfun(@times,fft(S,[],Dim),Taper)));

