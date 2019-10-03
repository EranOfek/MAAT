function [S,Sflux,IsLocalMax]=filter1_wbank(V,F)
% Filter a 1D data vector with a templates bank (multiple filters).
% Package: Util.filter
% Description: Filter a 1D data vector with a templates bank (multiple
%              filters). If required the templates are pad with zeros to
%              have the size of the data vector.
% Input  : - data vector.
%          - Matrix of templates bank, in which each column corresponds to
%            a template.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - Matrix of scores. Each column corresponds to a template, and
%            each row for a time.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [S]=Util.filter.filter1_wbank(V,F)
% Reliable: 
%--------------------------------------------------------------------------

%[Mag,Res]=AstroUtil.microlensing.microlens_ps([0 0.1 0.1 1 0],t);
%F=Res.Mu-1;
%[Mag,Res]=AstroUtil.microlensing.microlens_ps([0 0.2 0.1 1 0],t);
%F(:,2)=Res.Mu-1;   

%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% columns in F contains filters
% normalize filters
NormF = sum(F);
F = F./NormF;
% calculate sqrt(F^2) after normalization
NormF2 = sqrt(sum(F.^2));
RStdV = Util.stat.rstd(V);

NV = numel(V);
[NF,~] = size(F);

PadSize = (NV - NF)./2;
if (PadSize==floor(PadSize))
    PadSize1 = PadSize;
    PadSize2 = PadSize;
else
    PadSize1 = floor(PadSize);
    PadSize2 = ceil(PadSize);
end


F = padarray(F,PadSize1,0,'pre');
F = padarray(F,PadSize2,0,'post');

F = fftshift(F,1);

S = ifft(fft(V,[],1).*conj(fft(F,[],1)));
Sflux = S./(NormF2.^2).*max(F);
S = S./(NormF2.*RStdV);

if (nargout>2)
  IsLocalMax = islocalmax(S);
end


