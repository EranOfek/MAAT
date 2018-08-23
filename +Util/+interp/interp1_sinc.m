function Yxx=interp1_sinc(Y,XX)
% 1-D sinc interpolation
% Package: Util.interp
% Description: Interpolation of a 1-D array using the Whittaker–Shannon
%              interpolation formula (i.e., sinc interpolation).
% Input  : - Coloumn vector of values sampled at equal intervals.
%            The index (i.e., position) of this values is assumed to be
%            1 to N, where N is the vector length.
%          - Vector of positions in which to perform the sinc
%            interpolation.
% Output : - Interpolated values.
% Reference: https://en.wikipedia.org/wiki/Weighted_median
% License: GNU general public license version 3
% Tested : Matlab R2015a
%     By : Eran O. Ofek                    Jun 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: X=(1:1:100)'; Y=sin(2.*pi.*X./10); Yxx=Util.interp.interp1_sinc(Y,[2.1 3.4])
%          XX = [1:0.1:100]; Yxx=Util.interp.interp1_sinc(Y,XX)
% Reliable: 2
%--------------------------------------------------------------------------


N    = size(Y,1);
VecN = (1:1:N).'; 
%sum(Y(round(XX-VecN)));

% XX need to be a row vector - hence XX(:).'
XmN = bsxfun(@minus,XX(:).',VecN);

Yxx = sum(bsxfun(@times,Y,sinc(XmN)),1);

% alternatively: (but slower): Yxx = Y'*sinc(XmN);

