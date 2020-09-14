function [Mode,Std]=mode_density(Array,Dim,varargin)
% Calculate the mode by estimating the density of points
% Package: Util.stat
% Description: Estimate the mode by constructing equal number of points
%              histogram and choosing the bin with smallest range.
% Input  : - Array.
%          - Dimesnion over which to calculate the mode.
%            1,2,'all'.
%            Default is 'all'.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'NinBin' - Number of points in bin. Default is 30.
%            'JoinOut' - Join output parameters into one variable.
%                       Default is false.
% Output : - Mode estimator.
%          - Estimator for the Std of a Gaussian peak.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Mode]=Util.stat.mode_density(rand(10000,1),1)
% Reliable: 
%--------------------------------------------------------------------------


if nargin<2
    Dim = 'all';
end

DefV.NinBin               = 30;
DefV.JoinOut              = false;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if ischar(Dim)
    % calc mode over all dimensions
    Array = Array(:);
    Dim   = 1;
end

SortedArray = sort(Array,Dim);
if (Dim==2)
    SortedArray = SortedArray.';
end
Pos    = SortedArray(1:InPar.NinBin:end,:);
Ranges = diff(Pos,1,1);

%[[1:1:numel(Ranges)]',Ranges]
[MinRange,MaxI] = min(Ranges,[],1);
Mode = Pos(MaxI,:) - MinRange.*0.5;

if (nargout>1 || InPar.JoinOut)
    % for 1-sigma look for density which is factor of 1/exp(-1/2))=1.6487 higher
    DD = Ranges - MinRange./exp(-1./2);
    I1 = find(DD<0,1,'first');
    I2 = find(DD<0,1,'last');
    Low  = Pos(I1) - Ranges(I1).*0.5;
    High = Pos(I2) + Ranges(I2).*0.5;
    Std  = (High-Low).*0.5;
end

if InPar.JoinOut
    Mode = [Mode, Std];
end


