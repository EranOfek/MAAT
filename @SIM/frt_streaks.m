function []=frt_streaks(S,varargin)
% Applay the Fast Radon Transform to a set of images and find streaks.
% Package: @SIM
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Apr 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

ImageField = SIM.ImageField;
PSFField   = 

DefV.SubBack              = true;
DefV.BackPar              = {};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% subtract background
if (InPar.SubBack)
    S = background(S,InPar.BackPar{:});
    S = sub_background(S);
end


Nim = numel(S);

for Iim=1:1:Nim
    F = radon.Finder;
    
    F.input_psf(S(Iim).(PSFField));
    F.input(S(Iim).(ImageField));
    
end
