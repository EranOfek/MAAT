function [S]=frt(S,varargin)
% Applay the Fast Radon Transform to a set of images.
% Package: @SIM
% Description: Subtract image background and applay the Fast Radon
%              Transform to a set of images.
%              Either on the image or the image transpose.
% Input  : - A SIM object containing images.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Trans' - Applay to image transpose. Default is false.
%            'SubBack' - Subtract background. Default is true.
%            'BackPar' - Cell array of additional parameters to pass to the
%                      background method. Default is {}.
% Output : - A SIM object in which the Im field is replaced with the Radon
%            transform of the image.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Apr 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S=frt(S);
% Reliable: 
%--------------------------------------------------------------------------

ImageField = SIM.ImageField;

DefV.Trans                = false;
DefV.SubBack              = true;
DefV.BackPar              = {};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% subtract background
if (InPar.SubBack)
    S = background(S,InPar.BackPar{:});
    S = sub_background(S);
end

% execute FRT
Nim = numel(S);
for Iim=1:1:Nim
    S(Iim).(ImageField) = radon.frt(S(Iim).(ImageField),'Trans',InPar.Trans);
end