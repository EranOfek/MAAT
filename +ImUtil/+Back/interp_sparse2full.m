function [IntImage]=interp_sparse2full(X,Y,Image,OutSize,varargin)
% Interpolate sparse image to a full image
% Package: ImUtil.Back
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


DefV.Method               = 'interp';   % 'interp' | 'fit'
DefV.InterpMethod         = 'makima';   % 'linear' | 'nearest' | 'spline' | 'cubic' | 'makima'

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

[MatX,MatY] = meshgrid((1:1:OutSize(2)),(1:1:OutSize(1)));

switch lower(InPar.Method)
    case 'interp'
        
        if (any(size(Image)==1))
            % assume X,Y,Image are vectors and represent a scattered image
            % in this case use scattered interpolation
            
            
        
        IntImage = interp2(X,Y,Image,MatX,MatY,InPar.InterpMethod);
        
    case 'fit'
        
    otherwise
        error('Unknown Method option');
end
