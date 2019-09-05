function [IntImage]=interp_sparse2full(X,Y,Image,OutSize,varargin)
% Interpolate sparse image to a full image
% Package: ImUtil.Back
% Description: 
% Input  : - 
%          -
%          -
%          - Output size [X,Y]
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

[MatX,MatY] = meshgrid((1:1:OutSize(1)),(1:1:OutSize(2)));

switch lower(InPar.Method)
    case 'interp'
        
        if (any(size(Image)==1))
            % assume X,Y,Image are vectors and represent a scattered image
            % in this case use scattered interpolation
            warning('Switched to linear scattered interpolation');
            
            F = scatteredInterpolant(X,Y,Image,'linear');
            IntImage = F(MatX,MatY);
            
            
        else
            % regular grid - use interp2
            IntImage = interp2(X,Y,Image,MatX,MatY,InPar.InterpMethod);
        end
        
    case 'fit'
        error('fit option - not yet implemented');
    otherwise
        error('Unknown Method option');
end
