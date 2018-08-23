function [Color]=value2color(Vec,ColorMap)
%--------------------------------------------------------------------------
% value2color function                                            plotting
% Description: 
% Input  : - A coloumn vector of values.
%          - A color map. Default is 'jet'.
% Output : - A color corresponding to the value, by a linear
%            transformation between the ...
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Feb 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------
InterpMethod = 'linear';

Def.ColorMap = 'jet';

if (nargin==1),
    ColorMap = Def.ColorMap;
elseif (nargin==2),
    % do nothing
else
    error('Illegal number of input arguments');
end

CM = colormap(ColorMap);
MinVec = min(Vec);
MaxVec = max(Vec);
LinVec = (MinVec:(MaxVec-MinVec)./(size(CM,1)-1):MaxVec).';

% interpolate RGB
R = interp1(LinVec,CM(:,1).',Vec,InterpMethod);
G = interp1(LinVec,CM(:,2).',Vec,InterpMethod);
B = interp1(LinVec,CM(:,3).',Vec,InterpMethod);

Color = [R,G,B];