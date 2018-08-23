function Rot=rotm(p,RotType)
% Return a numeric or symbolic 3-D rotation matrix about the X, Y or Z axis
% Package: Util.Geom
% Description: Return a numeric or symbolic 3-D rotation matrix about the
%              X, Y or Z axis.
% Input  : - If numeric scalar than this is the rotation angle [radians].
%            Alternatively if this is a string then the function will
%            return a symbolic rotation matrix and the string will be the
%            symbolic angle to use in the rotation matrix.
%          - Rotation matrix type - options are:
%            {1|'x'} - for rotation about the X axis.
%            {2|'y'} - for rotation about the Y axis.
%            {3|'z'} - for rotation about the Z axis.
%            {'sx'}  - Return symbolic rotation matrix about the X axis.
%            {'sy'}  - Return symbolic rotation matrix about the Y axis.
%            {'sz'}  - Return symbolic rotation matrix about the Z axis.
% Output : - Rotation matrix.
% Tested : Matlab 5.0
%     By : Eran O. Ofek                    Feb 1994
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Rot=rotm(45./RAD,'x');   % 3-D rotation 40deg about z-axis.
%          Rot=rotm(45./RAD,'Z'); Rot(1:2,1:2); % 2-D rotation 45deg.
%          Rot=rotm('Theta','x');  % symbolic rotation matrix
% Reliable: 1
%--------------------------------------------------------------------------

if (nargin==~2),
   error('Illegal number of input arguments');
end

if (isnumeric(p)),
   % numeric rotation matrix
   switch RotType
    case {1,'x','X'}
       Rot = [1, 0, 0; 0, cos(p), -sin(p); 0, sin(p), cos(p)];
    case {2,'y','Y'}
       Rot = [cos(p), 0, sin(p); 0, 1, 0; -sin(p), 0, cos(p)];
    case {3,'z','Z'}
       Rot = [cos(p), -sin(p), 0; sin(p), cos(p), 0; 0, 0, 1];
    otherwise
       error('Unknown RotType option');
   end

elseif (ischar(p)),
   % symbolic rotation matrix
   eval(sprintf('syms %s',p));
   syms Rot;

   switch RotType
    case {1,'x','X'}
       eval(sprintf('Rot = [1, 0, 0; 0, cos(%s), -sin(%s); 0, sin(%s), cos(%s)];',p,p,p,p));
    case {2,'y','Y'}
       eval(sprintf('Rot = [cos(%s), 0, sin(%s); 0, 1, 0; -sin(%s), 0, cos(%s)];',p,p,p,p));
    case {3,'z','Z'}
       eval(sprintf('Rot = [cos(%s), -sin(%s), 0; sin(%s), cos(%s), 0; 0, 0, 1];',p,p,p,p));
    otherwise
       error('Unknown RotType option');
   end

else
   error('First input argument must be numeric or string');
end
