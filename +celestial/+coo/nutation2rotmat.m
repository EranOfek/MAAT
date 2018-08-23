function NutMatrix=nutation2rotmat(Nut,JD,MatType)
%--------------------------------------------------------------------------
% nutation2rotmat function                                           ephem
% Description: Given nutation in longitude and obliquity (in radians)
%              and JD, return the Nutation rotation matrix.
% Input  : - A two column matrix of nutation in [Long, Obliq].
%          - JD for which the nutations were calculated.
%          - Type of nutation matrix:
%            'f' : full precision, (default).
%            'l' : linearized matrix.
% Output : - A cube of Nutation rotation matrices. The 3rd Dim is for
%            various times.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jun 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: NutMatrix=nutation2rotmat(Nut,JD);
% Reliable: 2
%--------------------------------------------------------------------------

MatType = 'f';
if (nargin==2),
    MatType = Def.MatType;
end

DLon = Nut(:,1);
DObl = Nut(:,2);
NutMatrix = ones(3,3,length(JD));
Obl = obliquity(JD);
switch lower(MatType)
    case 'f'
       % Full precesion Nutation Matrix
       Obl1 = Obl + DObl;
       NutMatrix(1,1,:) = cos(DLon);
       NutMatrix(1,2,:) = -sin(DLon).*cos(Obl);
       NutMatrix(1,3,:) = -sin(DLon).*sin(Obl);
       NutMatrix(2,1,:) = sin(DLon).*cos(Obl1);
       NutMatrix(2,2,:) = cos(DLon).*cos(Obl1).*cos(Obl) + sin(Obl1).*sin(Obl);
       NutMatrix(2,3,:) = cos(DLon).*cos(Obl1).*sin(Obl) - sin(Obl1).*cos(Obl);
       NutMatrix(3,1,:) = sin(DLon).*sin(Obl1);
       NutMatrix(3,2,:) = cos(DLon).*sin(Obl1).*cos(Obl) - cos(Obl1).*sin(Obl);
       NutMatrix(3,3,:) = cos(DLon).*sin(Obl1).*sin(Obl) + cos(Obl1).*cos(Obl);
    case 'l'
       % Linearized Nutation matrix
       %NutMatrix(1,1,:) = 1;
       NutMatrix(1,2,:) = -DLon.*cos(Obl);
       NutMatrix(1,3,:) = -DLon.*sin(Obl);
       NutMatrix(2,1,:) = DLon.*cos(Obl);
       %NutMatrix(2,2,:) = 1;
       NutMatrix(2,3,:) = -DObl; 
       NutMatrix(3,1,:) = DLon.*sin(Obl);
       NutMatrix(3,2,:) = DObl;
       %NutMatrix(3,3,:) = 1;
    otherwise
       error('Unknown MatType option');
end


