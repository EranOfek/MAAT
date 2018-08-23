function Flag=in_halfspace(R,N,C,Dim,Crit)
% Is point in half space
% Package: celestial.htm
% Description: Given a unit vector R and half space (N,C) test if the
%              point is contained inside the half space (N dot R > C).
%              A halfspace is a plane that splits the sphere in two.
%              It is defined by a direction vector (N) and a signed scalar
%              (C), measured along the normal vector from the origin of
%              the sphere.
% Input  : - A three columns matrix in which each row is a unit vector (R).
%            Alternatively, by setting Dim=2, this can be a three rows
%            matrix in which each column is a unit vector (R).
%            If R is a two column matrix, then assume the columns are
%            [Long, Lat] in radians.
%            If Dim=2, then the input should be a two rows matrix with
%            [Long, Lat] in each column.
%          - A unit vector (3 elements) or [Long, Lat] coordinates
%            (two elements) specifing the half space position vector (N).
%            Alternatively, this can be a matrix which dimensions are
%            identical to the input position vectors (R). In this case,
%            each vector R will be compared with the corresponding vector
%            (N).
%          - A signed scalar (C) between -1 and 1, measured along the
%            normal vector from the origin of the sphere.
%            Alternatively, C may be a vector which length is identical
%            to the number of vectors in (R).
%            (C) is the cosine of the angular radius of the cap.
%          - Dimension along to operate.
%            If 1, then assume the position vectors in (R) and (C) are
%            in rows, while if 2, then assume they are in columns.
%          - Flag indicating if to use ">" or ">=".
%            If true then use (N dot R > C),
%            If flase then use (N dot R) >= C 
% Output : - A flag indicating if the vector is inside the halfspace (true)
%            or on or outside the halfspace (false).
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    Jul 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Flag=celestial.htm.in_halfspace([1 1 1],[1 1 1],0.5,1,1)
% Reliable: 2
%--------------------------------------------------------------------------

% Def.Dim  = 1;
% Def.Crit = true;
% if (nargin==3),
%     Dim  = Def.Dim;
%     Crit = Def.Crit; 
% elseif (nargin==4),
%     Crit = Def.Crit; 
% elseif (nargin==5),
%     % do nothing
% else
%     error('Illegal number of input arguments');
% end

if (Dim==2)
    R = R.';
    N = N.';
end

[~,Qr]= size(R);
[~,Qn]= size(N);
     
if (Qr==2)
   % convert [Long, Lat] to cosine directions
   [CD1, CD2, CD3] = celestial.coo.coo2cosined(R(:,1),R(:,2));
   R = [CD1, CD2, CD3];
end

if (Qn==2)
   % convert [Long, Lat] to cosine directions
   [CD1, CD2, CD3] = celestial.coo.coo2cosined(N(:,1),N(:,2));
   N = [CD1, CD2, CD3]; 
end


%Flag = dot(N,R,2)>C;   % slower


% if (Crit),
%    Flag = (N(:,1).*R(:,1) + N(:,2).*R(:,2) + N(:,3).*R(:,3))>C;
%    %Flag = sum(N.*R,2)>C; % somewhat slower
% else
%    Flag = (N(:,1).*R(:,1) + N(:,2).*R(:,2) + N(:,3).*R(:,3))>=C;
% end


Flag = [bsxfun(@times,N(:,1),R(:,1).') + bsxfun(@times,N(:,2),R(:,2).') + bsxfun(@times,N(:,3),R(:,3).')>C];
