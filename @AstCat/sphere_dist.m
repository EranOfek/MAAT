function [Dist]=sphere_dist(AstC,Long,Lat,ColNames,CatUnits)
% Spherical distance between a point and the positions in an AstCat object.
% Package: @AstCat
% Description: Calculate the spherical distance between a point and the
%              position of all sources in an AstCat object.
% Input  : - AstCat object.
%          - Longitude position of point
%            [radians, sexagesimal string or [H M S].
%          - Latitude position of point
%            [radians, sexagesimal string or [Sign D M S].
%          - Cell array of column names containing the X and Y position
%            in the AstCat object. Default is {'XWIN_IMAGE','YWIN_IMAGE'}.
%          - Units of the AstCat object catalog {'deg'|'rad'}.
%            If empty will try to read the units from the AstCat.ColUnits
%            If this field is empty then will set units to 'deg'.
%            Default is empty.
% Output : - A structure array (element per AstCat element) containing
%            .Dist - the distance [radians].
%            .PA   - The position angle [radians].
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Dist]=sphere_dist(AstC,1,1);
% Reliable: 2
%--------------------------------------------------------------------------
InvRAD = pi./180;

CatField = 'Cat';

Def.CatUnits = 'deg';
if (nargin==3),
    ColNames = {'ALPHAWIN_J2000','DELTAWIN_J2000'};
    CatUnits = [];
elseif (nargin==4),
    CatUnits = [];
elseif (nargin==5),
    % do nothing
else
    error('Illegal number of input arguments');
end

% convert Long/Lat to radians:
Long = celestial.coo.convertdms(Long,'gH','r');
Lat  = celestial.coo.convertdms(Lat,'gD','R');

Nc = numel(AstC);
Dist = struct('Dist',cell(size(AstC)),'PA',cell(size(AstC)));
for Ic=1:1:Nc,
    % for each catalog
    ColInd = colname2ind(AstC(Ic),ColNames);
    
    XY = AstC(Ic).(CatField)(:,ColInd);
    if (istable(XY)),
        XY = table2array(XY);
    end
    
    if (isempty(CatUnits)),
        % try to read units from AstCat
        if (~isempty(AstC(Ic).ColUnits)),
            CatUnits = AstC(Ic).ColUnits{ColInd(1)};
        else
            CatUnits = Def.CatUnits;
        end
    end
    
    switch lower(CatUnits)
        case {'rad','radian'}
            % do nothing
        case 'deg'
            XY = XY.*InvRAD;
        otherwise
            error('Unknown CatUnits option');
    end
            
    [Dist(Ic).Dist,Dist(Ic).PA] = celestial.coo.sphere_dist(Long,Lat,XY(:,1),XY(:,2));
    
end

