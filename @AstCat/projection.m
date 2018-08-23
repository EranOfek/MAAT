function [X,Y]=projection(AstC,ProjectionType,ColLongLat,Param,Units)
% Coordinate projection transformation for columns in AstCat object.
% Package: @AstCat
% Description: Apply coordinates projection transformation to coordinate
%              columns in an AstCat object.
% Input  : - A single element AstCat object.
%          - Projection type. One of the following:
%            'tan'|'gnomonic'  - Gnomonic projection.
%                       Parameters are: [Scale, Center_Long, Center_Lat].
%                       Default parameters: [1 0 0].
%            'itan'|'gnomonic' - Inverse Gnomonic projection.
%                       Parameters are: [Scale, Center_Long, Center_Lat].
%                       Default parameters: [1 0 0].
%          - Vector of indices of the Long and Lat (or X/Y) columns.
%            Alternativelly, a cell array of column names
%            (e.g., {'RA','Dec'})
%          - Vector of parameters to pass to the transformation function.
%            Units are always radians.
%          - Units of input coordinates ('rad'|'deg').
%            The program will first attempt to read the units
%            from the ColUnits field. If not available than
%            will use this argument. Default is 'rad'.
%            Long and Lat should have the same units.
% Output : - Vector of projected X (or Long) coordinates.
%            If only one output argument is requested than this is an
%            AstCat object in which the long/lat (or X/Y) columns contains
%            the projected coordinates.
%          - Vector of projected Y (or lat) coordinates.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [X,Y] = projection(AstC,'tan',{'RA','Dec'},[RAD.*3600,0 0],'rad')
% Reliable: 2
%--------------------------------------------------------------------------

import celestial.proj.*

CatField    = 'Cat';
CatColUnits = 'ColUnits';

DefParam.tan    = [1 0 0];
DefParam.itan   = [1 0 0];

            
if (numel(AstC)>1),
    error('projection input should be an AstCat with a single element');
end

if (nargin==3),
    Param = [];
    Units = 'rad';
elseif (nargin==4),
    Units = 'rad';
elseif (nargin==5),
    % do nothing
else
    error('Illegal number of input arguments');
end

[ColLong, ColLat] = colname2ind(AstC,ColLongLat);

% check units of input
if (~isempty(AstC.(CatColUnits))),
    Units = AstC.(CatColUnits){ColLong};
end
switch lower(Units)
    case {'rad','radian'}
        ConvFactor = 1;
    case 'deg'
        ConvFactor = pi./180;
    otherwise
        error('Unknown Units option');
end

if istable(AstC.(CatField)(:,ColLong)),
    Long = table2array(AstC.(CatField)(:,ColLong));
    Lat  = table2array(AstC.(CatField)(:,ColLat));
else
    Long = AstC.(CatField)(:,ColLong);
    Lat  = AstC.(CatField)(:,ColLat);
end

switch lower(ProjectionType)
    case {'tan','gnomonic'}
        if (isempty(Param)),
            Param = DefParam.tan;
        end
        [X,Y] = pr_gnomonic(Long.*ConvFactor,Lat.*ConvFactor, Param(1), Param(2:3) );
            
    case {'itan','ignomonic'}
        if (isempty(Param)),
            Param = DefParam.itan;
        end
        [X,Y] = pr_gnomonic(Long.*ConvFactor.*Param(1),Lat.*ConvFactor.*Param(1), Param(2:3) );
            
    otherwise
        error('Unknown projection type');
end


if (nargout<2),
    AstC.(CatField)(:,ColLong) = X;
    AstC.(CatField)(:,ColLat)  = Y;
    X = AstC;
end
    