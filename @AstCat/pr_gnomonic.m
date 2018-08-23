function [X,Y]=pr_gnomonic(AstC,ColLongLat,Scale,CenCoo,Units)
%--------------------------------------------------------------------------
% pr_gnomonic function                                       class/@AstCat
% Description: Project the long/lat coordinates in an AstCat
%              class object using gnomonic projection.
% Input  : - A single AstCat object.
%          - Vector of indices of the Long and Lat columns.
%            Alternativelly, a cell array of column names
%            (e.g., {'RA','Dec'})
%          - Scale radius. Default is 1.
%            For example, in order to convert from input
%            coordinates in radians, to ouput coordinates
%            with 1 arcsec scale use RAD.*3600.
%          - Projection central coordinates [Lon,Lat]
%            in radians.
%          - Units of input coordinates ('rad'|'deg').
%            The program will first attempt to read the units
%            from the ColUnits field. If not available than
%            will use this argument. Default is 'rad'.
%            Long and Lat should have the same units.
% Output : - Projected X coordinates.
%          - Projected Y coordinates.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [X,Y] = pr_gnomonic(AstC,{'RA','Dec'},RAD.*3600,[0 0],'rad')
% Reliable: 2
%--------------------------------------------------------------------------

import celestial.proj.*

            
if (numel(AstC)>1),
    error('pr_gnomonic input should be an AstCat with a single element');
end

if (nargin<5),
    Units = 'rad';
end

[ColLong, ColLat] = colname2ind(AstC,ColLongLat);

% check units of input
if (~isempty(AstC.ColUnits)),
    Units = AstC.ColUnits{ColLong};
end
switch lower(Units)
    case {'rad','radian'}
        ConvFactor = 1;
    case 'deg'
        ConvFactor = pi./180;
    otherwise
        error('Unknown Units option');
end

%what is going on here????  AstC(:,ColLong) doesn't work!!
if istable(AstC.Cat(:,ColLong)),
    Long = table2array(AstC.Cat(:,ColLong));
    Lat  = table2array(AstC.Cat(:,ColLat));
else
    Long = AstC.Cat(:,ColLong);
    Lat  = AstC.Cat(:,ColLat);
end
[X,Y] = pr_gnomonic(Long.*ConvFactor,Lat.*ConvFactor,Scale,CenCoo);

