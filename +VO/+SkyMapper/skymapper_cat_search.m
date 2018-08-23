function [Out,Link]=skymapper_cat_search(RA,Dec,varargin)
% Cone search the SkyMapper online catalog.
% Package: VO.SkyMapper
% Description: Cone search the SkyMapper online catalog.
% Input  : - Vector of J2000.0 R.A. [rad] or see celestial.coo.convertdms.
%          - Vector of J2000.0 Dec. [rad] or see celestial.coo.convertdms.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Radius'   - Search radius. Default is 30 [arcsec].
%            'RadiusUnits' - Search radius units. Default is 'arcsec'.
%            'Output'   - Catalog output format:
%                         'astcat' - An AstCat object (matrix catalog).
%                         'astcat_t' - An AstCat object (table catalog).
%                         'struct' - A structure array (matrix catalog).
%                         'struct_t' - A structure array (tanle catalog).
%            'BaseURL' - Base URL search. Default is:
%                        'http://skymapper.anu.edu.au/sm-cone/query?'.
% Output : - A structure or AstCat array in which each element is for
%            one search coordinate, and contains a table of entries within
%            the radius from the search coordinates.
%          - A cell array of links from which the data was retrieved.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Dec 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference: http://skymapper.anu.edu.au/how-to-access/#adql
% Example: [Out,Link]=VO.SkyMapper.skymapper_cat_search(1,-1)
% Reliable: 2
%--------------------------------------------------------------------------

CatField     = AstCat.CatField;
ColField     = AstCat.ColField;
ColCellField = AstCat.ColCellField;

DefV.Radius               = 30;
DefV.RadiusUnits          = 'arcsec';
DefV.Output               = 'astcat';   % 'astcat' | 'astcat_t' | 'struct' | 'struct_t'
DefV.BaseURL              = 'http://skymapper.anu.edu.au/sm-cone/query?';
DefV.TimeOut              = 30;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


RA  = celestial.coo.convertdms(RA,'gH','d');
Dec = celestial.coo.convertdms(Dec,'gD','d');
Ncoo = numel(RA);
Radius = convert.angular(InPar.RadiusUnits,'deg',InPar.Radius);
if (numel(Radius)==1)
    Radius = Radius.*ones(Ncoo,1);
end

RadiusStr = sprintf('SR=%f',Radius);
FormatStr = sprintf('RESPONSEFORMAT=CSV');
VerbStr   = sprintf('VERB=3');



Link = cell(Ncoo,1);
switch lower(InPar.Output)
    case {'astcat','astcat_t'}
        % define AstCat object
        Out = AstCat(Ncoo,1);
    case {'struct','struct_t'}
        % define structure
        Out = Util.struct.struct_def({CatField,ColCellField,ColField},Ncoo,1);
    otherwise
        error('Unknown Output option');
end
for Icoo=1:1:Ncoo
    PosStr     = sprintf('RA=%f&DEC=%f',RA(Icoo),Dec(Icoo));
    Link{Icoo} = sprintf('%s%s&%s&%s&%s',InPar.BaseURL,...
                    PosStr,RadiusStr,...
                    FormatStr, VerbStr);
                
    Options    = weboptions('timeout',InPar.TimeOut);
    Table      = webread(Link{Icoo},Options);
    
    switch lower(InPar.Output)
        case {'astcat','struct'}
            % save data as matrix
            % first delete object name (string)
            Table.smss_j = [];
        otherwise
            % do nothing
    end
    if (Icoo==1)
        ColCell    = Table.Properties.VariableNames;
    end
    
    Out(Icoo).(CatField)     = Table;
    Out(Icoo).(ColCellField) = ColCell;
    Out(Icoo).(ColField)     = cell2struct(num2cell(1:1:numel(Table.Properties.VariableNames)),Table.Properties.VariableNames,2);
end
    
   