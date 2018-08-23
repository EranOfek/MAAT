function [Cat,ColCell,ColUnits,Result]=cds_astcat_search(CatName,RA,Dec,varargin)
% Query a VizieR catalog using the cdsclient tools
% Package: VO.VizieR
% Description: Query a VizieR catalog using the cdsclient command line
%              tools. Allow to query a specific catalog by coordinates, or
%              object name and return the output in various formats.
%              Catalog.VizieR.cdsclient_prog_names can be used to identify
%              existing catalog excess programs.
%              In order to execuate a query you have to verify that the
%              catalog mapping exist in Catalog.VizieR.catalog.mapping.
% Input  : - Catalog name. E.g, 'finducac4'. See
%            Catalog.VizieR.cdsclient_prog_names for available catalog names.
%          - J2000.0 R.A. or object name.
%            If the third argument is empty then this is an object name.
%            RA and Dec Units are defined by 'CooUnits' argument.
%          - J2000.0 Dec.
%          * Arbitary number of pairs of ...,keyword,value,... arguments.
%            Possible keywords are:
%            'OutType' - Options are:
%                        'AstCat' - An AstCat object.
%                        'AstCatTable' - An AstCat object with a table
%                                   catalog.
%                        'cell' - Cell array.
%                        'table' - A Table.
%                        'mat' - A matrix.
%                        Default is 'AstCat'.
%            'MapName' - Mapping name. Default is empty.
%            'RemoveStringCol' - Remove all columns with string format from
%                        output catalog. This can be useful to present the
%                        output in an Matrix or AstCat formats.
%                        Default is true.
%            'Replace999' - Replace 999 with NaN. Default is true.
%            'CooUnits'- Coordinate units. Default is 'deg'.
%            'Radius'  - Search radius. Default is 10.
%                        For box search this can be [X,Y] box size.
%            'RadiusUnits' - Search radius units. Default is 'arcmin'.
%                        See convert.angular for options.
%                        Default is 'arcmin'.
%            'ObjName' - Object name. If provided, then will search by
%                        object name. Default is empty.
%            'RegionType' - Search region type: 'circ'|'box'.
%                        Default is 'circ'.
%            'MaxRecord' - Max. number of records. Default is 100000.
%            'PathPar'   - Additional arguments to pass to
%                          Catalog.VizieR.cdsclient_path.
%                          Default is {}.
%            'ExpandTabs'- Exapnd tabs to spaces. Default is true.
% Output : - Output catalog
%          - Cell array of column names.
%          - Cell array of column units.
%          - The complete result returned from CDS in a string format.
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cat=VO.VizieR.cds_astcat_search('finducac4',1,1,'OutType','table')
%          Cat=VO.VizieR.cds_astcat_search('finducac4',1,1,'OutType','astcat','CooUnits','rad')
% Reliable: 2

import VO.VizieR.*

DefV.OutType             = 'AstCat';  % AstCat | AstCatTable | cell | table | mat
DefV.MapName             = [];
DefV.RemoveStringCol     = true;
DefV.Replace999          = true;    % repalce 999 with NaN
DefV.CooUnits            = 'deg';
DefV.Radius              = 10;   % or [x,y]
DefV.RadiusUnits         = 'arcmin';
DefV.ObjName             = '';
DefV.RegionType          = 'circ';  % 'circ' | 'box'
DefV.MaxRecord           = 100000;
DefV.PathPar             = {};
DefV.ExpandTabs          = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% construct cdsclient search command
if (isempty(Dec))
    Command = construct_vizquery('CatName',CatName,'ObjName',RA,...
                                 'CooUnits',InPar.CooUnits,...
                                 'Radius',InPar.Radius,...
                                 'RadiusUnits',InPar.RadiusUnits,...
                                 'RegionType',InPar.RegionType,...
                                 'MaxRecord',InPar.MaxRecord,...
                                 'PathPar',InPar.PathPar);
else
    Command = construct_vizquery('CatName',CatName,'RA',RA,'Dec',Dec,...
                                 'CooUnits',InPar.CooUnits,...
                                 'Radius',InPar.Radius,...
                                 'RadiusUnits',InPar.RadiusUnits,...
                                 'RegionType',InPar.RegionType,...
                                 'MaxRecord',InPar.MaxRecord,...
                                 'PathPar',InPar.PathPar);
        
end
%Command
% execuate command
[Status,Result] = system(Command);
%Result

if (Status<0)
    % execuation failed
    fprintf('--- Result of failed execuation ---\n');
    fprintf('%s\n',Result);
    error('VizieR cdsclient command execuation failed');
end

% get catalog mapping
Map = catalog_mapping(CatName);
% If more than one map exist - choose first
Map = Map(1);

% Read catalog
Ncol = numel(Map.Cols);
Format = cell(Ncol,3);
for Icol=1:1:Ncol
    Format(Icol,:) = {Map.Pos{Icol}(1), Map.Pos{Icol}(2), '%s'};
end
    

if (InPar.ExpandTabs)
    TempFile = tempname;
    FID = fopen(TempFile,'w');
    fprintf(FID,'%s',Result);
    fclose(FID);
    [~,Result] = system(sprintf('expand %s',TempFile));
    delete(TempFile);
end

% note that the 'skip' 3 is required because of some wget errors!
% need to verify this?!

Data = Util.string.read_str_formatted(Result,Format,'comment','#','skip',3);
for Icol=1:1:Ncol
    switch Map.Format{Icol}
        case '%s'
            % do nothing
        case '%f'
            Data(:,Icol) = num2cell(cellfun(@str2double,Data(:,Icol)));
            if (InPar.Replace999)
                Flag999 = cell2mat(Data(:,Icol))==999;
                Data(Flag999,Icol) = {NaN};
            end
        otherwise
            error('Unknown Mapping format');
    end
end

if (InPar.RemoveStringCol)
    % remove columns containing strings
    Flag = ~strcmp(Map.Format,'%s');
    Data = Data(:,Flag);
    Map.Cols  = Map.Cols(Flag);
    Map.Units = Map.Units(Flag);
end

ColCell  = Map.Cols;
ColUnits = Map.Units;
switch lower(InPar.OutType)
    case 'cell'
        % do nothing - already in cell format
        Cat = Data;
    case 'mat'
        Cat = cell2mat(Data);
    case 'table'
        Cat = cell2table(Data,'VariableNames',ColCell);
        Cat.Properties.VariableUnits = ColUnits;
    case 'astcat'
        Cat = AstCat;
        Cat.(AstCat.CatField) = cell2mat(Data);
        Cat.(AstCat.ColCellField) = ColCell;
        Cat = colcell2col(Cat);
        Cat.ColUnits = ColUnits;
    case 'astcattable'
        Data = cell2table(Data,'VariableNames',ColCell);
        Data.Properties.VariableUnits = ColUnits;
        Cat = AstCat;
        Cat.(AstCat.CatField) = Data;
        Cat.(AstCat.ColCellField) = ColCell;
        Cat = colcell2col(Cat);
        Cat.ColUnits = ColUnits;
    otherwise
        error('Unknown OutType option');
end

