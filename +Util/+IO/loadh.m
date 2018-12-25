function [Data,Att]=loadh(FileName,VarName,GetAtt)
% Load a matrix from HDF5 file. 
% Package: Util.IO
% Description: Load a matrix from HDF5 file. 
%              If dataset name is not provided than will read all
%              datasets into a structure. This function doesn't support
%              groups.
%              This is becoming faster than matlab (2014a) for matices with
%              more than ~10^4 elements.
% Input  : - File name.
%          - variable name (dataset). If empty or not provided than will
%            attempt to read all datasets.
%          - Get attribute. If empty then do not get attributes.
%            'h' - store attributes in an HEAD object. 
%            's' - store attributes in a structure array.
%            Default is empty.
% Output : - Datasets.
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    May 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Data=Util.IO.loadh('R.hd5','R');
% Reliable: 2
%--------------------------------------------------------------------------

Def.VarName = [];
Def.GetAtt  = [];
if (nargin==1)
    VarName = Def.VarName;
    GetAtt  = Def.GetAtt;
elseif (nargin==2)
    GetAtt  = Def.GetAtt;
else
    % do nothing
end

if (isempty(VarName))
    % check for available datasets in HD5 file
    Info = h5info(FileName);
    if (~isempty(Info.Datasets))
        Nds = numel(Info.Datasets);
        DatasetNames = Util.external.natsortfiles.natsort({Info.Datasets.Name}); % natural order sort (alphanumeric)
        for Ids=1:1:Nds
            Data.(DatasetNames{Ids}) = Util.IO.loadh(FileName,DatasetNames{Ids},GetAtt);
        end   
        Ind = (1:1:Nds);  % indices of all datasets
    end
else
    DataSet = sprintf('/%s',VarName);
    Data    = h5read(FileName,DataSet);
    
    if (GetAtt)
        Info = h5info(FileName);
        DatasetNames = Util.external.natsortfiles.natsort({Info.Datasets.Name}); % natural order sort (alphanumeric)
        Ind  = strcmp(DatasetNames,VarName);
    end
end

if (~isempty(GetAtt))
    % Get all attributes for datasets
    Nds = numel(Ind);
    [~, SortIdx] = Util.external.natsortfiles.natsort({Info.Datasets.Name}); % natural order sort (alphanumeric)
    switch lower(GetAtt)
        case 'h'
            % attribute will be stored in an HEAD object
            Att = HEAD(Nds,1);
            for Ids=1:1:Nds
                Att(Ids).Header = [{Info.Datasets(SortIdx(Ids)).Attributes.Name}.', {Info.Datasets(SortIdx(Ids)).Attributes.Value}.'];
            end          
        case 's'
            % attributes will be stored in a structure
            for Ids=1:1:Nds
                Att(Ids) = cell2struct({Info.Datasets(SortIdx(Ids)).Attributes.Value},{Info.Datasets(SortIdx(Ids)).Attributes.Name},2);
            end
        otherwise
            error('Unknown attributes output type');
    end
end

    