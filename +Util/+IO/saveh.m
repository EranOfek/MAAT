function VarName=saveh(FileName,Data,VarName,AttribCell,varargin)
% Save a matrix into HDF5 file.
% Package: Util.IO
% Description: Save a matrix into HDF5 file. If file exist then will
%              add it as a new dataset to the same file.
%              This is becoming faster than matlab (2014a) for matices with
%              more than ~10^5 elements.
% Input  : - File name.
%          - Data to save (e.g., matrix).
%          - Dataset (variable) name under which to store the data
%            in the HDF5 file. If varaible already exist, then function
%            will fail. Default is 'V'.
%            If empty use 'V'.
%          - Cell array of attributes {key,val,...}.
%            Default is {}.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            to pass to the h5create. See h5create.m for details.
% Output : - Variable name (dataset) in which the matrix was stored.
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    May 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: VarName=Util.IO.saveh('R.hd5',T,'G/A');
% Reliable: 2
%--------------------------------------------------------------------------

Def.VarName = 'V';
Def.AttribCell = {};
if (nargin==2)
    VarName    = Def.VarName;
    AttribCell = Def.AttribCell;
elseif (nargin==3)
    AttribCell = Def.AttribCell;
else
    % do nothing
end
if (isempty(VarName))
    VarName = Def.VarName;
end

DataSet = sprintf('/%s',VarName);
if (isnumeric(Data))
    % save numeric data  
    h5create(FileName,sprintf('/%s',VarName),size(Data),varargin{:},'Datatype',class(Data)); %,'Deflate',1);
    h5write(FileName,DataSet,Data);
else
    error('Non numeric datatype are not supheported yet');
end

if (~isempty(AttribCell))
    % populate attributes
    Natt = numel(AttribCell);
    for Iatt=1:2:Natt-1
        h5writeatt(FileName,sprintf('/%s',VarName),AttribCell{Iatt},AttribCell{Iatt+1});
    end
end
