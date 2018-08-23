function Mat=load_htmcat(BaseName,HTMind,varargin)
%--------------------------------------------------------------------------
% load_htmcat function                                           Catalogue
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


NumVarargs = length(varargin);
if NumVarargs > 3
     errId = 'load_htmcat.m:TooManyInputArguments';
     errMsg = 'BaseName, HTMind, [GroupSize, FileType, Path]';
     error(errId, errMsg);
end
Gaps = cellfun(@isempty, varargin);
DefArgs = {100, ''};    % default input arguments
Suboptargs = DefArgs(1 : NumVarargs);
varargin(Gaps) = Suboptargs(Gaps);
DefArgs(1 : NumVarargs) = varargin;
[GroupSize, FileType, Path] = DefArgs{:};


PWD = pwd;
if (~isempty(Path)),
    cd(Path);
end

Nhtm = numel(HTMind);
for Ihtm=1:1:Nhtm,
    switch lower(FileType)
        case 'mat'
            FileName = sprintf('%s%
            load(


