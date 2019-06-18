function fun_template(FunName,ToolBox,Path)
% Generate a function file template
% Pacakge: Util.code
% Description: Generate a functionm template with help section and basic
%              optional commands.
% Input  : - Function name (e.g., 'my_fun1.m').
%          - Package name (e.g., '+Util/+string/').
%          - Toolbox path. Default is '~/matlab/fun/'.
% Output : null
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: fun_template('trya1.m','+Util/+code/');
% Reliable: 2
%--------------------------------------------------------------------------


if (nargin==2)
    Path = '~/matlab/MAAT/';
end
FullPath = sprintf('%s%s%s%s',Path,ToolBox,filesep,FunName);

ReSplit = regexp(strrep(ToolBox,'+',''),filesep,'split');
Nsp = numel(ReSplit);
PackageName = ReSplit{1};
for Isp=2:1:Nsp-1
    PackageName = sprintf('%s.%s',PackageName,ReSplit{Isp});
end


FID = fopen(FullPath,'w');

fprintf(FID,'function []=%s(varargin)\n',FunName(1:end-2));
fprintf(FID,'%% SHORT DESCRIPTION HERE\n');
fprintf(FID,'%% Package: %s\n',PackageName);
fprintf(FID,'%% Description: \n');
fprintf(FID,'%% Input  : - \n');
fprintf(FID,'%%          * Arbitrary number of pairs of arguments: ...,keyword,value,...\n');
fprintf(FID,'%%            where keyword are one of the followings:\n');
fprintf(FID,'%% Output : - \n');
Version = regexp(version,'\(','split');
Version = Version{2}(1:end-1);
fprintf(FID,'%% License: GNU general public license version 3\n');
Date   = date;
Month  = Date(4:6);
Year   = Date(8:11);
AuN = sprintf('By : Eran O. Ofek');
fprintf(FID,'%%     %s                    %s %s\n',AuN,Month,Year);
fprintf(FID,'%%    URL : http://weizmann.ac.il/home/eofek/matlab/\n');
fprintf(FID,'%% Example: \n');
fprintf(FID,'%% Reliable: \n');
fprintf(FID,'%%--------------------------------------------------------------------------\n');
fprintf(FID,'\n\n');
% input args
fprintf(FID,'NumVarargs = length(varargin);\n');
fprintf(FID,'if NumVarargs > 3\n');
fprintf(FID,'     errId = ''%s:TooManyInputArguments'';\n',FunName);
fprintf(FID,'     errMsg = ''InPar1, [InPar2, InPar3]'';\n');
fprintf(FID,'     error(errId, errMsg);\n');
fprintf(FID,'end\n');
fprintf(FID,'Gaps = cellfun(@isempty, varargin);\n');
fprintf(FID,'DefArgs = {InPar1Def InPar2Def InPar3Def};    %% default input arguments\n');
fprintf(FID,'Suboptargs = DefArgs(1 : NumVarargs);\n');
fprintf(FID,'varargin(Gaps) = Suboptargs(Gaps);\n');
fprintf(FID,'DefArgs(1 : NumVarargs) = varargin;\n');
fprintf(FID,'[Par1 Par2 Par3] = DefArgs{:}\n');
% key,val...
fprintf(FID,'\n\n');
fprintf(FID,'DefV. = \n');
fprintf(FID,'InPar = InArg.populate_keyval(DefV,varargin,mfilename);\n');

fclose(FID);

