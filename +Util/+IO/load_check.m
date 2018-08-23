function Var=load_check(File,StoreWS,WS)
% Load, but check if variable exist in workspace.
% Package: Util.IO
% Description: Load a matlab variable or file from disk (similar to the
%              load.m command). However, before the variable is loaded the
%              function checks if the variable with name identical to the
%              file name is already exist in the matlab main workspace.
%              If it is exist it will copy
%              the variable from the workspace.
%              If variable does not exist it will load it in the usual way
%              and store it in the main workspace.
%              This is usefull when you want to load big variables
%              in multiple function calles.
% Input  : - String of file name to load
%          - Store the variable in the main work space, if doesn't exist,
%            {true|false}. Default is true.
%          - Workspace {'base'|'caller'}, Default is 'base'.
% Output : - Loaded variable.
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    May 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Var=Util.IO.load_check('Tycho2.mat');
% Reliable: 1
%--------------------------------------------------------------------------

Def.StoreWS = true;
Def.WS = 'base';
if (nargin==1)
   StoreWS = Def.StoreWS;
   WS      = Def.WS;
elseif (nargin==2)
   WS      = Def.WS;
elseif (nargin==3)
   % do nothing
end

RE = regexp(File,filesep,'split');
FileNM = RE{end}; 
Ind = strfind(FileNM,'.mat');
FileNM = FileNM(1:Ind-1);

if (evalin(WS,sprintf('exist(''%s'')',FileNM))==1)
   % variable exist in workspace
   Var = evalin(WS,sprintf('%s;',FileNM));
else
   Var = Util.IO.load2(File);
   %eval(sprintf('Var=%s;',FileNM));
   if (StoreWS)
        assignin(WS,FileNM,Var);
   end
end
