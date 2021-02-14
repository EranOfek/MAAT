function FullPath=run_ciao_command(varargin)
% RUN CIAO command on single or multiple, or all Chandra directories
% Package: +VO.Chandra
% Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Dir'   - If empty then use all directories in cats.X.ChnadraObs
%                   If an AstCat object, then assume this is has the same
%                   format as cats.X.ChnadraObs.
%                   Default is [].
%            'BasePath' - Base path in which the Chandra data is located.
%                   Default is '/euler/eran/work/Chandra'.
%            'CIAO' - CIAO program name.
%                   Default is 'ciao412'.
%            'Command' - Command to run on each data directory.
%                   Default is 'chandra_repro indir=./ outdir=./'.
% Output : - A cell array of full path on which the commad was executed.
%     By : Eran Ofek                       Aug 2020
% Example: FullPath=VO.Chandra.run_ciao_command



InPar = inputParser;
addOptional(InPar,'Dir',[]); 
addOptional(InPar,'BasePath','/euler/eran/work/Chandra'); 
%addOptional(InPar,'CIAO','unset PYTHONPATH; ciao412 -o');  % bash
addOptional(InPar,'CIAO','unsetenv PYTHONPATH; ciao412 -o');  % bash
addOptional(InPar,'Command','chandra_repro indir=./ outdir=./');
parse(InPar,varargin{:});
InPar = InPar.Results;


if isempty(InPar.Dir)
    % run over all directories
    InPar.Dir = cats.X.ChandraObs;
end

if isa(InPar.Dir,'AstCat')
    % run over ChandraObs AstCat object directories
    Obs = InPar.Dir;
    
    Nfile = size(Obs.Cat,1);
    FullPath = cell(Nfile,1);
    for Ifile=1:1:Nfile
        FullPath{Ifile} = sprintf('%s%s%s%s%s%s%d%s',InPar.BasePath,filesep,...
                                     Obs.Cat.AO{Ifile},filesep,...
                                     Obs.Cat.Cat{Ifile},filesep,...
                                     Obs.Cat.ObsID(Ifile),filesep);
    end
    
else
    FullPath = InPar.Dir;
end

if ~iscell(FullPath)
    FullPath = {FullPath};
end

PWD = pwd;

Nfile = numel(FullPath);
% run for all dir
for Ifile=2933:1:2933
    [Ifile,Nfile]
    tic;
   cd(FullPath{Ifile});
   
   CmdStr = sprintf('%s; %s',InPar.CIAO,InPar.Command);
   [Status,Result] = system(CmdStr);
    toc
end
      
cd(PWD);
