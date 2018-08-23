function SimOut=cast(Sim,Type,varargin)
% Convert all the image elements in a SIM array into a specific class
% Package: @SIM
% Description: Convert all the image elements in a SIM array into
%              a specific class (e.g., 'int64') using the cast command.
%              By default this operates on the image field. In order to
%              run this on additional fields see ufun2sim.m
% Input  : - A SIM object.
%          - A class to into which to convert the SIM object elements.
%          * Additional arguments to pass to ufun2sim.m
% Output : - A SIM object array with the results.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: cast(S,'logical')
% Reliable: 2
%--------------------------------------------------------------------------

SimOut=ufun2sim(Sim,@cast,varargin{:},'FunAddPar',{Type});
