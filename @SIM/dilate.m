function SimOut=dilate(Sim,Filter,varargin)
% Dilate images in a SIM object
% Package: @SIM
% Description: Dilate images in a SIM object using imdilate.
% Input  : - A SIM object.
%          - The imdilate filter.
%          * Additional arguments to pass to ufun2sim.m
% Output : - A SIM object with the results.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: dilate(Sim,Filter)
% Reliable: 
%--------------------------------------------------------------------------

SimOut=ufun2sim(Sim,@imdilate,varargin{:},'FunAddPar',{Filter});
