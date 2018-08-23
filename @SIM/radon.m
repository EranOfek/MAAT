function SimOut=radon(Sim,varargin)
% Execute the Radon transform on SIM object images.
% Package: @SIM
% Description: Execute the Radon transform on SIM object images.
% Input  : - A SIM object.
%          * Additional arguments to pass to ufun2sim.m
% Output : - A SIM object containing the Radon transform.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S=radon(S);
%          S=radon(S,'FunAddPar',{5}); % radon transform for specific theta
% Reliable: 2
%--------------------------------------------------------------------------

SimOut=ufun2sim(Sim,@radon,varargin{:});
