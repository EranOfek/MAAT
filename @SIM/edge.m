function Sim=edge(Sim,Method,Thresh,varargin)
% Detect edges in images in a SIM objectusing edge.m.
% Package: @SIM
% Description: Detect edges in images in a SIM object. This is done by
%              executing edge.m.
% Input  : - A SIM object.
%          - The edge detection method. See edge.m for options.
%            Default is 'canny'.
%          - Edge detection threshold. See edge.m for details.
%          * A cell array of additional arguments to pass to ufun2sim.m.
%            Note that 'FunAddPar' will override the edge detection method
%            and threshold.
% Output : - A SIM object with the edge detected maps.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=edge(S(1),'Sobel',500);
% Reliable: 2
%--------------------------------------------------------------------------

Def.Method  = 'canny';
Def.Thresh  = {}; 
if (nargin==1),
    Method = Def.Method;
    Thresh = Def.Thresh;
elseif (nargin==2),
    Thresh = Def.Thresh;
else
    % do nothing
end

Sim = ufun2sim(Sim,@edge,'FunAddPar',{Method, Thresh},varargin{:});

