function Arc=spec_get_arc(ArcName)
%--------------------------------------------------------------------------
% spec_get_arc function                                             ImSpec
% Description: Get a spectroscopic arc (template and lines list) from the
%              SpecArcs.mat database, given the arc name.
% Input  : - String containing the arc name.
%            If empty matrix, then will return the SpecArcs.mat database.
%            Default is empty matrix.
%            Available arc includes:
%            'SkyLow'  - Low resolution sky spectrum (R~300).
%            'SkyHigh' - High resolution sky spectrum (R~600).
%            'Ar'      - Argon lamp
%            'Cd'      - Cadmium lamp
%            'Hg'      - Hg lamp
%            'Ne'      - Neon lamp
%            'Zn'      - Zn lamp.
%            'Ne+Ar'   - Ne+Ar lamps.
% Output : - Structure containing the arc information, or a structure array
%            if multiple arcs are returned.
%            The following fields are available:
%            .Name   - Arc name.
%            .Spec   - Arc spectrum [Wavelength(Ang), Intensity].
%            .Lines  - List of lines, in which the first column is
%                      the line wavelength (Ang).
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Jan 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S=spec_get_arc('Ar');
%          S=spec_get_arc; {S.Name}   % will show names of all available arcs
% Reliable: 2
%--------------------------------------------------------------------------

Def.ArcName = [];
if (nargin==0),
   ArcName = Def.ArcName;
end

SpecArcs=load2('SpecArcs.mat');

if (isempty(ArcName)),
   Arc = SpecArcs;
else
   %I   = find(strcmpi({SpecArcs.Name},ArcName));
   Arc = SpecArcs(strcmpi({SpecArcs.Name},ArcName));
end


