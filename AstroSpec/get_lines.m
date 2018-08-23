function [Out,All]=get_lines(Wave,MaxDist);
%------------------------------------------------------------------------------
% get_lines function                                                 AstroSpec
% Description: Search spectral lines by wavelength.
%              (for search by name see: get_lines1.m)
%              The list of lines contains 46663 spectral lines
%              (Reader et al. 1980; 1981) for 99 atomic species.
%              Neutral through quadruply ionized atoms are tabulated. 
% Input  : - Wavelength.
%          - Max dist [Ang], default is 5A.
% Output : - Structure containing all the lines found within MaxDist Ang.
%            from the specified wavelength.
%            (see LinesDB.head for information).
%          - Cell array in which each cell contains information
%            about the lines that have been found.
%            Each cell contain a structure of all the other spectral
%            lines of the same element and ionization state.
%            (see LinesDB.head for information).
% Tested : Matlab 7.0
%     By : Eran O. Ofek                       May 2006
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See Also: get_lines1.m; LinesDB*
% Needed: LinesDB.mat (database of spectral lines).
% Example: [Out,All]=get_lines(6564,10);
% Reliable: 2
%------------------------------------------------------------------------------
MaxDist_Def = 5;

if (nargin==1),
   MaxDist  = MaxDist_Def;
elseif (nargin==2),
   % do nothing
else
   error('Illegal number of input arguments');
end

load LinesDB.mat;

I = find(abs(LinesDB.Wave-Wave)<=MaxDist);
size(I)

Out.Z     = LinesDB.Z(I);
Out.Ion   = LinesDB.Ion(I);
Out.Int   = LinesDB.Int(I);
Out.Wave  = LinesDB.Wave(I);
Out.Vacum = LinesDB.Vacum(I);
for J=1:1:length(I),
   Out.Name{J}  = LinesDB.Name{I(J)};

   % search for other lines with the same Ion and Z
   K = find(LinesDB.Z==Out.Z(J) & LinesDB.Ion==Out.Ion(J));

   All{J}.Z     = LinesDB.Z(K);
   All{J}.Ion   = LinesDB.Ion(K);
   All{J}.Int   = LinesDB.Int(K);
   All{J}.Wave  = LinesDB.Wave(K);
   All{J}.Vacum = LinesDB.Vacum(K);

   for L=1:1:length(K),
      All{J}.Name{L} = LinesDB.Name{K(L)};
   end
end

