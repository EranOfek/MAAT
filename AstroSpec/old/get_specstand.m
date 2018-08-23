function [Spec,Sp,Ind]=get_specstand(Name,Units);
%-------------------------------------------------------------------------
% get_specstand function                                        AstroSpec
% Description: Get the spectrum of a spectrophotometric standard.
% Input  : - Standard name. If empty matrix (i.e., []), then list
%            all the available spectrophotometric standards and their
%            coordinates.
%          - Units of output spectrum:
%            'fl'   - Flux in erg/cm^2/s/A (default).
% Output : - Spectrum [Wavelength(A), Flux] of the requested
%            spectrophotometric standard.
%            Return empty matrix of not found.
%          - Structure of all the availablde spectrophotometric standards.
%          - Index of selected standard in the structure of standards.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                           April 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-------------------------------------------------------------------------
ColW   = 1;
ColF    = 2;

DefUnits   = 'fl';
if (nargin==1),
   Units    = DefUnits;
elseif (nargin==2),
   % do nothing
else
   error('Illegal number of input arguments');
end


%---------------------------------------------------
I = 1;
Sp.Name{I}  = 'BD+25 4655';
Sp.File{I}  = 'fbd25d4655.dat';
Sp.RA(I)    = convertdms([21 59 41.9749],'H','r');
Sp.Dec(I)   = convertdms([1 26 25 57.388],'D','R');
Sp.V(I)     = 9.68;
Sp.Type{I}  = 'WD';
Sp.Units(I) = 1e16;   % divide by this to get flux in erg/cm^2/s/A

I = 2;
Sp.Name{I}  = 'HZ44';
Sp.File{I}  = 'fhz44.dat';
Sp.RA(I)    = convertdms([13 23 35.2581],'H','r');
Sp.Dec(I)   = convertdms([1 36 07 59.514],'D','R');
Sp.V(I)     = 11.673;
Sp.Tyep{I}  = 'B2';
Sp.Units(I) = 1e16;


I = 3;
Sp.Name{I}  = 'BD+17 4708';
Sp.File{I}  = 'bd17_4708.dat';
Sp.RA(I)    = convertdms([22 11 31.3737],'H','r');
Sp.Dec(I)   = convertdms([1 18 05 34.174],'D','R');
Sp.V(I)     = 9.464;
Sp.Tyep{I}  = 'sdF8';
Sp.Units(I) = 1;

I = 4;
Sp.Name{I}  = 'HZ2';
Sp.File{I}  = 'fhz2.dat';
Sp.RA(I)    = convertdms([04 12 43.55],'H','r');
Sp.Dec(I)   = convertdms([1 11 51 49.0],'D','R');
Sp.V(I)     = 13.877;
Sp.Tyep{I}  = 'DA:';
Sp.Units(I) = 1e16;


I = 5;
Sp.Name{I}  = 'HZ21';
Sp.File{I}  = 'fhz21.dat';
Sp.RA(I)    = convertdms([12 13 56.27],'H','r');
Sp.Dec(I)   = convertdms([1 32 56 31.4],'D','R');
Sp.V(I)     = 14/688;
Sp.Tyep{I}  = 'DA';
Sp.Units(I) = 1e16;

%---------------------------------------------------

Flag = isempty_cell(strfind(Sp.Name,Name));

Ind  = find(Flag==0);

if (isempty(Ind)==1),
   % not found
   Spec = [];
else
   Spec = load(Sp.File{Ind});
   Spec = [Spec(:,ColW), Spec(:,ColF)./Sp.Units(Ind)];
end

