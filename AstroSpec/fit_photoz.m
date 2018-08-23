function [Chi2,Dof]=fit_photoz(Spec,Z,Ebv,MagSys,varargin);
%--------------------------------------------------------------------------
% fit_photoz function    
% Description: Photometric redshift \chi^2 estimator. Given a photometric
%              measurment of a source and a spectrum, calculate the \chi^2
%              fit between the spectrum as a function of redshift and the
%              photometric measurements.
% Input  : - Spectrum name (see get_spectra.m for options),
%            or spectrum matrix [Wavelength[A], flux].
%          - Redshift vector.
%          - E_{B-V} [mag] extinction, in observer frame.
%            Assumes R=3.08 (possible to change within the code).
%          - Magnitude system: {'AB' | 'Vega'}
%          - Arbitrary number of pairs of parameters:
%            ...,{FilterSystem, Filter},[Mag Err],...
%            or
%            ...,Filter,[Mag Err],...
%            in which FilterSystem is set to 'SDSS'.
% Output : - Vector of Chi2 for each redshift.
%          - Dof
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Feb 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Gal=get_spectra('Gal_Sa');
%          Gal=[Gal(:,1).*1.5, Gal(:,2)];
%          g=synphot(Gal,'SDSS','g','AB');
%          r=synphot(Gal,'SDSS','r','AB');
%          i=synphot(Gal,'SDSS','i','AB');
%          z=synphot(Gal,'SDSS','z','AB');
%          Z = [0:0.01:1]';
%          [Chi2,Dof]=fit_photoz('Gal_Sa',Z,0.016,'AB',...
%                    {'SDSS','u'},[20.96 0.195],...
%                    'g',[g 0.01],...
%                    'r',[r 0.01],...
%                    'i',[i 0.01],...
%                    'z',[z 0.01]);
%          plot(Z,Chi2);
% Reliable: 2
%--------------------------------------------------------------------------

R = 3.08;

Nz = length(Z);

DefFilterSys = 'SDSS';



Narg = length(varargin);
If = 0;
for Iarg=1:2:Narg-1,
   If = If + 1;
   if (iscell(varargin{Iarg})==0),
      FilterSys  = DefFilterSys;
      FilterName = varargin{Iarg};
   else
      FilterSys  = varargin{Iarg}{1};
      FilterName = varargin{Iarg}{2};
   end
   Filt{If} = get_filter(FilterSys,FilterName);
   ObsMag{If}  = varargin{Iarg+1};
   if (length(ObsMag{If})==1),
      ObsMag{If}(2) = NaN;
   end

   for Iz=1:1:Nz,
      Spectrum = get_spectra(Spec,[0 Ebv],[R R],Z(Iz));

      [SynMag,Dm,Ex]=basic_synthetic_photometry(Spectrum,Filt{If}.nT{1},MagSys,[],[0 1]);
      Diff(Iz,If) = SynMag - ObsMag{If}(1);
      Err(Iz,If)  = sqrt(Dm.^2 + ObsMag{If}(2).^2);
   end
end

Nfilt = If;
Dof   = Nfilt - 1;

Chi2 = sum(((Diff-mean(Diff,2)*ones(1,Nfilt) )./Err).^2,2);
