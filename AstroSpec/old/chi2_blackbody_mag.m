function [Chi2,Npar,DeltaMag]=chi2_blackbody_mag(System,Temp,Radius,Dist,Ebv,R,varargin)
%--------------------------------------------------------------------------
% chi2_blackbody_mag function                                    AstroSpec
% Description: Given a list of observed magnitudes in several bands,
%              perform a \chi^2 fit with a black-body spectrum with a
%              given temperatures and extinctions.
% Input  : - Magnitude system: {'A' | 'V'} for AB or Vega, respectively.
%          - Temperature [K].
%          - Radius [cm].
%          - Distance [pc].
%          - E_{B-V}
%          - R
%          * Arbitrary number of triplets of arguments:
%            ...,Filter,Mag,Err,Filter,Mag,Err,...
%            See blackbody_mag for filter options.
% Output : - Chi2
%          - Number of constraints (i.e., number of filters)
%          - Best fit delta magnitude.
% Tested : Matlab 7.0
%     By : Eran O. Ofek           March 2006
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: [Chi2,Npar,DeltaMag]=fit_blackbody_mag('V',14200,Rad,Dist,Ebv,3.1,...
%                               'V',5,0.01,...
%                               'R',5,0.01,...
%                               'I',5',0.01)
%--------------------------------------------------------------------------

Narg = length(varargin);
Npar = Narg./3;
if ((Npar)~=floor(Npar)),
   error('Illegal number of input arguments');
end

Ifilt = 0;
for Iarg=1:3:Narg,
   Ifilt = Ifilt + 1;
   CalcMag(Ifilt,:)     = blackbody_mag(Temp,varargin{Iarg},System,Radius,Dist).';
   A_3                  = optical_extinction(Ebv,'B','V',varargin{Iarg},'C',R);
   CalcMag(Ifilt,:)     = CalcMag(Ifilt,:) - A_3;  % intinsic magnitude
   ObsMag(Ifilt,:)      = varargin{Iarg+1}.';
   ObsMagErr(Ifilt,:)   = varargin{Iarg+2}.';
end

DeltaMag = mean(ObsMag - CalcMag);

Chi2 = sum(((ObsMag - (CalcMag + ones(Npar,1)*DeltaMag))./ObsMagErr).^2);
