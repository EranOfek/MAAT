function [Chi2,DeltaM,Npar]=fit_blackbody_mag(VecT,VecEbv,VecR,Radius,Dist,System,varargin);
%--------------------------------------------------------------------------
% fit_blackbody_mag function
% Description: 
% Input  : - Vector of temperature [K]
%          - Vector of E_{B-V} [mag]
%          - Vector of R []
%          - blackbody radii [cm] for scaling
%          - blackbody distances [pc] for scaling
%          - Magnitude system: {'A' | 'V' | 'F'} for AB, Vega or flux, respectively.
%          * Arbitrary number of triplets of arguments:
%            ...,Filter,Mag,Err,Filter,Mag,Err,...
%            See blackbody_mag for filter options.
%            If flux then Filter is in wavelength in ang.
% Output : - Cube of \chi^2
%          - Cube of best fit scaling (in magnitude) for each parameter 
%          - Number of constraints (i.e., number of filters)
% Tested : Matlab 7.0
%     By : Eran O. Ofek           March 2006
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------------
Pc = get_constant('pc','cgs');

NT    = length(VecT);
NEbv  = length(VecEbv);
NR    = length(VecR);

Chi2   = zeros(NT,NEbv,NR);
DeltaM = zeros(NT,NEbv,NR);
for IT=1:1:NT,
   for IEbv=1:1:NEbv,
      for IR=1:1:NR,
         switch System
          case 'F'
   	     %error('Flux option not working')
             % flux
	     J = 0;
 	     for I=1:3:(length(varargin)-2),
	        J = J + 1;
	        FluxBB(J)      = black_body(VecT(IT),varargin{I}).*1e-8.*Radius.^2./((Dist.*Pc).^2);
                Extin          = optical_extinction(VecEbv(IEbv),'B','V',varargin{I}.*1e-4,'C',VecR(IR));
                FluxBB(J)      = FluxBB(J)./(2.512.^Extin);
                FluxObs(J)     = varargin{I+1};
                FluxErrObs(J)  = varargin{I+2};
             end
             [DM, ErrDM] = wmean([FluxBB.'./FluxObs.',FluxErrObs.']);

             C2 = sum((FluxBB.'./DM - FluxObs.').^2./(FluxErrObs.'.^2 + ErrDM.^2));
             DM = -2.5.*log10(DM);
             Npar = length(FluxBB);
          otherwise
             [C2,Npar,DM] = chi2_blackbody_mag(System,VecT(IT),Radius,Dist,VecEbv(IEbv),VecR(IR),varargin{:});
         end
         Chi2(IT,IEbv,IR)   = C2;
         DeltaM(IT,IEbv,IR) = DM;
      end
   end
end
