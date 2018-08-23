function [Chi2_BBS,Dof,SigmaCounter]=period_events(Events,BackgroundEvents,AreaRatio,Nbin,Period,MeanMethod,Nsigma);
%-----------------------------------------------------------------------------
% period_events function                                           timeseries
% Description: Search for periodicity in a list of "time tagged"
%              events. The search is done by folding the events
%              and background event to each trial period. In each
%              trial period the data is binned and the \chi^2
%              for a constant rate is calculated.
% Input  : - Vector of events [time tags].
%          - Vector of background events [time tags].
%          - Ratio between the area of aperture from which the background
%            events were extracted and the area of the aperture from which
%            the events were extracted (=BackgroundArea/SourceArea).
%          - Number of bins in each trial period.
%          - Vector of trial periods to test.
%          - Method in which to calculate mean count rate:
%            {'wmean'|'mean'|'median}, default is 'median'.
%          - Number of sigmas from mean which will be counted
%            (see ouput: Sigma counter), default is 2.
% Output : - Chi^2 after subtracting the actual background rate in each bin.
%          - Number of degrees of freedom.
%          - Sigma Counter: Number of bins which falls [below, above] Nsigma
%            from mean.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                       Feb 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%-----------------------------------------------------------------------------
DefMeanMethod   = 'median';
DefNsigma       = 2;

if (nargin==5),
   MeanMethod   = DefMeanMethod;
   Nsigma       = DefNsigma;
elseif (nargin==6),
   Nsigma       = DefNsigma;
elseif (nargin==7),
   % do nothing
else
   error('Illegal number of input arguments');
end

Epoch0 = 0;

Np  = length(Period);
%Npd = length(PeriodDeriv);

Nev     = length(Events);
NevBack = length(BackgroundEvents);


Dof = Nbin - 1;

Chi2_MBS     = zeros(Np,1);
Chi2_BBS     = zeros(Np,1);
SigmaCounter = zeros(Np,2);
for Ip=1:1:Np,

   PerPhase     = (Events - Epoch0)./Period(Ip);
   Phase        = PerPhase - floor(PerPhase);

   PerPhaseBack = (BackgroundEvents - Epoch0)./Period(Ip);
   PhaseBack    = PerPhaseBack - floor(PerPhaseBack);

   [Xbin,Nbin_ev]     = realhist(Phase,    [0 1 Nbin]);
   [Xbin,Nbin_evBack] = realhist(PhaseBack,[0 1 Nbin]);
%Nbin_evBack = 0;   % no background subtraction
   % [N, Err] in bin
   %-----------------
   % number of event in each bin corrected for the Mean background rate
   %Nbin_evMBS         = [Nbin_ev - Nev./(Nbin.*AreaRatio),...
   %			 sqrt(Nbin_ev + Nev./(Nbin.*AreaRatio.^2))];

   % number of event in each bin corrected for the background rate in each bin
   Nbin_evBBS         = [Nbin_ev - Nbin_evBack./AreaRatio,...
     		         sqrt(Nbin_ev + Nbin_evBack./(AreaRatio.^2))];

   
   switch lower(MeanMethod)
    case 'wmean'
       [MeanBBS,ErrMeanBBS] = wmean(Nbin_evBBS);
    case 'mean'
       MeanBBS              = mean(Nbin_evBBS(:,1));
    case 'median'
       MeanBBS              = median(Nbin_evBBS(:,1));
    otherwise
       error('Unknown MeanMethod option');
   end

   %[MeanMBS,ErrMeanMBS] = wmean(Nbin_evMBS);
   %Chi2_MBS(Ip) = sum( (Nbin_evMBS(:,1) - MeanMBS).^2./(Nbin_evMBS(:,2).^2 + ErrMeanMBS.^2) );
   % non bias corrected
   %Chi2_BBS(Ip) = sum( (Nbin_evBBS(:,1) - MeanBBS).^2./(Nbin_evBBS(:,2).^2 + ErrMeanBBS.^2) );

   % bias corrected
   Chi2_BBS(Ip) = sum( (Nbin_evBBS(:,1) - MeanBBS).^2./abs(MeanBBS) );


   I_sigUp  = find((Nbin_evBBS - MeanBBS)>( Nsigma.*sqrt(MeanBBS)));
   I_sigLow = find((Nbin_evBBS - MeanBBS)<(-Nsigma.*sqrt(MeanBBS)));

   SigmaCounter(Ip,:) = [length(I_sigLow), length(I_sigUp)];

end
