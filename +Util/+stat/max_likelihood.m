function [ML,Chi2,MLRT,SimML,P_SimML]=max_likelihood(ParNumDist, Events, DoF, Nsim, Psig, NFactor,AddOffset,InterpMethod)
% Likelihood from observations and numerical probability distribution.
% Package: Util.stat
% Description: Given a numerical probability distribution and list of
%              'events', calculate the 'likelihood' for the events given
%              the probability distribution. In addition, calculate the
%              ML ratio test and preform Monte-Carlo simulations by
%              generating realizations of events given the probability
%              distribution and calculate the likelihood-probability
%              distribution.
% Input  : - Two column matrix of parent numerical distribution.
%            [X, P], where P is the probabiliy density.
%            The program will normalize P if necessary.
%            If one column matrix is given, then the numerical probability
%            distribution is calculated by an histogram.
%          - List of event [X].
%          - Number of Degree's of Freedom. Default is 1
%            (Free parameters of the model).
%          - Number of Monte-Carlo simulations. Default is 1000.
%          - Column vector of Probabilities for the ML ratio-test.
%            Return the MLRT for each probability.
%            Default is [0.6827; 0.9545; 0.9973].
%            Usefull when ML is calculated as function of a free parameters.
%          - Optional value: If (ParNumDist) is single-column, then
%            generate ParNumDist so that each bin will have (Nfactor) points.
%            Default is 100.
%          - Optional offset to add to the objects and to the
%            Monte-Carlo objects.
%            In each MC simulation, add an offset to the simulated events
%            before calculating the ML.
%            Default is 0 - no offset.
%          - Interpolation method in calculating the probability from
%            the numerical function. Default is 'linear'.
% Output : - log Likelihood.
%          - Chi2 per DoF.
%          - ML Ratio-test.
%          - Vector of simulated MLs.
%          - Given the ML and Simulated ML, return the probability to
%            get the Events ML from the Parent Simulated ML distribution:
%            P(SimML>ML).
% Tested : Matlab 5.3
%     By : Eran O. Ofek          September 2001
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: % calculate likelihood for events from a gaussian distribution
%          x=[-10:0.01:10]';
%          y=exp(-x.^2./(2))./sqrt(2.*pi);
%          [ML,Chi2,MLRT,SimML,P_SimML]=max_likelihood([x,y],randn(100,1));
% Reliable: 2
%------------------------------------------------------------------------------
NsimDef    = 1000;   % default number of simulations
PsigDef    = [0.6827; 0.9545; 0.9973];
NFactorDef = 100;
InterpMethodDef = 'linear';
ProbType   = 'd';
NewSeed    = 'n';

if (nargin==2)
   DoF     = 1;
   Nsim    = NsimDef;
   Psig    = PsigDef;
   NFactor = NFactorDef;
   AddOffset = 0;
   InterpMethod = InterpMethodDef;
elseif (nargin==3)
   Nsim    = NsimDef;
   Psig    = PsigDef;
   NFactor = NFactorDef;
   AddOffset = 0;
   InterpMethod = InterpMethodDef;
elseif (nargin==4)
   Psig    = PsigDef;
   NFactor = NFactorDef;
   AddOffset = 0;
   InterpMethod = InterpMethodDef;
elseif (nargin==5)
   NFactor = NFactorDef;
   AddOffset = 0;
   InterpMethod = InterpMethodDef;
elseif (nargin==6)
   AddOffset = 0;
   InterpMethod = InterpMethodDef;
elseif (nargin==7)
   InterpMethod = InterpMethodDef;
elseif (nargin==8)
   % do nothing
else
   error('Illegal number of input arguments');
end


SizeNp = size(ParNumDist);
% number of points in numerical distribution.
Np     = SizeNp(1);
if (SizeNp(2)==1)
   %--- Single column probability ---
   Nh     = floor(Np./NFactor);
   %[min(ParNumDist), max(ParNumDist), Nh]
   [X, P] = realhist(ParNumDist,[min(ParNumDist), max(ParNumDist), Nh]);
   ParNumDist = [X, P];
   MinProb    = 0.5./Np;
elseif (SizeNp(2)==2)
   %--- Two column probability ---
   MinProb    = 0.5.*min(ParNumDist(:,2));
else
   error('parent numerical distribution should be two column matrix');
end

Ne      = length(Events);

% sorted events vector
SEvents = sort(Events);

% normalize the parent numerical distribution
ParNumDist(:,2) = ParNumDist(:,2)./trapz(ParNumDist(:,1), ParNumDist(:,2));

% calculate the probability density per event
SEvents = SEvents + AddOffset;
PE      = interp1(ParNumDist(:,1), ParNumDist(:,2), SEvents, InterpMethod);

%
%PE = (1./1.02).*exp(-1.02.*SEvents);


K       = find(isnan(PE)==1);
PE(K)   = MinProb;
K       = find(PE==0);
PE(K)   = MinProb;
% calculate the ML
ML      = sum(log(PE));


% ML ratio test
Chi2 = chi2inv(Psig,DoF);
MLRT = -0.5.*Chi2;


% Monte-Carlo simulation
%rand('state',sum(100*clock));
if (nargout>3)
   SimML = zeros(Nsim,1);
   for I=1:1:Nsim

      %SimEvents  = randgen([ParNumDist(:,1), ParNumDist(:,2)./max(ParNumDist(:,2))],Ne,ProbType,InterpMethod,NewSeed);
      SimEvents  = randgen([ParNumDist(:,1), ParNumDist(:,2)],Ne,ProbType,InterpMethod,NewSeed);

%
%SimEvents = exprnd((1./1.02),Ne,1);

      SimSEvents = sort(SimEvents);

      % optional offset
      SimSEvents = SimSEvents + AddOffset;

      % calculate the probability density per simulated event
      PSE      = interp1(ParNumDist(:,1), ParNumDist(:,2), SimSEvents,'linear');
      %
      %PSE = (1./1.02).*exp(-1.02.*SimSEvents);


      K        = find(isnan(PSE)==1);
      PSE(K)   = MinProb;
      K        = find(PSE==0);
      PSE(K)   = MinProb;
      SimML(I) = sum(log(PSE));

   end

   % The probability to get the ML from the Monte-Carlo.
   P_SimML = length(find(SimML>ML))./Nsim;

end


