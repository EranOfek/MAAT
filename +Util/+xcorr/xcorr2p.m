function [CC,SimCorr]=xcorr2p(List1,List2,AnnuVec,Geom,CorrType,Nsim,RandFun,varargin);
%------------------------------------------------------------------------------
% xcorr2p function                                                   AstroStat
% Description: Given two lists of coordinates (selected in the same region),
%              calculate the cross correlation function as a function of
%              distance. The cross-correlation can be calculated using the
%              DR, RR or DR2 methods.
% Input  : - First list [X,Y] or [RA,Dec] in radians.
%          - Second list [X,Y] or [RA,Dec] in radians.
%          - Correlation annuli [Start:BinSize:End] in radians.
%          - Distance method:
%            'sphere' : using spherical coordinates, default.
%            'plane'  : Using plane coordinates.
%          - Type of correlation:
%            'DR'     : 2*DD/DR - 1      (default)
%            'RR'     : (DD-DR+RR)/RR
%            'DR2'    : 4*DD*RR/(DR*DR) - 1
%          - Number of bootstrap simulations for error estimation,
%            default is 0.
%          - Name of function [Coo]=fun(N), that return N
%            random coordinates in the area in which the
%            coordinates in List1/2 are taken from.
%            Default is 'cel_coo_rnd'.
%          * Arbitrary number of input arguments that
%            passed to the function.
% Output : - Cross correlation [BinCenter, Correlation, Error],
%            where error is the analytic error.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                  November 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------

RAD = 180./pi;

Stat          = 'corr';

DistMethodDef = 'sphere';
CorrTypeDef   = 'DR';
NsimDef       = 0;
RandFunDef    = 'cel_coo_rnd';


if (nargin==3),
   Geom         = DistMethodDef;
   CorrType     = CorrTypeDef;
   Nsim         = NsimDef;
   RandFun      = RandFunDef;
elseif (nargin==4),
   CorrType     = CorrTypeDef;
   Nsim         = NsimDef;
   RandFun      = RandFunDef;
elseif (nargin==5),
   Nsim         = NsimDef;
   RandFun      = RandFunDef;
elseif (nargin==6),
   RandFun      = RandFunDef;
else (nargin==7),
   % do nothing
end

N1 = size(List1,1);
N2 = size(List2,1);

SimCorr = zeros(Nsim,length(AnnuVec)-1);

for Isim=0:1:Nsim,

   if (Isim==0),
      % Real sample (no bootstrap)
      SimList1 = List1;
      SimList2 = List2;
   else
      % Bootstrap
      % select random from lists
      SimInd1 = unidrnd([1:1:N1]');
      SimInd2 = unidrnd([1:1:N2]');
      SimList1 = List1(SimInd1,:);
      SimList2 = List2(SimInd2,:);
   end


   if (isempty(varargin)==1),
      RandList2 = feval(RandFun,N2);
   else
      RandList2 = feval(RandFun,N2,varargin{:});
   end
   
   
   [X,N_DD] = hist_dist2p(SimList1,SimList2    ,Geom,AnnuVec,Stat);
   [X,N_DR] = hist_dist2p(SimList1,RandList2,Geom,AnnuVec,Stat);
   DD     = N_DD./sum(N_DD);            % normalize...
   DR     = N_DR./sum(N_DR);            % normalize...
   ErrDD  = sqrt(N_DD)./sum(N_DD);
   ErrDR  = sqrt(N_DR)./sum(N_DR);
   
   
   switch CorrType
    case 'DR'
       % do nothing
       % 'DR'     : 2*DD/DR - 1
       Corr    = 2.*DD./DR - 1;
       ErrCorr = sqrt((2.*ErrDD./DR).^2 + (2.*DD.*ErrDR./(DR.^2)).^2);
   
    otherwise
       % RR correlation
       RandList1 = eval(sprintf('%s(%d)',RandFun,N1));
       RandList2 = eval(sprintf('%s(%d)',RandFun,N2));
       [X,N_RR,C_RR] = hist_dist2p(RandList1,RandList2,Geom,AnnuVec,Stat);
       RR     = N_RR./sum(N_RR);
       ErrRR  = sqrt(N_RR)./sum(N_RR);
   
       switch CorrType
          case 'RR' 
          % 'RR'     : (DD-DR+RR)/RR
          Corr    = (DD - DR + RR)./RR;
          ErrCorr = sqrt((ErrDD./RR).^2 + (ErrDR./RR).^2 + ((DD-DR).*ErrRR./(RR.^2)).^2);
   
          case 'DR2'
          % 'DR2'    : 4*DD*RR/(DR*DR) - 1
          Corr    = 4.*DD.*RR./(DR.*DR) - 1;
          ErrCorr = sqrt((4.*ErrDD.*RR./(DR.*DR)).^2 + (4.*DD.*ErrRR./(DR.*DR)).^2 + (2.*4.*DD.*RR.*ErrDR./(DR.^3)).^2);
   
   
          otherwise
             error('Unknown CorrType Option');
       end
   end

   if (Isim==0),
      % save non-bootstrap result
      CC = [X, Corr, ErrCorr];
   else
      % bootstrap result
      % each raw of SimCorr contains the coorelations for a single realization
      SimCorr(Isim,:) = Corr.';   

   end

end



