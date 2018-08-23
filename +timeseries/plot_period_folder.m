function VecPeriod=plot_period_folder(PowerSpec,LC,FoldBinSize);
%------------------------------------------------------------------------------
% plot_period_folder function                                       timeseries
% Description: Given the power spectrum and light curve, let the user to
%              interactively select peaks in the power spectrum, and
%              present the folded (and binned) light curve for the
%              selected periods.
% Input  : - Power spectrum [frequency, power].
%          - Light curve [Time, Magnitude, Error], where Error is
%            optional.
%            If only one column is given [Time], then assume
%            the light curve is photon time tags, and present
%            an histogram of counts as a function o phase.
%          - Bin size for folding [phase units], default is 0.2.
% Output : - Vector of all selected periods.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                       Feb 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------

RangePS = 10;   % Range of power spectrum window in which to search for local maximum

if (nargin==2),
   FoldBinSize = 0.2;
elseif (nargin==3),
   % do nothing
else
   error('Illegal number of input arguments');
end

if (size(LC,2)==1),
   Type = 'Count';
else
   Type = 'Mag';
end

Hfig = figure;
plot(PowerSpec(:,1),PowerSpec(:,2));   % plot power spectrum
set(gca,'FontSize',14);
H=xlabel('Frequency'); set(H,'FontSize',18);
H=ylabel('Power'); set(H,'FontSize',18);

HfigFold = figure;

%[Peaks,PeaksInd,WMeanPos]=fmaxs(PowerSpec,1,2);

Ip = 0;
RightClick = 0;
while (RightClick==0),
   figure(Hfig);
   zoom on;
   R1 = input('Zoom in is on - press any key to continue','s');
   zoom off;
   disp('Select a peak in the power spectrum [select by left click, abort by right click]')
   figure(Hfig);
   [Freq,Pow,But] = ginput(1);

   if (But==1),
      % search for nearest local maximum in power spectrum
      [Min,Ind] = min(abs(Freq-PowerSpec(:,1)));
      [Peaks,PeaksInd,WMeanPos]=fmaxs(PowerSpec(Ind-RangePS:Ind+RangePS,:),1,2);
      [MinDist,MinInd] = min(abs(Freq-Peaks(:,1)));
   
      Period = 1./Peaks(MinInd,1);
      Ip = Ip + 1;
      VecPeriod(Ip) = Period;

      disp(sprintf('Selected maximum - Frequency: %f, Power: %f, Period: %f',Peaks(MinInd,1),Peaks(MinInd,2),Period))

      switch Type
       case 'Mag'      
          % LC is [Time, Mag, Err]
          FoldedLC = folding(LC,Period,1);
          BinFoldLC=binning(FoldedLC,FoldBinSize,0,1);
      
          figure(HfigFold),
          clf(HfigFold);
          plot(FoldedLC(:,1),FoldedLC(:,2),'r.');
          hold on;
          if (size(LC,2)==2),
	     plot(BinFoldLC(:,1),BinFoldLC(:,2),'o');
          else
             errorxy(BinFoldLC,'o');
          end
          set(gca,'FontSize',14);
          H=xlabel('Phase'); set(H,'FontSize',18);
          H=ylabel('Intensity'); set(H,'FontSize',18);
          invy;
          hold off;
       case 'Count'
          % LC is photon time tags
          Phase = LC./Period;
          Phase = Phase - floor(Phase);
          [P,PC]=realhist(Phase,[0 1 1./FoldBinSize]);

          figure(HfigFold),
          clf(HfigFold);
          bar(P,PC);
          set(gca,'FontSize',14);
          H=xlabel('Phase'); set(H,'FontSize',18);
          H=ylabel('Counts'); set(H,'FontSize',18);

       otherwise
          error('Unknown Type Option');
      end
   else
      % abort
      RightClick = 1;
   end
end
