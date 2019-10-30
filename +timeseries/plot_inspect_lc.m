function [Prop,PS]=plot_inspect_lc(LC,varargin)
% plot and inspect light curve and power spectrum
% Package: timeseries
% Description: 
% Input  : - Light curve [time, mag, error]
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - Structure of light curve properties.
%          - Power spectrum.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Prop,PS]=plot_inspect_lc(LC)
% Reliable: 
%--------------------------------------------------------------------------


DefV.PS                   = [];
DefV.MaxFreq              = 200;
DefV.TypePS               = 'normnl';
DefV.PhaseBinSize         = 0.1;
DefV.BootStrapQuantil     = 0.99;
DefV.ColT                 = 1;
DefV.ColM                 = 2;
DefV.ColE                 = 3;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

LC = LC(:,[InPar.ColT, InPar.ColM, InPar.ColE]);
% calculate LC properties
Prop.Std       = std(LC(:,2));
Prop.RStd      = Util.stat.rstd(LC(:,2));
Prop.TimeRange  = range(LC(:,1));
Prop.TimeStd    = std(LC(:,1));

FreqVec = (0:1./(2.*Prop.TimeRange):InPar.MaxFreq).';
[PS,Peaks] = timeseries.period(LC,FreqVec,'Type',InPar.TypePS);
Ipeak = numel(Peaks.Freq);

%PSb = timeseries.period_norm_bootstrap(LC,FreqVec);
PSw = timeseries.period(LC,FreqVec,'Type','WinNL');

% LC
figure(1);
cla;
plot.errorxy(LC);

% power spectrum
figure(2);
cla;
plot(PS(:,1),PS(:,2),'k-');
hold on;
% mark highest peak
Hp = plot(Peaks.Freq(Ipeak),Peaks.PS(Ipeak),'ro');
plot(PSw(:,1),PSw(:,2),'r-')
%plot(FreqVec,quantile(PSb,InPar.BootStrapQuantil,2),'-','Color',[0.8 0.8 0.8]);


Cont = true;
while Cont
    % folded LC
    FLC = timeseries.folding(LC,Peaks.Per(Ipeak));
    fprintf('Folded period: %f\n',Peaks.Per(Ipeak));
    
    figure(3);
    cla;
    B   = timeseries.binning(FLC,InPar.PhaseBinSize,[0 1],{'MidBin',@median,@std,@numel});
    plot(FLC(:,1),FLC(:,2),'o','Color',[0.8 0.8 0.8]);
    hold on;
    B(:,3) = B(:,3)./sqrt(B(:,4));
    plot.errorxy(B)


    % select another period
    R = input('select another period (n) or continue [Y/n]','s');
    switch lower(R)
        case 'n'
            fprintf('Click on power spectrum to select another peak');
            [x,y]=ginput(1);
            [~,Ipeak] = min(abs(x-Peaks.Freq));
            
            figure(2);
            delete(Hp);
            Hp = plot(Peaks.Freq(Ipeak),Peaks.PS(Ipeak),'ro');
            
            Cont = true;
        otherwise
            % do nothing
            Cont = false;
    end
end





