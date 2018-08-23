function [Time,fftSourceF,SourceF,ShiftMatrix]=recover_source_flux(LC,TimeDelay,FluxRatio,varargin)
% Recover the source flux from the combined flux of multiple images.
% Package: AstroUtil.TimeDelay
% Description: Recover the source flux from the combined flux of multiple
%              images with known specified time delays and flux ratio.
% Input  : - The total flux light curve of all time-delayed lensed images.
%            This is a three column matrix of [Time, Flux, Error].
%          - A vector of time delays.
%          - A vector of flux ratio of the images relative to the first
%            image.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Resample'
%            'Taper'
%            'ColT'
%            'ColF'
%            'ColE'
% Output : - The vector of times (resampled).
%          - The fft of the recovered source flux.
%          - Thr recovered source flux.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Aug 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: T=timeseries.random_time_sequence;
%          [T,F,TotalF,X,Y]=AstroUtil.TimeDelay.generate_timedelay_lc(T);
%          LC=[T,TotalF, 0.05.*ones(size(T))];
%          LC=LC(~isnan(TotalF),:);
%          [TT,~,SourceF]=AstroUtil.TimeDelay.recover_source_flux(LC,[0 50],[1 0.8]);
%          plot(T,F(:,1)); hold on; plot(TT,SourceF);
% Reliable: 
%--------------------------------------------------------------------------


DefV.Resample             = false;
DefV.ResamplePar          = {};
DefV.HighestFreq          = 1;  % 1/day
DefV.InterpMethod         = 'linear';
DefV.InterpMethodErr      = 'nearest';
DefV.Taper                = [];
DefV.ColT                 = 1;
DefV.ColF                 = 2;
DefV.ColE                 = 3;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% make TimeDelay and FluxRatio to column vectors
%TimeDelay = TimeDelay(2:end);
%FluxRatio = FluxRatio(2:end);
TimeDelay = TimeDelay(:).';
FluxRatio = FluxRatio(:).';

LC   = sortrows(LC,InPar.ColT);

Time  = LC(:,InPar.ColT);
Freq  = 1./Time;
TotF  = LC(:,InPar.ColF);
ErrF  = LC(:,InPar.ColE);
Nobs  = numel(Time);
TimeRange = range(Time);
MeanF = 0; %mean(TotF);
TotF  = TotF - MeanF;


if (~isempty(InPar.Taper))
    Frac      = (TimeRange - max(TimeDelay))./TimeRange;
    [CosB]    = timeseries.cosbell(Frac,Time);
    TotF      = TotF.*CosB;
end


if (~InPar.Resample)
    % evenly spaced series
    fftTotF = fft(TotF);
    fftErrF = fft(ErrF);
else
    % time series is not evenly spaced
%     if (isempty(InPar.Resample))
    % do not resample light curve - use unevely spaced Fourier
    % Transform

%     LowFreq = 0.5./TimeRange;
%     Freq    = (LowFreq:LowFreq:InPar.HighestFreq)';
% 
%     fftTotF = timeseries.period_complex([Time,TotF],Freq,'amp');
%     fftErrF = timeseries.period_complex([Time,ErrF],Freq,'amp');
    %Time    = 1./Freq;
%         
%     else
        % resample light curve to a uniform grid
       
    [IntVal,Time]=resample_uniform(Time,[TotF,ErrF],InPar.ResamplePar{:});
    TotF = IntVal(:,1);
    ErrF = IntVal(:,2);
    
    fftTotF = fft(TotF);
    fftErrF = fft(ErrF);
    %end
    
end
DeltaT = Time(2) - Time(1);

ShiftMatrix = bsxfun(@times,FluxRatio,exp(2.*pi.*1i.*bsxfun(@times,Freq,TimeDelay)));
ShiftFactor = sum(ShiftMatrix,2);

%ShiftFactor = (1 + sum(FluxRatio.*exp(2.*pi.*1i.*Freq.*TimeDelay)));

fftSourceF = (fftTotF - fftErrF)./ShiftFactor;

if (nargout>2)
    SourceF = (real(ifft(fftSourceF)) + MeanF); %./sum(FluxRatio);
end


    
    


