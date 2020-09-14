function [Flux,Res]=calib_flux_sum(Flux,varargin)
% Normalize (calibrate) a matrix of fluxes using sum of flux in each epoch.
% Package: timeseries
% Description: Normalize (calibrate) a matrix of fluxes (Stars, Epochs)
%              using sum of flux in each epoch.
%              Select stars that appear in all epochs.
%              Calculate the sum of flux of these stars per epoch.
%              Normalize by the mean sum of fluxes.
%              Divide the stars fluxes by the normaliztion.
%              Remove stars based on quantile and max rms and iterate.
% Input  : - Matrix of flux measurments (Flux_ij), for star i at epoch j.
%            THe values should be in units of electrons (see 'Gain').
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Gain' - Gain by which to multiply to convert to electrons.
%                     Default is 1.
%            'MeanNormFun' - Function to use for the sum flux normaliztion.
%                     Default is @mean.
%            'MaxRMS' - Maximum relative rms of calibration stars.
%                     Default is 0.1.
%            'QuantileRMS' - Fraction of best rms stars to use in
%                     calibration. Default is 0.5.
%            'Niter' - Number of iterations. Default is 3.
% Output : - Matrix of calibrated flux.
%            A structure containing the following fields:
%            'MeanFlux' - Vector of mean flux per star.
%            'StarsRelRMS' - Vector of stars relative rms after
%                       calibration.
%            'FlagCS' - Flag indicating stars used for calibration.
%            'rmsZP'  - rms of ZP (normalization) at each epoch, as
%                       100measured from the mean-normalized calibrated stars.
%            'rmsflux_selectPar' - Cell array of key,val parameters to pass
%                       to timeseries.rmsflux_select.m.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: MeanFlux = rand(100,1).*1e5;
%          Norm     = 1+0.05.*randn(1,200); 
%          Flux     = poissrnd(MeanFlux.*Norm);
%          [CalFlux,Res] = timeseries.calib_flux_sum(Flux);
%          loglog(Res.MeanFlux,Res.StarsRelRMS,'.')
% Reliable: 2
%--------------------------------------------------------------------------



DefV.Gain                = 1;
DefV.MeanNormFun         = @mean;
DefV.MaxRMS              = 0.1;
DefV.QuantileRMS         = 0.5;
DefV.Niter               = 3;
DefV.StarsInBin          = 30;
DefV.rmsflux_selectPar   = {};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

[Nst,Nep] = size(Flux);

Flux = Flux.*InPar.Gain;

% select stars that have measurments in all epochs
IsN   = isnan(Flux);
NepNN = sum(~IsN,2);

% selected Calibration Stars
FlagNN = NepNN==Nep;
FlagCS = FlagNN;
if (sum(FlagCS)==0)
    error('No stars that appear in all epochs');
end

% Flux matrix of calibrated stars
FluxCS = Flux(FlagCS,:);

for Iter=1:1:InPar.Niter
    % sum flux of calibration stars in each epoch
    SumFluxEpoch = sum(FluxCS,1);

    Norm = SumFluxEpoch./InPar.MeanNormFun(SumFluxEpoch);

    % Calibrated Flux matrix
    Flux = Flux./Norm;

    % Mean flux of each star over all epochs
    MeanFlux = mean(Flux,2);

    StarsRelRMS = std(Flux./MeanFlux,[],2);

    Q = quantile(StarsRelRMS,InPar.QuantileRMS);
    
    
    

    % Select stars by removing outliers from the RMS vs. Flux diagram.
    %[SortedMeanFlux,SI] = sort(MeanFlux);
    %SortedRMS           = StarsRelRMS(SI);
    
    [FlagRMS,ResS]=timeseries.rmsflux_select(MeanFlux,StarsRelRMS,InPar.rmsflux_selectPar{:}); %'Plot',true);
    
    
    
    if (Iter<InPar.Niter)
        FlagCS   = FlagNN & StarsRelRMS<InPar.MaxRMS  & StarsRelRMS<Q & FlagRMS;
        FluxCS   = Flux(FlagCS,:);
    end
end

Res.MeanFlux    = MeanFlux;
Res.StarsRelRMS = StarsRelRMS;
Res.FlagCS      = FlagCS;

NormFlux        = Flux(FlagCS,:)./MeanFlux(FlagCS);
Res.rmsZP       = std(NormFlux,[],1);










