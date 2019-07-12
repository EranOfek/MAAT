function [Flag,Res]=rmsflux_select(Flux,RMS,varargin)
% Select points in an flux/mag vs. rms like diagram
% Package: timeseries
% Description: Select points in an rms vs. flux like diagram.
%              Bin the flux/mag vs. rms diagram by flux bins or by constant
%              number of points per bin. Calculate the median and rstd in
%              each bin and flag sources which are below the median +
%              N*rstd, where N is the 'SigClip' parameter.
%              I.e., good sources are flaged.
% Input  : - A vector of flux/mag -like parameter.
%          - A vector of rms-like parameter.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'BinMethod' - 'Flux': bin by constant flux bins.
%                           'Number' (default) - by constant number of
%                           sources in bin.
%            'NinBin' - Number in bin (default is 30) for the 'Number'
%                           method.
%            'FluxBinSize' - Bin size for the 'Flux' method. Default is 0.5
%                           (good for magnitude units).
%            'MinInBin' - Minimum niber of sources in bin. Default is 3.
%            'MinGoodBins' - Required number of good bins (that contains
%                           more than 'MinInBin' stars. Default is 2.
%                           If not will end with error.
%            'SigClip' - Sigma clipping parameter. Default is 3.
%            'StdFun'  - Default is @Util.stat.rstd
%            'MeanFun' - Default is @nanmedian
%            'InterpMethod' - Default is 'linear'
% Output : - Flag of sources which RMS is below the rms vs. flux envelope.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Flag,Res]=timeseries.rmsflux_select(Flux,RMS);
% Reliable: 
%--------------------------------------------------------------------------



DefV.BinMethod            = 'Number';  % 'Flux' | 'Number'
DefV.NinBin               = 30;
DefV.FluxBinSize          = 0.5;  % usually in magnitude
DefV.MinInBin             = 3;
DefV.MinGoodBins          = 2;
DefV.SigClip              = 3;
DefV.StdFun               = @Util.stat.rstd;
DefV.MeanFun              = @nanmedian;
DefV.InterpMethod         = 'linear';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);



N = numel(Flux);

switch lower(InPar.Method)
    case 'number'
        SortedFlux = sort(Flux);

        Vec = (1:InPar.NinBin:N);
        if (Vec(end)<N)
            Vec = [Vec, N];
        end
        Edges = SortedFlux(Vec);
        Edges(1)   = Edges(1)   - eps.*10;
        Edges(end) = Edges(end) + eps.*10;

        Bin = timeseries.binning([Flux,RMS],Edges,[],{'MidBin',@numel,InPar.MeanFun,InPar.StdFun});
        Nbin = size(Bin,1);
        
        FlagGoodBins = true(Nbin,1);
    case 'flux'
        Bin = timeseries.binning([Flux,RMS],InPar.FluxBinSize,[],{'MidBin',@numel,InPar.MeanFun,InPar.StdFun});

        % interpolate over bad bins (too few sources)
        FlagGoodBins = Bin(:,2)>=InPar.MinInBin;
        
    otherwise
        error('Unknown Method option');
end

if (sum(FlagGoodBins)<InPar.MinGoodBins)
    error('Not enough good bins');
end

Envelope =  Bin(FlagGoodBins,3) + InPar.SigClip.*Bin(FlagGoodBins,4);
InterpEnvelope = interp1(Bin(FlagGoodBins,1),Envelope,Flux,InPar.InterpMethod);
% extrapolate using edges
InterpEnvelope(isnan(InterpEnvelope) & Flux<Bin(1,1))   = Envelope(1);
InterpEnvelope(isnan(InterpEnvelope) & Flux>Bin(end,1)) = Envelope(end);

Flag = RMS<InterpEnvelope;











