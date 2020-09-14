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
%            'BinMethod' - 'Fit' : Fit a polynomial to the flux vs. rms.
%                                  This is more reliable when only small
%                                  number of sources are available.
%                          'Flux': bin by constant flux bins.
%                          'Number' (default) - by constant number of
%                                  sources in bin.
%            'PolyOrder' - Order of polynomial fitted when 'BinMethod' is
%                          'Fit'.
%            'PolyFluxLog' - {true|false} When using the 'Fit' BinMethod,
%                            apply log10 to flux before fit. Default is
%                            true.
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
%            'Plot'    - Plot rms vs. flux with selected and unselected
%                        sources.
% Output : - Flag of sources which RMS is below the rms vs. flux envelope.
%          - Structure containing additional data:
%            'Bin' - Matrix of binned rms vs, flux.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Flag,Res]=timeseries.rmsflux_select(Flux,RMS);
%          % A full example:
%          MeanFlux = rand(100,1).*1e5;
%          Norm     = 1+0.05.*randn(1,200); 
%          Flux     = poissrnd(MeanFlux.*Norm);
%          [CalFlux,Res] = timeseries.calib_flux_sum(Flux);
%          loglog(Res.MeanFlux,Res.StarsRelRMS,'.')
%          [Flag,ResS]=timeseries.rmsflux_select(Res.MeanFlux,Res.StarsRelRMS,'Plot',true);
%          set(gca,'XS','log','YS','log')
% Reliable: 2
%--------------------------------------------------------------------------

DefV.BinMethod            = 'number';  % 'fit' | 'Flux' | 'Number'
DefV.PolyOrder            = 2;
DefV.PolyFluxLog          = true;
DefV.NinBin               = 30;
DefV.FluxBinSize          = 0.5;  % usually in magnitude
DefV.MinInBin             = 3;
DefV.MinGoodBins          = 2;
DefV.SigClip              = 3;
DefV.StdFun               = @Util.stat.rstd;
DefV.MeanFun              = @nanmedian;
DefV.InterpMethod         = 'linear';
DefV.Plot                 = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);



N = numel(Flux);

switch lower(InPar.BinMethod)
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
        
        Envelope =  Bin(FlagGoodBins,3) + InPar.SigClip.*Bin(FlagGoodBins,4);
        InterpEnvelope = interp1(Bin(FlagGoodBins,1),Envelope,Flux,InPar.InterpMethod);
        % extrapolate using edges
        InterpEnvelope(isnan(InterpEnvelope) & Flux<Bin(1,1))   = Envelope(1);
        InterpEnvelope(isnan(InterpEnvelope) & Flux>Bin(end,1)) = Envelope(end);
        
        if (sum(FlagGoodBins)<InPar.MinGoodBins)
            error('Not enough good bins');
        end

    case 'flux'
        Bin = timeseries.binning([Flux,RMS],InPar.FluxBinSize,[],{'MidBin',@numel,InPar.MeanFun,InPar.StdFun});

        % interpolate over bad bins (too few sources)
        FlagGoodBins = Bin(:,2)>=InPar.MinInBin;
        
        Envelope =  Bin(FlagGoodBins,3) + InPar.SigClip.*Bin(FlagGoodBins,4);
        InterpEnvelope = interp1(Bin(FlagGoodBins,1),Envelope,Flux,InPar.InterpMethod);
        % extrapolate using edges
        InterpEnvelope(isnan(InterpEnvelope) & Flux<Bin(1,1))   = Envelope(1);
        InterpEnvelope(isnan(InterpEnvelope) & Flux>Bin(end,1)) = Envelope(end);

        if (sum(FlagGoodBins)<InPar.MinGoodBins)
            error('Not enough good bins');
        end

        
    case 'fit'
        
        if (InPar.PolyFluxLog)
            PolyP = polyfit(log10(Flux),RMS,InPar.PolyOrder);
            Y     = polyval(PolyP,log10(Flux));
        else
            PolyP = polyfit(Flux,RMS,InPar.PolyOrder);
            Y     = polyval(PolyP,Flux);
        end
        Std   = InPar.StdFun(RMS-Y);
        InterpEnvelope = Y + InPar.SigClip.*Std;
        
        Bin = [];
        
    otherwise
        error('Unknown Method option');
end




Flag = RMS<InterpEnvelope;

Res.Bin = Bin;

if (InPar.Plot)
    plot(Flux(Flag),RMS(Flag),'k.');
    hold on;
    plot(Flux(~Flag),RMS(~Flag),'.');
    
    switch lower(InPar.BinMethod)
        case 'fit'
            [~,SI] = sort(Flux);
            plot(Flux(SI),InterpEnvelope(SI),'k-','Color',[0.8 0.8 0.8]);
            plot(Flux(SI),10.^polyval(PolyP,Flux(SI)),'k-');
        otherwise
            plot(Bin(FlagGoodBins,1),Envelope,'k-','Color',[0.8 0.8 0.8]);
            plot(Bin(:,1),Bin(:,3),'k-');
    end
    set(gca,'YS','log');
    set(gca,'XS','log');
end








