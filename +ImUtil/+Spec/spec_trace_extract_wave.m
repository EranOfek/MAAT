function Info=spec_trace_extract_wave(Sim,StartPos,SimArc,varargin)
%--------------------------------------------------------------------------
% spec_trace_extract_wave function                                  ImSpec
% Description: A wrapper around spec_trace.m, spec_extract 2/1-d and
%              automatic wavelength calibration.
% Input  : - Sim image of spectra to extract.
%          - Start position of trace [X,Y].
%          - Simimage of corresponding wavelength calibration arc image.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'GoodRange'   - Good range of pixels [min max] along the
%                            dispersion axis. If empty will attempt
%                            to find automatically. Default is empty.
%            'MasterTrace' - Master trace (see spec_trace.m).
%                            Default is empty.
%            'BackRegion'  - background region (see spec_skysub.m),
%                            relative to peak position.
%            'WaveCalibType'- Cell array of types of wave calib to
%                            perform. Default is {'arc','sky'}.
%            'WaveCalibSourceSky' - Source for sky wave calib.
%                            Default is 'SkyLow'.
%            'WaveCalibSourceArc' - Source for arc wave calib.
%                            Default is 'Ne+Ar'.
%            'MinDN'         - Minimum number of counts in trace.
%                            Regions near the edges with
%                            counts below this value will be
%                            ignored.
%            'ThreshBackRMS' - Threshold for background RMS.
%                            Default is 10. Regions near the edges with
%                            background RMS above this value will be
%                            ignored.
%            'MaxTraceRange' - maximum allowed range in trace Y axis.
%                            Default is 100. If range is larger than
%                            the source will not be extracted.
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
%            spec_trace.m, spec_extract_2d.m, spec_skysub.m, 
%            spec_fit1d_psf.m, spec_wavecalib_xcorr.m,
%            spec_wavecalib_lines.m
% Output : - A structure containing the trace, extraction, and wave
%            calib information.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%-----------------------------------------------------------------------------


ImageField  = 'Im';
%HeaderField = 'Header';
%FileField   = 'ImageFileName';
%MaskField   = 'Mask';
BackImField = 'BackIm';
%ErrImField  = 'ErrIm';

DefV.GoodRange          = [];
DefV.MasterTrace        = [];
DefV.BackRegion         = [];
DefV.WaveCalibType      = {'arc','sky'};
DefV.WaveCalibSourceSky = 'SkyLow';
DefV.WaveCalibSourceArc = 'Ne+Ar';
DefV.ThreshBackRMS      = 10;
DefV.MinDN              = 10;
DefV.ExtSemiW           = 50;
DefV.MaxTraceRange      = 100;


InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

  %AllPeaks{Is}(Ipeak)
          
           %StartPos = [MeanCollapseRegion,AllPeaks{Is}(Ipeak).X];
           % populate Sim
           %Sim(Is).AllTraces(Ipeak).Peak = AllPeaks{Is}(Ipeak);
           
           %BackRegion = AllPeaks{Is}(Ipeak).Back-AllPeaks{Is}(Ipeak).X+InPar.ExtSemiW;

           
Info.StartPos = StartPos;
           
           
% Trace the spectrum
Info.Trace    = spec_trace(Sim,Info.StartPos,varargin{:});

% check if Trace failed
if (range(Info.Trace.Y)>InPar.MaxTraceRange || range(Info.Trace.SmY)>InPar.MaxTraceRange),
    % failed
    Info.ExtractedSpec = [];
    Info.BackSubSpec   = [];
    Info.BackInf       = [];
    Info.FitPSF        = [];
    Info.FitWave       = [];
    Info.MatchedList   = [];
else
    % Extract the spectrum
    Info.ExtractedSpec = spec_extract_2d(Sim,Info.Trace,varargin{:});


    % Subtract background
    [Info.BackSubSpec,Info.BackInfo]=spec_skysub(Info.ExtractedSpec,varargin{:});

    % Extract 1D spectrum
    Info.FitPSF = spec_fit1d_psf(Info.BackSubSpec.(ImageField),Info.BackSubSpec.(BackImField),varargin{:},'Mask',Info.BackSubSpec.Mask);

    % wavelength calibration observation
    % search continus regions near the edge which are unreliable
    if (isempty(InPar.GoodRange)),
        SmBackRMS = medfilt1(Info.BackInfo.RMS,101);
        %plot(SmBackRMS)
        MinPix = find(SmBackRMS<InPar.ThreshBackRMS,2,'first');
        MinPix = MinPix(2);
        MaxPix = find(SmBackRMS<InPar.ThreshBackRMS,2,'last');
        MaxPix = MaxPix(1);
    else
       MinPix = InPar.GoodRange(1);
       MaxPix = InPar.GoodRange(2);
    end

    Nwc = numel(InPar.WaveCalibType);
    for Iwc=1:1:Nwc,
        InPar.WaveCalibType{Iwc}
        switch lower(InPar.WaveCalibType{Iwc})
             case 'sky'
                  % background vector
                  WaveCalibSource = InPar.WaveCalibSourceSky;
                  BackVec = Info.BackSubSpec.(BackImField)(InPar.ExtSemiW+1,:);
             case 'arc'
                  % arc vector along trace

                  %IndArc = find(IsArc,1,'first')
                  ArcSpec = [Info.Trace.X, interp2(SimArc.(ImageField),Info.Trace.X,Info.Trace.SmY,'linear')];
                  
                  WaveCalibSource = InPar.WaveCalibSourceArc; 
                  BackVec = ArcSpec(:,2).';

            otherwise
                  error('Unknwon WaveCalibType option');
        end

        BackVec(1:MinPix) = 0;
        BackVec(MaxPix:end) = 0;

        [FitWave,MatchedList] = spec_wavecalib_xcorr(BackVec.',WaveCalibSource,varargin{:});
    
        Info.FitWave.(InPar.WaveCalibType{Iwc}) = FitWave;
        Info.MatchedList.(InPar.WaveCalibType{Iwc}) = MatchedList;
    end
end
