function [SegmentsInfo]=sedm_plot_seg_int(SegmentsInfo,varargin)
%--------------------------------------------------------------------------
% sedm_plot_seg_int function                                          SEDM
% Description: Calculate for each spexcell its median value within
%              a givenb wavelength range and store it in SegmentsInfo,
%              and plot the median intensity of each spexcell as a function
%              of its mean X/Y position.
% Input  : - SegmentsInfo.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'WaveRange' - Wavelength range [nm] in which to calculate
%                          the median. Defauly is [450 850].
%            'SpecField' - From which field to get the intensity.
%                          Default is 'SpexSpecFit'.
%            'MeanFun'   - Median function. Default is @median.
%            'MarkerSize'- Marker size to use in scatter plot.
%                          Default is 30.
%            'Plot'      - Plot {true|false}. Default is true.
% Output : - SegmentsInfo with an additional field 'SpecMedianVal'.
% Tested : Matlab R2013a
%     By : Eran O. Ofek                    Aug 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: sedm_plot_seg_int(SegmentsInfoIm);
% Reliable: 2
%--------------------------------------------------------------------------



DefV.WaveRange  = [450 850];
DefV.SpecField  = 'SpexSpecFit';
DefV.MeanFun    = @median;
DefV.MarkerSize = 30;
DefV.Plot       = true;
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

Nseg = numel(SegmentsInfo);
for Iseg=1:1:Nseg,
     FlagW = SegmentsInfo(Iseg).WaveCalib>InPar.WaveRange(1) & SegmentsInfo(Iseg).WaveCalib<InPar.WaveRange(2);
     
     SegmentsInfo(Iseg).SpecMedianVal   = InPar.MeanFun(SegmentsInfo(Iseg).(InPar.SpecField)(FlagW));
end
     
if (InPar.Plot),
    H = scatter([SegmentsInfo.MeanX].',[SegmentsInfo.MeanY].',InPar.MarkerSize,[SegmentsInfo.SpecMedianVal].','filled');
    set(H,'Marker','p')
    
end
     
