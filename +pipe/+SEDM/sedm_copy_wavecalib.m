function SI=sedm_copy_wavecalib(CalSI,SI)
%--------------------------------------------------------------------------
% sedm_copy_wavecalib function                                        SEDM
% Description: Copy wavelength calibration information from a
%              SegmentsInfo structure array of an Arc
%              (i.e., after sedm_wavecalib.m) to that of a science
%              SegmentsInfo structure array.
% Input  : - SegmentsInfo structure array of an Arc with wavelength
%            calibration information returned by sedm_wavecalib.m
%          - A SegmentsInfo structure array to copy the wavelength
%            calibration information into.
% Output : - A SegmentsInfo structure array with the wavelength
%            calibration information.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Sep 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: sedm_abswavecalib('SI',ArcSegmentsInfo(1130))
% Reliable: 2
%--------------------------------------------------------------------------

Fields2Copy = {'WaveCalib','Arc_BestCorr'};
Nf          = length(Fields2Copy);

Nseg = length(SI);
for Iseg=1:1:Nseg,
    for If=1:1:Nf,
       SI(Iseg).(Fields2Copy{If}) = CalSI(Iseg).(Fields2Copy{If});
    end
end

