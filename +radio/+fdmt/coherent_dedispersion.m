function CD=coherent_dedispersion(raw_signal,d, f_min, f_max, alreadyFFTed)
% Perform coherent dedispersion using convolution
% Package: radio.fdmt 
% Description: Perform coherent dedispersion using convolution.
% Input :  - raw signal - is assumed to be a one domensional signal
%            of voltage measurments as a function of time.
%          - d - is the dispersion measure for which to perform the
%            coherent dedispersion. units: pc*cm^-3.
%          - f_min - the minimum freq, given in Mhz.
%          - f_max - the maximum freq, given in Mhz.
%          - A flag (true|false) indicating if the raw_signal is already
%            FFTed (to reduce complexity).
% Output : - Coherent dedispersion.
% Notes  : For future improvements:
%          1) Signal partition to chunks of length N_d is not applied.
%          2) No use of packing is done.
%            (either in the coherent stage (and take special care of
%             the abs()^2 operation done by other functions) or in the
%             incoherent stage).
% Tested : Matlab R2014a
%     By : Barak Zackay                    Nov 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------
    DispersionConstant = 4.148808*10^9;

    N_total = length(raw_signal);
    
    practicalD = DispersionConstant .* d;
    f = linspace(0,f_max-f_min,N_total);
    
    %# The added linear term makes the arrival times of the highest frequencies be 0
    H = exp(-(2.*pi.*1i .* practicalD ./(f_min + f) + 2.*pi.*1i .* practicalD.*f ./(f_max.^2)));
    
    if alreadyFFTed,
        CD = ifft(raw_signal .* H);
    else
        CD = ifft(fft(raw_signal) .* H);
    end
end