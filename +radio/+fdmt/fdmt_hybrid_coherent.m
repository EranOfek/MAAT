function [d_0, DMT, DispersionAxis]=fdmt_hybrid_coherent(raw_signal,N_p, f_min, f_max, D_max, SigmaBound)
% coherent hybrid - FDMT transform
% Package: radio.fdmt
% Description: Performs the coherent hybrid - FDMT transform for incoherent 
%            dedispersion. For further details and explanations,
%            See Zackay and Ofek 2015. Algorithm 3.
% Input  : - 1 dimensional signal, dimension should be a power of 2,
%            f_min, f_max - frequency band edges, with units of MHz
%            maxD         - Maximum dispersion measure, in pc*cm^-3
% Output : - The FDMT of the range of dispersions that finds the brightest
%            curve.
%            The used coarse trial dispersion for the CoherentDedispersion
%          - Dispersion measure time.
%          - DispersionAxis - the dispersion measure axis of the output.
%            The units for this vector are (pc cm^-3 dT^-1 s)
%            So, to get the answer in pc cm^-3 user needs to multiply the 
%            result by (time delay between subsequent spectra)/s (see
%            example).
% Reference: Zackay & Ofek 2015 (Algorithm 3).
% See also: fdmt.m, fdmt_fft.m
% Tested : Matlab R2014a
%     By : Barak Zackay                    Nov 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

import radio.fdmt.*

DispersionConstant = 4.148808*10^9;
ConversionConst = DispersionConstant *(1./f_min^2 - 1./f_max^2) * (f_max - f_min);

tau = 1/((f_max - f_min)*10^6);

dataType = 'uint32';

N_d = D_max * ConversionConst;

n_coherent = ceil(N_d/(N_p^2));
disp('number of coherent iterations:');
disp(n_coherent);

% In order to FFT only once...
ffted_signal = fft(raw_signal);

% In order to normalize the FDMT, we need to know how many numbers were
% summed to each dispersion bin. this is most easily traced by this method
FDMT_normalization = fdmt(ones([length(raw_signal)/N_p,N_p]),f_min,f_max,N_p, 'double');
for outer_d = 0:n_coherent-1,
    disp('coherent iteration'); 
    disp(outer_d);
    cur_coherent_d = outer_d .* (D_max./double(n_coherent));
    disp('cur_coherent_d = ');
    disp(cur_coherent_d);
    FDMT_input = abs(stfft(CoherentDedispersion(ffted_signal, cur_coherent_d, f_min, f_max, true), N_p)).^2;
    FDMT_input = FDMT_input - mean(FDMT_input(:));
    FDMT_input = FDMT_input ./ (0.25.*std(FDMT_input(:)));
    V = var(FDMT_input(:));
    [FDMT_output,DispersionAxisOutput] = fdmt(FDMT_input, f_min, f_max, N_p, dataType);
    FDMT_output = double(FDMT_output) ./ sqrt(double(FDMT_normalization .* V + 0.000001));
    if max(FDMT_output(:)) > SigmaBound,
            SigmaBound = max(FDMT_output(:));
            d_0 = cur_coherent_d;
            disp('achieved score with'); disp(SigmaBound); disp('sigmas');
            DMT = FDMT_output;
            DispersionAxis = DispersionAxisOutput * (tau * N_p) + d_0;
            
            
    end

end
