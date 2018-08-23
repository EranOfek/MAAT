function [DMT,DispersionAxis]=fdmt_fft(Image, f_min,f_max, maxDT, dataType)
% Incoherent Fast Disperssion Measure-FFT Transform (FDMT-FFT)
% Package: radio.fdmt
% Description: Performs the incoherent Fast Disperssion Measure-FFT
%              Transform (FDMT-FFT) (Zackay & Ofek 2015).
% Input  : - Two dimensional matrix [Time,Frequency]. Both dimensions
%            should be powers of 2, but not necessarily the same power.
%          - f_min - minimum frequency stored in Input(:,1) [MHz].
%          - f_max - maximum frequency stored in Input(:,end) [MHz].
%          - MaxDT - Maximum dispersion measure [units of bins].
%            Translation between Dt and dispersion measure (D) can be 
%            calculated by inverting the formula 
%            Dt * tau = D*d * (1/f_min^2 - 1/f_max^2).
%            Here, tau is time difference between adjacent spectra,
%            D is the dispersion constant (=4.15 ms Mhz^2 cm^3 pc^-1),
%            and d is the dispersion measure (in units of pc^cm^-3).
%          - The data type of the output.
% Output : - Dispersion measure transform of the input matrix.
%            The first dimension is time which corresponds to the 
%            arrival time of the lowest frequency (t0), while the second
%            dimension (Dt) corresponds to time delay between the lowest
%            and highest frequency.
%            The time have the same units as the input time.
%            t0 is in the range of XXX to XXX while Dt is in the range
%            of 0 to maxDT.
%            Each coordinate is the integral (sum) of the input over
%            a (t-t0) = K*(1/f^2 - 1/f^2) curve.
%            All curves that intersect the f_max and f_min axes are
%            calculated.  
%          - DispersionAxis - the dispersion measure axis of the output.
%            The units for this vector are (pc cm^-3 dT^-1 s)
%            So, to get the answer in pc cm^-3 user needs to multiply the 
%            result by (time delay between subsequent spectra)/s (see example)
% Reference: Zackay & Ofek 2015 (Algorithm 2).
% See also: fdmt.m, fdmt_hybrid_coherent.m
% Tested : Matlab R2014a
%     By : Barak Zackay                    Nov 2014
%    URL : https://sites.google.com/site/barakzackayhomepage/4-code-and-products
%          http://weizmann.ac.il/home/eofek/matlab/
% Example:  import radio.fdmt.*
%           f_min = 1200; f_max = 1600;
%           N_p = 2^10; N_s = 2^20;N_f = 2^10;N_d = 2^10;
%           tau = 1/ ((f_max - f_min)*10^6);
%           PulsePos = 123456;
%           PulseStrength = 1.5; % in sigmas/bin
%           d = 1; % pc cm^-3
%           signal = normrnd(0,1,[1,N_s]);
%           signal(:,PulsePos:PulsePos + N_p-1) = normrnd(0,PulseStrength + 1 ,[1,N_p]);
%           signal = CoherentDedispersion(signal,-d, f_min, f_max, false);
%           FDMT_input = abs(stfft(signal, N_f)).^2;
%           FDMT_input = FDMT_input - median(FDMT_input(:));
%           [FDMT_output,DispersionAxis] = fdmt_fft(FDMT_input, f_min, f_max, 1024, 'single');
%           DispersionAxis = DispersionAxis * (tau* N_f); % Converting the dispersionAxis to units of pc ^ cm^-3
%           surface(double(FDMT_output)); shading interp
% Reliable: 2
%--------------------------------------------------------------------------

[T,F] = size(Image);
[t,f] = deal(log2(T),log2(F));

if nargin == 4,
    dataType = 'single';
end

State3 = FDMT_initialization_FFT(Image,f_min,f_max,maxDT,dataType);

for i_t = 1:f,
    State3 = FDMT_iteration_FFT(State3,maxDT,F,f_min,f_max,i_t,dataType);
end

[T_tmp,F_tmp,dT_tmp] = size(State3);

DMT = reshape(ifft(State3,[],1),T_tmp,dT_tmp);

if nargout == 2,
    DispersionConstant = 4.148808 * 10^3 ; %Mhz^2 * s * pc^-1 * cm^3
    D_max = maxDT ./ ((1./f_min.^2 - 1./f_max.^2)*DispersionConstant);
    DispersionAxis = linspace(0,D_max,maxDT);
end


end


function [Output] = FDMT_iteration_FFT(Input,maxDT,F,f_min,f_max,iteration_num,dataType)
% internal function that does one iteration of dispersion measure transform on the Input
% Each iteration, the frequency dimension gets smaller and DM dimension
% gets larger

input_dims = size(Input);
output_dims = input_dims;

deltaF = 2^(iteration_num) * (f_max - f_min)/F;
dF = (f_max - f_min)/F;

% the maximum deltaT needed to calculate at the i'th iteration
deltaT = ceil(maxDT *(1/f_min^2 - 1/(f_min + deltaF)^2) / (1/f_min^2 - 1/f_max^2));

% number of subbands in the output is divided by 2
output_dims(2) = output_dims(2)/2;

% The number of DM's relevant to this iteration
output_dims(3) = deltaT + 1; 
Output = zeros(output_dims,dataType);

% No negative DM's are calculated 
% If you want negative DM's, this will have to change to 1+deltaT,
% 1+deltaTOld
ShiftOutput = 1;
ShiftInput = 1;
T = output_dims(1);
F_jumps = output_dims(2);

% No point in preparing this array for all shifts.
deltaTShift = ceil(maxDT *(1/f_min^2 - 1/(f_min + deltaF/2 + dF/2)^2) / (1/f_min^2 - 1/f_max^2))+2;

% If this turns out slow, or memory intensive, it is possible (and easy) to
% calclate this in a closed formula on the moment of multiplication
ShiftRow = fft(eye(T,deltaTShift),[],1);

% allowing for a shift to occur between layers
correction = dF/2;

for i_F = 1:F_jumps,
    
    % for convenience
    f_start = (f_max - f_min)/F_jumps * (i_F-1) + f_min;
    f_end = (f_max - f_min)/F_jumps *(i_F) + f_min;
    f_middle = (f_end - f_start)/2 + f_start;
    
    % to save work, each band's dT limit is calculated differently
    deltaTLocal = ceil(maxDT *(1/f_start^2 - 1/(f_end)^2) / (1/f_min^2 - 1/f_max^2));
    
    % going over all dT shifts in the output array
    for i_dT = 0:deltaTLocal,
        
        %calculating all the indices
        dT_middle = round(i_dT * (1/(f_middle - correction)^2 - 1/f_start^2)/(1/f_end^2 - 1/f_start^2));
        dT_middle_index = dT_middle + ShiftInput;
        dT_middle_larger = round(i_dT * (1/(f_middle + correction)^2 - 1/f_start^2)/(1/f_end^2 - 1/f_start^2));

        dT_rest = i_dT - dT_middle_larger;
        dT_rest_index = dT_rest + ShiftInput;
        
        % Addition rule of FDMT-FFT
        if dT_middle_larger == 0,
            Output(:,i_F,i_dT + ShiftOutput) = Input(:,2*i_F-1, dT_middle_index) + Input(:,2*i_F, dT_rest_index);
        else
            Output(:,i_F,i_dT + ShiftOutput) = Input(:,2*i_F-1, dT_middle_index) + Input(:,2*i_F, dT_rest_index) .* ShiftRow(:,dT_middle_larger+1);
        end    
        
    end
end

end


function [Output] = FDMT_initialization_FFT(Image,f_min,f_max, maxDT, dataType)
% Internal function that takes care at:
% 1) Build the DM dimension, copies the data
% 2) fill dT!=0 bins with the proper partial sum of shifts within a freq
% bin
% 3) fft'ing the time axis, as required by the FDMT-FFT algorithm.
[T,F] = size(Image);

deltaF = (f_max - f_min)/F;
deltaT = ceil(maxDT *(1/f_min^2 - 1/(f_min + deltaF)^2) / (1/f_min^2 - 1/f_max^2));

Output = zeros([T,F,deltaT+1],dataType);
Output(:,:,1) = fft(Image,[],1);
I = fft(cumsum(eye(T,deltaT+1),2),[],1);

% shifting is like complex multiplication with the FFT of the eye rows => 
% cumsum is like complex multiplication with the FFT of the cumsum of the eye rows 
for i_dT = 2:deltaT+1,
    Output(:,:,i_dT) = bsxfun(@times,Output(:,:,1), I(:,i_dT));
end


end


