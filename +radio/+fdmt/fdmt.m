function [DMT, DispersionAxis]=fdmt(Input, f_min,f_max, maxDT, dataType)
% Incoherent Fast Disperssion Measure Transform (FDMT)
% Package: radio.fdmt
% Description: Performs the incoherent Fast Disperssion Measure
%              Transform (FDMT; Zackay & Ofek 2015).
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
% Reference: Zackay & Ofek 2015 (Algorithm 1).
% See also: fdmt_fft.m, fdmt_hybrid_coherent.m
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
%           [FDMT_output,DispersionAxis] = fdmt(FDMT_input, f_min, f_max, 1024, 'single');
%           DispersionAxis = DispersionAxis * (tau* N_f); % Converting the dispersionAxis to units of pc ^ cm^-3
%           surface(double(FDMT_output)); shading interp
% Reliable: 2
%--------------------------------------------------------------------------

[T,F] = size(Input);
[~,f] = deal(log2(T),log2(F));

if nargin == 4,
    dataType = 'single';
end


% FDMT initialization
State3 = FDMT_initialization(Input,f_min,f_max,maxDT,dataType);

% main iteration - log N iterations
for i_t = 1:f,
    State3 = FDMT_iteration(State3,maxDT,F,f_min,f_max,i_t,dataType);
end

[T_tmp,~,dT_tmp] = size(State3);
DMT= reshape(State3,T_tmp,dT_tmp);

if nargout == 2,
    DispersionConstant = 4.148808 * 10^3 ; %Mhz^2 * s * pc^-1 * cm^3
    D_max = maxDT ./ ((1./f_min.^2 - 1./f_max.^2)*DispersionConstant);
    DispersionAxis = linspace(0,D_max,maxDT);
end

end


function [Output] = FDMT_iteration(Input,maxDT,F,f_min,f_max,iteration_num,dataType)
% internal function that does one iteration of dispersion measure transform on the Input
% Each iteration, the frequency dimension gets smaller and DM dimension
% gets larger

input_dims = size(Input);
output_dims = input_dims;

% frequency difference between adjacent frequency sub-bands
deltaF = 2^(iteration_num) * (f_max - f_min)/F;

% frequency difference between adjacent frequency bins
dF = (f_max - f_min)/F;

% the maximum deltaT needed to calculate at the i'th iteration
deltaT = ceil(maxDT *(1/f_min^2 - 1/(f_min + deltaF)^2) / (1/f_min^2 - 1/f_max^2));

% frequency dimension shrinks by a factor of 2
output_dims(2) = output_dims(2)/2;

% The size of the DM simension
output_dims(3) = deltaT + 1; 

% Allocate the output array
Output = zeros(output_dims,dataType);

% Matlab counts from one
ShiftOutput = 1;
ShiftInput = 1;
T = output_dims(1);
F_jumps = output_dims(2);

% allowing for a shift to occur between subbands
correction = dF/2;
%counter = 0;
for i_F = 1:F_jumps,
    f_start = (f_max - f_min)/F_jumps * (i_F-1) + f_min;
    f_end = (f_max - f_min)/F_jumps *(i_F) + f_min;
    f_middle = (f_end - f_start)/2 + f_start;
        
    deltaTLocal = ceil(maxDT *(1/f_start^2 - 1/(f_end)^2) / (1/f_min^2 - 1/f_max^2));
    for i_dT = 0:deltaTLocal,
        
        % determining all the needed shift constants for the iteration
        dT_middle = round(i_dT * (1/(f_middle - correction )^2 - 1/f_start^2)/(1/f_end^2 - 1/f_start^2));
        dT_middle_index = dT_middle + ShiftInput;
        dT_middle_larger = round(i_dT * (1/(f_middle + correction )^2 - 1/f_start^2)/(1/f_end^2 - 1/f_start^2));

        dT_rest = i_dT - dT_middle_larger;
        dT_rest_index = dT_rest + ShiftInput;
        
        i_T_min = 1;
        i_T_max = dT_middle_larger+1;
        
        % Alternative addition rule!
        Output(i_T_min:i_T_max,i_F,i_dT + ShiftOutput) = Input(i_T_min:i_T_max,2*i_F-1, dT_middle_index);
        
        
        i_T_min = dT_middle_larger+2;
        i_T_max = T;
            
        % Addition rule!
        Output(i_T_min:i_T_max,i_F,i_dT + ShiftOutput) = Input(i_T_min:i_T_max,2*i_F-1, dT_middle_index) + Input(i_T_min - dT_middle_larger:i_T_max-dT_middle_larger,2*i_F, dT_rest_index);
                
    
    end
end
%counter
end

function [Output] = FDMT_initialization(Image,f_min,f_max,maxDT,dataType)
% Internal function that takes care at:
% 1) Build the DM dimension, copies the data
% 2) fill dT!=0 bins with the proper partial sum of shifts within a freq
% bin
[T,F] = size(Image);

deltaF = (f_max - f_min)/F;
deltaT = ceil(maxDT *(1/f_min^2 - 1/(f_min + deltaF)^2) / (1/f_min^2 - 1/f_max^2));

Output = zeros([T,F,deltaT+1],dataType);
ExtendedInput = zeros(T+deltaT+1,F,dataType);
ExtendedInput(1:T,:) = Image;
Output(:,:,1) = Image;

for i_dT = 2:deltaT+1,
    Output(:,:,i_dT) = Output(:,:,i_dT-1) + ExtendedInput(i_dT:i_dT + T-1,:);
end

end

