function Phase=phase_complex(Matrix)
% Return the phase of complex numbers.
% Package: telescope.Optics
% Description: Given an N-D array of complex numbers, return the phase of
%              each complex number
%              (i.e., atan2(imag(Matrix),real(Matrix))).
% Input  : - An array of complex numbers..
% Output : - The phases of the complex numbers.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Phase=phase_complex(Matrix);
% Reliable: 2
%--------------------------------------------------------------------------

Phase = atan2(imag(Matrix),real(Matrix));
