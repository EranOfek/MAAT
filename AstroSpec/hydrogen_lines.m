function L=hydrogen_lines(M,N);
%------------------------------------------------------------------------------
% hydrogen_lines function                                            AstroSpec
% Description: Calculate the vacum wavelength of Hydrogen lines, given their
%              shell numbers.
% Input  : - Inner shell number (scalar or a vector).
%            (e.g., 2 for Balmer series).
%          - Outer shell number (scalar or a vector).
% Output : - Wavelengths of the corresponding Balmer lines [A].
% Tested : Matlab 7.3
%     By : eran O. Ofek                       Feb 2008
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: L=hydrogen_lines(2,3);  % Balmer H\alpha line
% Reliable: 1
%------------------------------------------------------------------------------

B = 3645.6;  % A
R = 10967758.341.*1e-10;  % [A^-1] Rydberg constant

L = 1./( R .* (1./(M.^2) - 1./(N.^2)));
