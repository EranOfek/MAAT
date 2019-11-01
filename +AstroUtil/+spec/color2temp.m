function temp = color2temp(color, filter1, filter2, filter_system, mag_system)
% Usage: temp = color2temp(color, filter1, filter2, filter_system='GAIA', mag_system='AB')
% Use the difference between two color magnitudes to estimate the
% temperature, using Eran's AstroUtil.spec.blackbody_mag_c.
%
% This function tries to reverse the blackbody by minimizing the difference
% between the given color term and the one you get from blackbody_mag_c. 
%
% Inputs: -The color magnitude difference, e.g., B-V
%         -The first filter (in the above example, B). 
%         -The second filter (in the above example, V).
%         -Filter system (see get_filter.m for details). Default is GAIA. 
%         -The magnitude system used for the input color. Default is AB. 
%
% Output: the temperature that most closely resembles the color given. 
%
% Example : temp = color2temp(1.3, 'B', 'V', 'GAIA', 'AB');
%
% Written by: Guy Nir 31/10/2019

    import AstroUtil.spec.blackbody_mag_c;
    
    if nargin==0, help('AstroUtil.spec.color2temp'); return; end
    
    if nargin<4 || isempty(filter_system)
        filter_system = 'GAIA';
    end
    
    if nargin<5 || isempty(mag_system)
        mag_system = 'AB';
    end
    
    func = @(T) abs(blackbody_mag_c(T, filter_system, filter1, mag_system) - blackbody_mag_c(T, filter_system, filter2, mag_system) - color);

    temp = fminsearch(func, 5000, optimset('TolX', 1)); 


end