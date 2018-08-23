function [ Data ] = getZodiac( lon,lat,filter_family,filter_name)
%
% newGetZodiac function                                                      LIM
% Description: Given the a position on the sky with respect to the zodiac
%               based on HST WFC3 handbook 
%              
% Input  : Position in Heliocentric coordinates
%          - Longtitude - [degrees]
%          - Latitude   - [degrees]
%          Spectral Band
%           -filter_family (see get_filter)
%           -filter_name
%           or
%           -Band upper limit in Angstrum
%           -Band lower limit in Angstrum
% Output : 
%          - AB mag in the new band at the position  
%          - # of photons / sec cm^2 within the new band
%          - color - 
%          - Magnitude in V-band in the position
% Tested : 
%     By : Ilan Sagiv                    Sep 2014
%    
% Example: ZodBack = newGetZodiac([180],[55],[2000],[2700])
%   

h=get_constant('h','cgs');
c=get_constant('c','cgs');
hc=h*c*1e8; %[erg]*[Ang]
% data from HST WFC3 handbook (table 9.3)
%flux units are erg/sec cm sec A arcsec^2
data=[
1000 2.41e-23 9.69e-29 2.41e-23
1100 4.38e-22 1.04e-26 4.38e-22
1200 4.01e-23 1.08e-25 4.03e-23
1300 7.41e-25 6.59e-25 1.40e-24
1400 4.29e-25 2.55e-24 2.98e-24
1500 4.16e-25 9.73e-24 1.01e-23
1600 2.55e-25 2.35e-22 2.35e-22
1700 7.89e-25 7.21e-21 7.21e-21
1800 9.33e-23 1.53e-20 1.54e-20
1900 4.39e-22 2.25e-20 2.29e-20
2000 1.01e-21 3.58e-20 3.68e-20
2100 1.60e-21 1.23e-19 1.25e-19
2200 7.49e-22 2.21e-19 2.22e-19
2300 3.32e-22 1.81e-19 1.81e-19
2400 2.50e-22 1.83e-19 1.83e-19
2500 2.39e-22 2.53e-19 2.53e-19
2600 5.62e-22 3.06e-19 3.06e-19
2700 6.77e-21 1.01e-18 1.02e-18
2800 2.03e-21 2.88e-19 2.90e-19
2900 4.32e-20 2.08e-18 2.12e-18
3000 9.34e-20 1.25e-18 1.35e-18
3100 2.07e-19 1.50e-18 1.70e-18
3200 3.60e-19 2.30e-18 2.66e-18
3300 4.27e-19 2.95e-18 3.38e-18
3400 6.40e-19 2.86e-18 3.50e-18
3500 8.20e-19 2.79e-18 3.61e-18
3600 1.06e-18 2.74e-18 3.80e-18
3700 1.22e-18 3.32e-18 4.54e-18
3800 1.23e-18 3.12e-18 4.35e-18
3900 1.52e-18 3.34e-18 4.86e-18
4000 2.38e-18 4.64e-18 7.01e-18
4250 2.38e-18 4.65e-18 7.03e-18
4500 2.86e-18 5.58e-18 8.44e-18
4750 2.79e-18 5.46e-18 8.25e-18
5000 2.63e-18 5.15e-18 7.77e-18
5250 2.67e-18 5.37e-18 8.04e-18
5500 2.58e-18 5.34e-18 7.92e-18
5750 2.54e-18 5.40e-18 7.94e-18
6000 2.42e-18 5.25e-18 7.67e-18
6250 2.26e-18 5.02e-18 7.28e-18
6500 2.17e-18 4.92e-18 7.09e-18
6750 2.07e-18 4.79e-18 6.87e-18
7000 1.93e-18 4.55e-18 6.48e-18
7250 1.85e-18 4.43e-18 6.29e-18
7500 1.74e-18 4.23e-18 5.97e-18
7750 1.63e-18 4.04e-18 5.67e-18
8000 1.56e-18 3.92e-18 5.49e-18
8250 1.48e-18 3.76e-18 5.23e-18
8500 1.35e-18 3.50e-18 4.85e-18
8750 1.31e-18 3.43e-18 4.74e-18
9000 1.22e-18 3.23e-18 4.44e-18
9250 1.15e-18 3.07e-18 4.21e-18
9500 1.10e-18 2.98e-18 4.08e-18
9750 1.04e-18 2.86e-18 3.91e-18
10000 1.00e-18 2.78e-18 3.78e-18
10250 9.54e-19 2.67e-18 3.63e-18
10500 9.04e-19 2.56e-18 3.46e-18
10750 8.41e-19 2.41e-18 3.25e-18
11000 8.03e-19 2.31e-18 3.11e-18];

wavelength=data(:,1);
Zodiac    =data(:,3);
Spectrum=[wavelength,Zodiac];

if(0)% plot Spectrum
	figure()
	semilogy(wavelength,Zodiac,'r');
	xlabel('Wavelength [A]')
	ylabel('Flux [erg / sec\cdotcm^2 \cdot A]')
    axis([1000 12e3 1e-25 1e-16])
end

EClon=0:15:180;
EClat=0:15:90;
% approximate zodiacal sky background as a function of
% Heliocentric ecliptic longitude and ecliptic latitude
% (in Vmag / arcsec^2)
% from HST WFC3 Handbook; (table 9.4)

%lon-->
%v lat
zmag=[...
nan ,nan ,nan ,nan ,22.6,23.0,23.3;
nan ,nan ,nan ,nan ,22.6,23.1,23.3;
nan ,nan ,nan ,22.3,22.7,23.1,23.3;
nan ,nan ,22.1,22.5,22.9,23.1,23.3;
21.3,21.9,22.4,22.7,23.0,23.2,23.3;
21.7,22.2,22.6,23.9,23.1,23.2,23.3;
22.0,22.3,22.7,23.0,23.2,23.3,23.3;
22.2,22.5,22.9,23.1,23.3,23.3,23.3;
22.4,22.6,22.9,23.2,23.3,23.3,23.3;
22.4,22.6,22.9,23.2,23.3,23.3,23.3;
22.4,22.6,22.9,23.1,23.3,23.3,23.3;
22.3,22.5,22.8,23.0,23.2,23.4,23.3;
22.1,22.4,22.7,23.0,23.2,23.4,23.3];

%%
if(ischar(filter_family))
    Filter = get_filter(filter_family,filter_name);
    bandcenter=Filter.eff_wl{1};
    Filter = Filter.nT{1};
else
    band_l=filter_family;
    band_h=filter_name;
    %make filter-------
    Filter=[(band_l:50:band_h)',0*(band_l:50:band_h)'+1];  %create filter
    Filter(:,2)=Filter(:,2)./trapz(Filter(:,1),Filter(:,2)); %normalize
    bandwidth=max(wavelength(wavelength<=band_h))...
             -min(wavelength(wavelength>=band_l));
    bandcenter=band_l+bandwidth/2;
end

%The V mag from the Zodiac Spectrum given for Vmag=22.1
% this is just a sanity check: gives 22.087
Vmag=synphot([wavelength,Zodiac],'Johnson','V','Vega');
%The difference in Vmag due to sun position
Vmag_out=interp2(EClat,EClon,zmag,lat,lon);
if(isnan(Vmag_out))
    Data='Too close to the Sun bubileh';
    return;
end
deltaMag=Vmag_out-Vmag;

mag=synphot([wavelength,Zodiac],Filter,'NUV','AB');
mag=mag+deltaMag;

color=deltaMag;

%get total number of photos in band
phFluxA=convert_flux(mag,'AB','ph/A',bandcenter,'A');
phFlux=trapz(Filter(:,1),Filter(:,2)/max(Filter(:,2))*phFluxA);

Data.pos                = [lat,lon];
Data.Mag                = mag;
Data.phFlux             = phFlux;
Data.color              = color;
Data.VMag               = Vmag_out;
end

