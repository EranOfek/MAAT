function [Itot]=fresnel_occultation_filt(thetaStar,radiusOcculter,distanceToOcculter,vectorOfSamplePoints,impactParameter,spatialRes, spectrum, deltaLambda)
% Polychromatic difftraction pattern by finite soirce
% Package: AstroUtil.Occultation
% Description: Calculates the diffraction pattern caused by a finite
%              or point polychromatic source.
% Input  : - Angular size of the star [milli-arcseconds]. Default is 0.
%          - Radius of the occulting object [km]. Default is 1.
%          - Distance to the occulting object [AU]. Default is 40.
%          - A one dimensional row vector specifies the x-coordinates of
%            the occulter trajectory in which to calculate the diffraction.
%            [Fresnel radius units]. Default is (0:0.1:10).'.
%          - The impact parameter [km], in the plane of the occulter.
%            I.e., The minimum distance from the source center to the
%            occulter trajectory in km. Default is 0.
%          - The integration grid spatial resolution, in units of the
%            stellar radius. Default is 0.1. For a point source, enter 0.
%          - Either a two element cell array: {'family', 'band'} to be
%            used with the get_filter.m function, or a Nx2 ,matrix which
%            its columns are a wavelength in Angstrom and the corresponding 
%            transmission in that wavelength. Default is {'SDSS','r'}.
%          - The desired wavelength sampling rate (Angstrom). Default is
%            100.
%          - (CURRENTLY UNUSED) the componenet of the relative earth-KBO
%            velocity vector in units of the Fresnel scale. Default is 0.
% Output : - The total intensity, the sum of the intensity of each source.
% Tested : Matlab R2011b
%     By : Lior Rubanenko                  Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: fresnel_occultation.m; fresnel_occultation_fs.m
% Examples: I=AstroUtil.Occultation.fresnel_occultation_filt; plot((0:0.1:10),I)
% Reliable: 2
%--------------------------------------------------------------------------
RAD = 180./pi;

import AstroUtil.Occultation.*

AU = get_constant('au');

Def.thetaStar             = 0.01;
Def.radiusOcculter        = 1;
Def.distanceToOcculter    = 40;
Def.vectorOfSamplePoints  = (0:0.1:10);
Def.impactParameter       = 0;
Def.spatialRes            = 0.1;
Def.spectrum              = {'SDSS','r'};
Def.deltaLambda           = 100;

if (nargin==0),
    thetaStar              = Def.thetaStar;
    radiusOcculter         = Def.radiusOcculter;
    distanceToOcculter     = Def.distanceToOcculter;
    vectorOfSamplePoints   = Def.vectorOfSamplePoints;
    impactParameter        = Def.impactParameter;
    spatialRes             = Def.spatialRes;
    spectrum               = Def.spectrum;
    deltaLambda            = Def.deltaLambda;
elseif (nargin==1),
    radiusOcculter         = Def.radiusOcculter;
    distanceToOcculter     = Def.distanceToOcculter;
    vectorOfSamplePoints   = Def.vectorOfSamplePoints;
    impactParameter        = Def.impactParameter;
    spatialRes             = Def.spatialRes;
    spectrum               = Def.spectrum;
    deltaLambda            = Def.deltaLambda;
elseif (nargin==2),
    distanceToOcculter     = Def.distanceToOcculter;
    vectorOfSamplePoints   = Def.vectorOfSamplePoints;
    impactParameter        = Def.impactParameter;
    spatialRes             = Def.spatialRes;
    spectrum               = Def.spectrum;
    deltaLambda            = Def.deltaLambda;
elseif (nargin==3),
    vectorOfSamplePoints   = Def.vectorOfSamplePoints;
    impactParameter        = Def.impactParameter;
    spatialRes             = Def.spatialRes;
    spectrum               = Def.spectrum;
    deltaLambda            = Def.deltaLambda;
elseif (nargin==4),
    impactParameter        = Def.impactParameter;
    spatialRes             = Def.spatialRes;
    spectrum               = Def.spectrum;
    deltaLambda            = Def.deltaLambda;
elseif (nargin==5),
    spatialRes             = Def.spatialRes;
    spectrum               = Def.spectrum;
    deltaLambda            = Def.deltaLambda;
elseif (nargin==6),
    spectrum               = Def.spectrum;
    deltaLambda            = Def.deltaLambda;
elseif (nargin==7),
    deltaLambda            = Def.deltaLambda;
elseif (nargin==8),
    %   Do nothing
else
    error('Illegal number of input arguments');
    
end

%%  Prepare input:
if iscell(spectrum)
    filter = get_filter(spectrum{1}, spectrum{2});
    wavelength = filter.nT{1}(:,1);
    transmission = filter.nT{1}(:,2);
    
elseif  isnumeric(spectrum)
    if numel(spectrum) == 1
        error('Input must be a cell array or vector.');
    else
        wavelength  = spectrum(:,1);
        transmission = spectrum(:,2);
    end
else
    error('Illegal filter transmission format');
end


%%  Get transmission function T(lambda):
%   Verify input:
if (wavelength(end) < wavelength(1))
    wavelength = flipud(wavelength);
    transmission = flipud(transmission);
end



selectedWavelengths = wavelength(1):deltaLambda:wavelength(end);
T = interp1(wavelength, transmission, selectedWavelengths);

%%  Calculate a vector containing the Fresnel scale that corresponds to each selected wavelength:
%   Convert everything to micrometers:
distanceToOcculterInMicrometers = distanceToOcculter .* AU.*1e4;
selectedWavelengthsInMicrometers = selectedWavelengths .* 1e-4;
radiusOfOcculterInMicrometer = radiusOcculter.*1e9;   % km to micrometer
impactParameter = impactParameter.*1e9;


F = sqrt(selectedWavelengthsInMicrometers .* distanceToOcculterInMicrometers .*0.5);

%   Fresnel scale in units of milliarcseconds:
F_as = (F ./ distanceToOcculterInMicrometers) .* RAD.*3600.*1000;

%%  Calculate the diffraction:
rhoOcculter = radiusOfOcculterInMicrometer ./ F;

if thetaStar == 0 %Point source:
    I_rho = zeros(length(vectorOfSamplePoints), length(rhoOcculter));
    impactParameter = impactParameter ./ F;
    
    for Iw=1:length(selectedWavelengths),
        vectorOfDistances = sqrt(vectorOfSamplePoints.^2 + impactParameter(Iw).^2);
        I_rho(:,Iw) = fresnel_occultation_ps(rhoOcculter(Iw), vectorOfDistances);
    end  
    
    %   Create a matrix of the different fractions of I, according to the
    %   transmission function:
    Transmission_poly = bsxfun(@times, T, I_rho);
    Itot = sum(Transmission_poly, 2).*deltaLambda;
        
else
    %%  Set physical parameters:
    rhoStar = thetaStar ./ F_as;
    rhoOcculter = radiusOfOcculterInMicrometer ./ F;
    I_rho = zeros(length(vectorOfSamplePoints), length(rhoOcculter));

    impactParameter = impactParameter ./ F;
        
    for Iw=1:length(selectedWavelengths),
        
        I_rho(:,Iw) = fresnel_occultation_fs(rhoStar(Iw), rhoOcculter(Iw), vectorOfSamplePoints, impactParameter(Iw), spatialRes);
    end  
    
    %   Create a matrix of the different fractions of I, according to the
    %   transmission function:
    Transmission_poly = bsxfun(@times, T, I_rho);
    Itot = sum(Transmission_poly, 2).*deltaLambda;
end


end
