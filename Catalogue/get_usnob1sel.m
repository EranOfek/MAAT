function Cat=get_usnob1sel(RA,Dec,FOV,Shape)
%---------------------------------------------------------------------------
% get_usnob1sel function     
% Description: Get USNO-B1.0 catalog from local disk (Using Stephen Levine
%              C program: ubcone).
% Input  : - RA (J2000.0),    [HH MM SS] or [Radians]
%          - Dec (J2000.0), [Sign, DD MM SS] or [Radians]
%          - Width [RA Width, Dec height] or diameter of the Field of view [arcsec].
%          - Field of view shape:
%            'b' - box, default.
%            'c' - circle.
% Output : - USNO-B1.0 catalog sorted by Declination, columns:
%            1  - RA J2000.0 [radians]
%            2  - Dec J2000.0 [radians]
%            3  - Error in RA [mas]
%            4  - Error in Dec [mas]
%            5  - Epoch
%            6  - RA PM[mas/yr]
%            7  - Dec PM [mas/yr]
%            8  - PM probability [0 no PM .. 9 most probable]
%            9  - Error in RA PM [mas/yr]
%            10 - Error in Dec PM [mas/yr]
%            11 - The overall error in the fit in the RA components [mas]
%            12 - The overall error in the fit in the Dec components [mas]
%            13 - Number of detections
%            14 - Flag : 1 = indicates this object has been correlated
%                            with a known proper motion catalogue star
%                            (catalogues include LHS, NLTT, Giclas etc.)
%                        2 = indicates proximity to a diffraction spike
%                            or a bright star (ie a Tycho2 star)
%                        4 = indicates that this star is a YS4.0 star.
%                            (YS4.0 was the intermediate catalogue used
%                            to link the Schmidt plate detections to the
%                            Tycho2 astrometric reference frame.)
%                        5 = the star is near a known proper motion
%                            catalogue star, and is also in YS4.0.
%                        8 = indicates that this star is a Tycho2 star
%                            that has been added into the catalogue. 
%            15 - B1 magnitude
%            16 - Photometric calib source for B1 magnitude
%            17 - Survey ID, where first digit is:
%                 0 = POSS-I O 1 = POSS-I E 2 = POSS-II J 3 = POSS-II F
%                 4 = SERC-J or SERC-EJ 5 = ESO-R or SERC-ER 6 = AAO-R
%                 7 = POSS-II N 8 = SERC-I 9 = POSS-N
%                 and last three digits are the field ID.
%            18 - Star/galaxy classification:
%                 The range is from 0 to 11. 0 means quite dissimilar,
%                 11 means very similar.
%                 NaN means that there was no value computed. 
%            19 - The B1 plate, Xi residual in the tangent plane between
%                 the measured position on this plate, and the position
%                 computed from the best fit.This is given as a value
%                 in decimal arcseconds between -50.00 and +49.99 arcseconds.
%            20 - The B1 plate, Eta residual in the tangent plane between
%                 the measured position on this plate, and the position
%                 computed from the best fit. This is given as a value in
%                 decimal arcseconds between -50.00 and +49.99 arcseconds. 
%            21 - The B1 plate, look back index into the PMM raw scan files.
%            22..28 - like 15..21, but for the R1 plate
%            29..35 - like 15..21, but for the B2 plate
%            36..42 - like 15..21, but for the R2 plate
%            43..49 - like 15..21, but for the I2 plate
% Comments: /UBC/ubcone + USNO-B1.0 should be installed
% Tested : MATLAB 5.3
%     By : Eran O. Ofek                October 2003
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------
RAD = 180./pi;
USNO_PATH = '/home/castor3/usnob';

if (nargin==3),
   Shape = 'b';
elseif (nargin==4),
   % do nothing
else
  error('Illegal number of input arguments');
end

if (length(RA)==3),
   RA = convertdms(RA,'H','r');
end
if (length(Dec)==4),
   Dec = convertdms(Dec,'D','R');
end

% convert radians to hours/degrees
RA  = RA.*RAD./15;
Dec = Dec.*RAD;
if (RA>=24),
   RA = RA - 24;
end

ExtName  = '_usnob.asc';
TmpFile  = 'tmp';
FileName = [TmpFile,ExtName];
eval(['!rm -rf ',FileName,'*']);
switch Shape
 case 'b'
    if (length(FOV)==1),
       FOV = [FOV(1), FOV(1)];
    end
    RunStr1 = ['!setenv USNOBDIR ',USNO_PATH, '; ubcone -F R -O ',TmpFile,' -P ',sprintf('%10.6f',RA),' -p ',sprintf('%9.5f',Dec),' -S ',sprintf('%10.6f',FOV(1)./3600),' -s ',sprintf('%10.6f',FOV(2)./3600),' -Z L -z D'];
 case 'c'
    RunStr1 = ['!setenv USNOBDIR ',USNO_PATH, '; ubcone -F C -O ',TmpFile,' -P ',sprintf('%10.6f',RA),' -p ',sprintf('%9.5f',Dec),' -S ',sprintf('%10.6f',FOV./3600),' -Z L -z D'];
 otherwise
    error('Unknown Shape option');
end
RunStr1
eval(RunStr1);

switch Shape
 case 'b'
    TailStr = ['!tail +27 ',FileName, ' >> ',FileName,'.1'];
 case 'c'
    TailStr = ['!tail +26 ',FileName, ' >> ',FileName,'.1'];
 otherwise
    error('Unknown Shape option');
end
eval(TailStr);


Cat      = load([FileName,'.1']);
if (~isempty(Cat)),
   Cat      = Cat(:,2:50);
   Cat(:,1) = Cat(:,1).*15./RAD;
   Cat(:,2) = Cat(:,2)./RAD;
   % replace s/g class 19->NaN
   Col      = 18;
   I        = find(Cat(:,Col)==19);
   Cat(I,Col)   = NaN;
   Col      = 25;
   I        = find(Cat(:,Col)==19);
   Cat(I,Col)   = NaN;
   Col      = 32;
   I        = find(Cat(:,Col)==19);
   Cat(I,Col)   = NaN;
   Col      = 39;
   I        = find(Cat(:,Col)==19);
   Cat(I,Col)   = NaN;
   Col      = 46;
   I        = find(Cat(:,Col)==19);
   Cat(I,Col)   = NaN;
   % replace 0 mag with NaN
   Col      = 15;
   I        = find(Cat(:,Col)==0);
   Cat(I,Col) = NaN;
   Col      = 22;
   I        = find(Cat(:,Col)==0);
   Cat(I,Col) = NaN;
   Col      = 29;
   I        = find(Cat(:,Col)==0);
   Cat(I,Col) = NaN;
   Col      = 36;
   I        = find(Cat(:,Col)==0);
   Cat(I,Col) = NaN;
   Col      = 43;
   I        = find(Cat(:,Col)==0);
   Cat(I,Col) = NaN;


end
