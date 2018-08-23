function [Links,FitsName,Images,Headers]=get_dss(RA,Dec,FOV,Save,Filter)
% Get link to and the FITS image of the DSS
% Package: VO.POSS
% Description: Get link to and the FITS image of a digital sky survey
%              image (POSS-I/II, UKST survey). Read the FITS file into
%              matlab.
% Input  : - J2000 RA in [radians] or [H M S] or sexagesimal string.
%          - J2000 Dec in [radians] or [H M S] or sexagesimal string.
%          - Field of view [arcmin arcmin], default is [12 12].
%          - Save FITS images to disk {'gzip' | 'y' | 'n'},
%            default is 'gzip'.
%          - Filter/epoch to retrieve,
%            {'2'|'2R'|'2B'|'2I'|'DSS1'|'1R'|'1B'|'QV'|'HST-P2'},
% Output : - Cell array of links to each gif image.
%          - Cell array of fits images name.
%            (e.g., "DSS_2R_######.###+######.##.fits").
%            The fits images are retrieved only if the user ask for this
%            argument.
%          - Cell array of images matrix.
%            If cell is NaN - image can't retrieved.
%          - Cell array of images header.
%            Note that keywords are stored in:
%            Header{:,:}.PrimaryData.Keywords (see example).
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jul 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Link]=VO.POSS.get_dss(1,1,'n','2R');
% Reliable: 1
%---------------------------------------------------------------------------
DSS_Server = 'http://archive.stsci.edu/cgi-bin/dss_search';

RAD = 180./pi;

ColSep     = '%3A';   % : seperator
SignP      = '%2B';   % +
SignN      = '-';      % -
Equinox    = 'J2000';  

DefFOV    = [12 12];
DefSave   = 'gzip';
DefFilter = '2'; 

if (nargin==2)
   FOV          = DefFOV;
   Save         = DefSave;
   Filter       = DefFilter;
elseif (nargin==3)
   Save         = DefSave;
   Filter       = DefFilter;
elseif (nargin==4)
   Filter       = DefFilter;
elseif (nargin==5)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (length(FOV)==1)
   FOV = [FOV, FOV];
end

switch Filter
 case '2'
    Version  = 'poss2ukstu';
 case '2R'
    Version  = 'poss2ukstu_red';
 case '2B'
    Version  = 'poss2ukstu_blue';
 case '2I'
    Version  = 'poss2ukstu_ir';
 case 'DSS1'
    Version  = 'dss1';
 case '1R'
    Version  = 'poss1_red';
 case '1B'
    Version  = 'poss1_blue';
 case 'QV'
    Version  = 'quickv';
 case 'HST-P2'
    Version  = 'phase2';
 otherwise
    error('Unknown Filter Option');
end

RA  = celestial.coo.convertdms(RA,'gH','H');
Dec = celestial.coo.convertdms(Dec,'gD','D');


N = size(RA,1);


for I=1:1:N

   if (Dec(I,1)==+1)
      Sign = SignP;
   else
      Sign = SignN;
   end

   ImType = 'gif';
   URL = sprintf('%s?v=%s&r=%02d%s%02d%s%05.2f&d=%s%02d%s%02d%s%04.1f&e=%s&h=%5.2f&w=%5.2f&f=%s&c=none&fov=NONE&v3=',DSS_Server,Filter,RA(I,1),ColSep,RA(I,2),ColSep,RA(I,3),Sign,Dec(I,2),ColSep,Dec(I,3),ColSep,Dec(I,4),Equinox,FOV(1),FOV(2),ImType);

   ImType = 'fits';
   FitsURL = sprintf('%s?v=%s&r=%02d%s%02d%s%05.2f&d=%s%02d%s%02d%s%04.1f&e=%s&h=%5.2f&w=%5.2f&f=%s&c=none&fov=NONE&v3=',DSS_Server,Filter,RA(I,1),ColSep,RA(I,2),ColSep,RA(I,3),Sign,Dec(I,2),ColSep,Dec(I,3),ColSep,Dec(I,4),Equinox,FOV(1),FOV(2),ImType);

   Links{I} = URL;

   if (nargout>1)
      StrRA       = celestial.coo.convertdms(RA(I,:),'H','SHn');
      StrDec      = celestial.coo.convertdms(Dec(I,:),'D','SDn');

      FitsName{I} = sprintf('DSS_%s_%s%s.fits',Filter,StrRA,StrDec);
   
      RunStr    = sprintf('wget -q "%s" -O %s',FitsURL,FitsName{I});
      system(RunStr);
   end

   if (nargout>2)
      %--- read fits image ---
      Images{I}   = fitsread(FitsName{I});

      if (nargout>3)
         %--- read fits header ---
         Headers{I} = fitsinfo(FitsName{I});
      end

      switch Save
       case 'n'
          % delete fits image
          delete(FitsName{I});
       case 'y'
          % do nothing
       case 'gzip'
          % gzip fits image
          gzip(FitsName{I});
          % change FitsName accordingly
          FitsName{I} = sprintf('%s%s',FitsName{I},'.gz');
       otherwise
          error('Unknown Save Option');
      end

   end

end
