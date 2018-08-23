function [Sim,Center,FlagValid]=stamp_coo(Sim,RA,Dec,HalfSizeStampXY,varargin)
% Get stamp from SIM object array by WCS coordinates.
% Package: @SIM
% Description: Get stamp from SIM object array (i.e., trim images) by
%              RA/Dec (or WCS) coordinates.
%              This function supports two options: regions of the stamp
%              outside the image are padded by some value (e.g., NaN),
%              or the requested coordinates are shifted such that it will
%              not contain regions outside the image.
%              This function also convert the WCS to the new stamp.
%              For trimming by X/Y pixel coordinates see SIM/trim_image.m.
%              The stamps are always rounded to whole pixels values.
%              For sub pixel shift is needed use stamp_xy.m instead.
% Input  : - A SIM object of images.
%          - A single J2000.0 R.A., [radians, sexagesimal, or [H M S]].
%            See convertdms.m for options.
%          - A single J2000.0 Dec., [radians, sexagesimal, or [Sign D M S]].
%            See convertdms.m for options.
%          - Stamp size half size [X, Y] size. If scalar then X and Y size
%            are identical.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'StampMethod' - The method by which to extract the stamp:
%                            'pad'  - Requested coordinates are located
%                                     at the stamp center (i.e.,
%                                     HalfSizeStampXY+1). Regions outside
%                                     the valid images are padded.
%                            'valid'- Select only images which requested
%                                     coordinate is contained within the
%                                     MinDist from image boundries. Than
%                                     Shift the requested coordinate such
%                                     that the stamps will contains valid
%                                     regions of the images. The final
%                                     stamp do not contain NaNs.
%                            Default is 'valid'.
%            'SectionPad'  - Value for padding. Default is NaN.
%            'MinDist'     - Minimum distance from image boundry. In the
%                            StampMethod='valid' option, if the requested
%                            coordinate is less than MinDist from image
%                            edge, then this image is decalred as invalid
%                            and a stamp will not be created for this
%                            image. Default is 64 pix.
%            'TrimPar'     - Cell array of additional arguments to pass to
%                            SIM/trim_image.m.
%                            E.g., {'ExecField',{'Im'}}
% Output : - A SIM object containing the image stamps.
%            For the StampMethod='valid' the number of stamps images is not
%            necesserly equal the number of input images as some of the
%            images may be declared as invalid.
%          - The [X,Y] pixel coordinates of the requested RA/Dec
%            coordinates in the stamp. For StampMethod='pad', this is
%            always the stamp center. For the 'valid' option', this is not
%            necesserly the case.
%          - A vector of logical flags indicating if the images are valid
%            (i.e., a stamp was created).
% See also: SIM/trim_image.m, trim_image.m, SIM/stamp_xy.m.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [S2,C,F]=stamp_coo(S1,'09:54:45.0','+68:57:00.0',512);
% Reliable: 2
%--------------------------------------------------------------------------

DefV.StampMethod        = 'valid'; %'pad'|'valid'
DefV.SectionPad         = NaN;
DefV.MinDist            = 64;
DefV.TrimPar            = {};
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (numel(HalfSizeStampXY)==1)
    HalfSizeStampXY = [HalfSizeStampXY, HalfSizeStampXY];
end

% go over all images and get the X/Y coordinates corresponding to the
% RA/Dec positions:
Nsim     = numel(Sim);
X        = zeros(Nsim,1);
Y        = zeros(Nsim,1);
for Isim=1:1:Nsim
    [X(Isim),Y(Isim)]  = coo2xy(Sim(Isim),RA,Dec);
end
RoundX = round(X);
RoundY = round(Y);

% get the image size
SizeIm = imagesize(Sim);  % size [X,Y]

FlagValid = true(Nsim,1);
RetEmpty  = false;
Coo = [RoundX, RoundY, HalfSizeStampXY(1).*ones(Nsim,1), HalfSizeStampXY(2).*ones(Nsim,1)];
switch lower(InPar.StampMethod)
    case 'pad'
        % trim and pad regions out of boundries by pad value:
        Center = repmat((HalfSizeStampXY+1),Nsim,1);
        SectionMethod = 'center';
    case 'valid'
        % Select only images which requiested coordinate is within image
        % and some minimum distance from bounderies:
        
        % Search for images which requested coordinate is within the image
        % and at least MinDist pixels from image edges:
        FlagValid = RoundX>InPar.MinDist & RoundX<(SizeIm(:,1)-InPar.MinDist) & ...
                    RoundY>InPar.MinDist & RoundY<(SizeIm(:,2)-InPar.MinDist);
        
        if (all(~FlagValid))
            warning('The requested central location is outside the image+MinDist boundry - return empty');
            RetEmpty = true;
        else
            % The logic here is suppose to keep the stamp size the same as
            % requested byut to shift the requested position such that the
            % stamp bounderies will contain the images.
            Xmin      = min(RoundX(FlagValid) - HalfSizeStampXY(1));
            Xmax      = max(RoundX(FlagValid) + HalfSizeStampXY(1));
            SubX = min(Xmin,0);
            AddX = min(min(SizeIm(:,1))-Xmax,0);
            Coo(:,1) = Coo(:,1) - SubX + AddX;
            CenterX = HalfSizeStampXY(1) - AddX + SubX;

            Ymin      = min(RoundY(FlagValid) - HalfSizeStampXY(2));
            Ymax      = max(RoundY(FlagValid) + HalfSizeStampXY(2));
            SubY = min(Ymin,0);
            AddY = min(min(SizeIm(:,2))-Ymax,0);
            Coo(:,2) = Coo(:,2) - SubY + AddY;
            CenterY = HalfSizeStampXY(2) - AddY + SubY;

            Center = [CenterX, CenterY];
            SectionMethod = 'center';
        end
    otherwise
        error('Unknown StampMethod option');
end

if (RetEmpty)
    Sim    = SIM(Nsim,1);
    Center = nan(Nsim,2);
else
    Sim = trim_image(Sim(FlagValid),Coo,'SectionMethod',SectionMethod,'SectionPad',InPar.SectionPad,InPar.TrimPar{:});
end   
