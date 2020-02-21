function [URL,Spec,ColCell]=wget_sdss_spec(ID,varargin)
% wget SDSS FITS spectra and links
% Package: VO.SDSS
% Description: wget SDSS FITS spectra and links
% Input  : - A 3 column matrix of [Plate, MJD, Fiber] SDSS spectra ID
%            to retrieve.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'OutType' - Options: 'AstCat' | 'cellmat' | 'AstSpec'
%                        Default is 'AstCat'.
%            'Format' - Default is 'fits'.
%            'Survey' - 'boss' | ['eboss'] | 'sdss'
%            'DR'     - Default is 'dr16'.
%            'SpecType' - 'lite' | 'full'. Default is 'lite'.
%            'Run2D'    - Reduction version. Default is 'v5_10_0'.
%            'SaveFile' - Save FITS file to local disk.
%                         Default is true.
%            'UseLink'  - Which link to use for file retrieval.
%                         'link' - direct link (faster).
%                         'api' - API.
%            'WgetOption' - Retrieval method:
%                         'pwget' - Use www.pwget. Default.
%                         'websave' - Use websave.
%            'PwgetPar' - Additionl aparmeters to pass to www.pwget.
%            'NP'       - Number of parallel retrival in pwget.
% Output : - A structure array of URL/file names, with the following
%            fields:
%            'SpecView'  - Link to spectrum viewer.
%            'SpecAPI'   - Link to spectra retireval API.
%            'SpecLink'  - Direct link to spectra retrival.
%            'FileName'  - File name.
% License: GNU general public license version 3
% Reference: https://dr16.sdss.org/optical/spectrum/view/data/access
%     By : Eran O. Ofek                    Sep 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [URL,Spec,ColCell]=VO.SDSS.wget_sdss_spec([4055 55359 596])
% Reliable: 2
%--------------------------------------------------------------------------

DefV.OutType              = 'AstCat';
DefV.Format               = 'fits';
DefV.SpecType             = 'lite'; % 'lite | 'full'
DefV.Run2D                = 'v5_10_0';
DefV.DR                   = 'dr16';
DefV.Survey               = 'eboss'; % 
DefV.SaveFile             = true;
DefV.UseLink              = 'link';   % 'link' | 'api'
DefV.WgetOption           = 'pwget';  % 'pwget' | 'websave'
DefV.PwgetPar             = '';
DefV.Pwget.NP             = 10;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% instructions page:
%https://dr15.sdss.org/optical/spectrum/view/data/access

% optical spectrum viewer:
%https://dr15.sdss.org/optical/spectrum/view/?plateid=PLATEID&mjd=MJD&fiberid=FIBERID
% e.g.,  https://dr15.sdss.org/optical/spectrum/view/?plateid=4055&mjd=55359&fiberid=596
%        https://dr15.sdss.org/optical/spectrum/view?mjd=55359&fiberid=596&plateid=4055

% download using API:
%https://dr15.sdss.org/optical/spectrum/view/data/format=FORMAT/spec=SPEC?plateid=PLATEID&mjd=MJD&fiberid=FIBERID
% e.g.,  https://dr15.sdss.org/optical/spectrum/view/data/format=fits/spec=lite?plateid=4055&mjd=55359&fiberid=596

% direct download
%https://dr15.sdss.org/sas/dr15/SURVEY/spectro/redux/RUN2D/spectra/SPEC/PLATE4/spec-PLATE4-MJD-FIBERID4.fits
% e.g., https://dr15.sdss.org/sas/dr15/eboss/spectro/redux/v5_10_0/spectra/lite/4055/spec-4055-55359-0596.fits

% SURVEYShould be replaced by the appropriate survey (currently either 'sdss' for DR8, 'boss' for DR9 through DR12, or 'eboss' for DR13 and later)
% RUN2DShould be replaced by the reduction number.
% PLATE4Should be replaced by the zero-padded, 4-digit plate number.
% MJDShould be replaced by the MJD number.
% FIBERID4Should be replaced by the zero-padded, 4-digit fiber number.

PlateID = ID(:,1);
MJD     = ID(:,2);
FiberID = ID(:,3);

N = numel(PlateID);
for I=1:1:N
    % create viewr link:
    URL(I).SpecView = sprintf('https://%s.sdss.org/optical/spectrum/view?plateid=%d&mjd=%d&fiberid=%d',...
        InPar.DR,PlateID(I),MJD(I),FiberID(I));
    % create API link:
    URL(I).SpecAPI  = sprintf('https://%s.sdss.org/optical/spectrum/view/data/format=%s/spec=%s?plateid=%d&mjd=%d&fiberid=%d',...
        InPar.DR,InPar.Format,InPar.SpecType,PlateID(I),MJD(I),FiberID(I));
    % Create direct link:
    URL(I).SpecLink = sprintf('https://%s.sdss.org/sas/%s/%s/spectro/redux/%s/spectra/%s/%04d/spec-%04d-%d-%04d.fits',...
        InPar.DR,InPar.DR, InPar.Survey,  InPar.Run2D,InPar.SpecType,PlateID(I),PlateID(I),MJD(I),FiberID(I));

    URL(I).FileName = sprintf('spec-%04d-%d-%04d.fits',...
        PlateID(I),MJD(I),FiberID(I));
end

if (InPar.SaveFile)
    switch lower(InPar.UseLink)
        case 'link'
            Link = {URL.SpecLink};
        case 'api'
            Link = {URL.SpecAPI};
        otherwise
            error('Unknown UseLink option');
    end
    
    N = numel(Link);
    switch lower(InPar.WgetOption)
        case 'websave'
            for I=1:1:N
                websave(URL(I).FileName,Link{I});
            end
        case 'pwget'
            www.pwget(Link,InPar.PwgetPar,InPar.Pwget.NP);
        otherwise
            error('Unknown WgetOption option');
    end
end
    
% Read spectra to output argument
if (nargout>1)
    [Spec,ColCell] = FITS.read_sdss_spec({URL.FileName},'OutType',InPar.OutType);
end