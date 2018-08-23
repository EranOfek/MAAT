function Spec=read_spec(FileName,FileType,varargin)
%------------------------------------------------------------------------------
% read_spec function                                                    ImSpec
% Description: Read a spectrum from a FITS file or ASCII file.
% Input  : - File name.
%          - File type. One of the following:
%            'text' - A text ascii file in which the first column is
%                     wavelength and the second is flux.
%            'fits' - FITS file (default).
%          * Arbitrary number of pairs of ...,keyword,value,...
%            The following keywords are available:
%            'FluxCor' - Multiply the output flux by this factor.
%                        Default is 1.
%            'WaveCor' - Multiply the output wavelength by this factor.
%                        Default is 1.
%            'ExtType' - FITS extension type (see fitsread.m).
%                        Default is 'primary'.
%            'Ext'     - FITS file extension number.
%                        Default is 1.
%            'Comment' - Comment style in ascii spectrum. Default is '#'.
%            'Header'  - Number of header lines to ignore in ascii spectrum.
%            'Keys'    - Cell array of additional header keywords to extract
%                        from FITS file.
%                        These keywords will be stored in the structure
%                        output parameter.
%                        Defaults {}.
% Output : - Structure containing the spectrum.
%            The fields depend on the content of the input file.
%            Default fields are:
%            .Wave  - Wavelength
%            .Flux  - Flux
%            Additional keywords
%            see read_sdss_spec.m foe SDSS keys.
%            or the keywords are specified by the COLUMN## header keywords.
% Tested : Matlab 7.8
%     By : Eran O. Ofek                    Dec 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Spec=read_fits_spec('ptf09uj.all.fits');
%------------------------------------------------------------------------------

Def.FileType = 'fits';
if (nargin==1),
   FileType = Def.FileType;
else
   % do nothing
end

DefV.FluxCor   = 1;
DefV.WaveCor   = 1;
DefV.ExtType   = 'primary';
DefV.Ext       = 1;
DefV.Comment   = '#';
DefV.Header    = 0;
DefV.Keys      = {};
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});


switch lower(FileType)
 case 'text'
    FID = fopen(FileName,'r');
    C   = textscan(FID,'%f %f %*[^\n]',...
                   'CommentStyle',InPar.Comment,...
                   'Headerlines',InPar.Header);
    fclose(FID);

    Spec.Wave = C{1}.*InPar.WaveCor;
    Spec.Flux = C{2}.*InPar.FluxCor;
 case 'fits'
    % read FITS header
    H = fitsinfo(FileName);

    % look for the TELESCOPE keyword
    [NewCellHead,Lines]=cell_fitshead_getkey(H.PrimaryData.Keywords,'TELESCOP');

    if (~isempty(findstr(NewCellHead{2},'SDSS'))),
       % SDSS observations
       S = read_sdss_spec(FileName);
       Spec = S{1};
    else
       % non-SDSS data
       KeywordVal = get_fits_keyword(FileName,{'CRVAL1','CDELT1','CRPIX1','CTYPE1','CD1_1'});
       S = fitsread(FileName,InPar.ExtType,InPar.Ext);
       if (sum(cellfun(@isnan,KeywordVal))>0),
           % one of the WCS keywords is missing
           % assume Wavelength is in the first column
           
           Spec.Wave = S(1,:).*InPar.WaveCor;
           Spec.Flux = S(2,:).*InPar.FluxCor;
           Spec.Pix  = [1:1:length(Spec.Wave)];

           % Check for extra column keywords
           % This extra columns will be identified using
           % an header keyword e.g., "COLUMN01" and kthe keyword value is the column name.
	   Icol = find(Util.cell.isempty_cell(strfind(lower(H.PrimaryData.Keywords(:,1)),'column'))==0);
           for I=1:1:length(Icol),
	      ColumnIndex = str2doubel(H.PrimaryData.Keywords{Icol(I),1}(7:8));
              ColumnName  = H.PrimaryData.Keywords{Icol(I),2};
              Spec.(ColumnName) = S(ColumnIndex,:);
           end
       else
	   if (isnan(CDELT1)==1),
	     CDELT1 = CD1_1;
           end
           switch lower(deblank(CTYPE1)),
            case 'linear'
                Wave = CRVAL1 + (VecPix - CRPIX1).*CDELT1;
                Spec.Wave = Wave.*InPar.WaveCor;
                Spec.Flux = S.*InPar.FluxCor;
                Spec.Pix  = VecPix;

             otherwise
                error('Unsported CTYPE1 value');
            end
       end
    end

    % get extra keywords from header
    [NewCellHead,Lines]=cell_fitshead_getkey(H.PrimaryData.Keywords,InPar.Keys);
    for I=1:1:size(NewCellHead,1),
       Spec.(NewCellHead{I,1}) = NewCellHead{I,2};
    end
 otherwise
    error('Unknown FileTyep option');
end




