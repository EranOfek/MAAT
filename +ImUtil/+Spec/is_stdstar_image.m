function [IsStd,ImageStdData,NumStd,IsBright]=is_stdstar_image(Sim,varargin)
%--------------------------------------------------------------------------
% is_stdstar_image function                                         ImSpec
% Description: Given a list of FITS images or SIM, look for standard
%              stars spectra. The search is done by comparing the RA/Dec
%              or/and object name keywords with the database of
%              spectroscopic standard stars.
% Input  : - List of images to check for std stars images.
%            The following inputs are possible:
%            (1) Cell array of image names in string format.
%            (2) String containing wild cards (see create_list.m for
%                option). E.g., 'lred00[15-28].fits' or 'lred001*.fits'.
%            (3) Structure array of images (SIM).
%                The image should be stored in the 'Im' field.
%                This may contains also mask image (in the 'Mask' field),
%                and an error image (in the 'ErrIm' field).
%            (4) A file contains a list of image (e.g., '@list').
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'KeyRA'   - Image header R.A. keyword. Default is 'RA'.
%                        Keyword value should be in deg or sexagesimal.
%            'KeyDec'  - Image header Dec. keyword. Default is 'Dec'.
%                        Keyword value should be in deg or sexagesimal.
%            'Equinox' - Equinox of image header coordinates.
%                        This can be a string containing an image
%                        header keyword of the Equinox, or a scalar
%                        indicating the Julian year of the equinox.
%                        Default is 2000.
%            'KeyObj'  - Image header object name keyword.
%                        Default is 'OBJECT'.
%            'CheckCoo'- Look for standard stars using the image
%                        coordinates {true|false}.
%                        Default is true.
%            'CheckName'- Look for standard stars using the image
%                        object name {true|false}.
%                        Default is false.
%            'Radius'   - Search radius [arcsec]. Default is 60.
%            'NameSearch'- Name search method:
%                          'exact' - to perform an exact name search
%                                    (default).
%                          'substr'- to search for substring within the
%                                    names.
%            'CheckImage'- Check image for the exsistence of a bright
%                         trace with high S/N {true|false}.
%                         Default is true.
%            'DispDir'   - Dispersion direction {'x'|'y'}. Default is 'x'.
%            'MinDN'     - Minimal counts required for the peak of the
%                          median collapsed std star observation.
%                          Default is 7000.
%            'FieldName' - Field name containing the header in the
%                        structure returned by fitsinfo.m.
%                        Default is 'PrimaryData'. If empty then use
%                        default. If NaN then will attempt to look for 
%                        the correct field.
%            'Ind'     - Index of image header, in header structure.
%                        Default is 1.
%            'Verbose' - Print progress messages {true|false}.
%                        Default is false.
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
%            images2sim.m
% Output : - A flag vector (IsStd) indicating if the image is a standard
%            star image or not.
%          - A cell vector in which the standard star name is listed.
%            If empty cell then the corresponding image is not a
%            standard star image.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [IsStd,ImageStdData,NumStd,IsBright]=is_stdstar_image(Sim);
% Reliable: 2
%--------------------------------------------------------------------------
import Util.string.*

RAD = 180./pi;


ImageField  = 'Im';
HeaderField = 'Header';
%FileField   = 'ImageFileName';
%MaskField   = 'Mask';
%BackImField = 'BackIm';
%ErrImField  = 'ErrIm';

% input parameters
DefV.KeyRA         = 'RA';
DefV.KeyDec        = 'Dec';
DefV.Equinox       = 2000.0;
DefV.KeyObj        = 'OBJECT';
DefV.CheckCoo      = true;
DefV.CheckName     = false;
DefV.Radius        = 60;
DefV.NameSearch    = 'exact';
DefV.CheckImage    = true;
DefV.DispDir       = 'x';
DefV.MinDN         = 1000;
DefV.FieldName     = 'PrimaryData';
DefV.Ind           = 1;
DefV.Verbose       = false;

InPar = set_varargin_keyval(DefV,'n','use',varargin{:});


if (isstruct(Sim)),
    % do nothing
    InputSim = true;
    Nim      = numel(Sim);
else
    [~,ImageListCell] = create_list(Sim,NaN);
    InputSim = false;
    Nim      = numel(ImageListCell);
end


% go over images
IsStd        = false(Nim,1);
IsBright     = false(Nim,1);
NumStd       = zeros(Nim,1);
ImageStdData = cell(Nim,1);

for Iim=1:1:Nim,
    % read image
    if (InputSim);
        Header = Sim(Iim).(HeaderField);
        if (InPar.CheckImage),
            Image = Sim(Iim).(ImageField);
        end
    else
        % read from FITS
        Header = fits_header_cell(ImageListCell{Iim},InPar.FieldName,InPar.Ind);
        if (InPar.CheckImage),
            Image = fitsread(ImageListCell{Iim},InPar.FitsPars{:});
        end
    end
    
    % search by coordinates
    if (InPar.CheckCoo),
        % Read RA
        NewCellHead = cell_fitshead_getkey(Header,{InPar.KeyRA},'NaN');
        RA          = NewCellHead{1,2};
        if (isnan(RA)),
            error('Image header (of image %d) RA keyword is not present',Iim);
        end
        if (isnumeric(RA)),
            % assume already in deg
        elseif (ischar(RA)),
            RA = spacedel(RA);
            if (isempty(strfind(RA,':'))),
                % deg in string
                RA = str2double(RA);
            else
                % sexagesimal
                RA = convertdms(RA,'SH','d');  % deg
            end
        else
            error('Image header (of image %d) RA keyword is not valid',Iim);
        end
        
        % Read Dec
        NewCellHead = cell_fitshead_getkey(Header,{InPar.KeyDec},'NaN');
        Dec          = NewCellHead{1,2};
        if (isnan(Dec)),
            error('Image header (of image %d) Dec keyword is not present',Iim);
        end
        if (isnumeric(Dec)),
            % assume already in deg
        elseif (ischar(Dec)),
            Dec = spacedel(Dec);
            if (isempty(strfind(Dec,':'))),
                % deg in string
                Dec = str2double(Dec);
            else
                % sexagesimal
                Dec = convertdms(Dec,'SD','d');  % deg
            end
        else
            error('Image header (of image %d) Dec keyword is not valid',Iim);
        end
        
        % Equinox
        if (ischar(InPar.Equinox)),
            % read Equinox from image header
            NewCellHead = cell_fitshead_getkey(Header,{InPar.Equinox},'NaN');
            Equinox          = str2num_nan(NewCellHead{1,2});
            if (isnan(Equinox)),
                error('Image header (of image %d) Equinox keyword is not present',Iim);
            end
        else
            Equinox = InPar.Equinox;
        end
        
        % precess to J2000.0
        if (Equinox~=2000),
            Coo = coco([RA, Dec]./RAD, sprintf('j%6.1f',Equinox),'j2000');
            RA  = Coo(:,1).*RAD;  % deg
            Dec = Coo(:,2).*RAD;  % deg
        end
        
        Res = search_specphot_stand(RA./RAD,Dec./RAD,InPar.Radius);
        
        NumStd(Iim) = numel(Res);  % number of Std stars found
        if (NumStd(Iim)>0),
            IsStd(Iim) = true;
            ImageStdData{Iim} = Res;
        end
    end
    
    % search by name
    if (InPar.CheckName), 
            
        % Read object name
        NewCellHead = cell_fitshead_getkey(Header,{InPar.KeyObj},'NaN');
        ObjectName  = NewCellHead{1,2};
        if (isnan(ObjectName)),
            error('Image header (of image %d) %s keyword is not present',Iim,InPar.KeyObj);
        end
        ObjectName = spacedel(ObjectName);
         
        % search by object name
        Res = search_specphot_stand(ObjectName,InPar.NameSearch);
        
        NumStd(Iim) = numel(Res);  % number of Std stars found
        if (NumStd(Iim)>0),
            IsStd(Iim) = true;
            ImageStdData{Iim} = Res;
        end
    end
    
    % check image
    if (InPar.CheckImage),
        Collapse=spec_collapse_dispaxis(Image,'DispDir',InPar.DispDir);
        if (max(Collapse{1})>InPar.MinDN),
            IsBright(Iim) = true;
        end
    end
    
end

if (InPar.Verbose),
    fprintf('Standard star images search include total of %d images\n',Nim);
    if (InPar.CheckImage),
        fprintf('Found %d good standard star images out of %d candidate images\n',length(find(NumStd==1 & IsBright)),length(find(IsStd)));
    else
        fprintf('Found %d std star images with\n',length(find(NumStd==1)));
    end
end
    
