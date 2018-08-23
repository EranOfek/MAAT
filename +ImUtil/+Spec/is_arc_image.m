function [IsArc,ImageArcName]=is_arc_image(Sim,varargin)
%--------------------------------------------------------------------------
% is_arc_image function                                            ImBasic
% Description: Given a list of FITS images or SIM, look for arc
%              (wavelength calibration) images. The search is done by
%              looking for specific keyword values in the image headers.
% Input  : - List of images to check for bias images.
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
%            'Instrum' - Instrument name. See def_instrument_arc.m
%                        for options. Default is 'P200-DBSP'.
%            'Compare' - How to compare the keyword value with the
%                        keyword values map:
%                        'strcmpi' - string comparison using strcmpi
%                        'bit'     - assumes that the keyword value is
%                                    a string of bits (e.g., '00101'),
%                                    and if one or more of this bits is
%                                    open then the image is an arc image.
%                                    This may return multiple arc names.
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
% Output : - A flag vector (IsArc) indicating if the image is an arc
%            image or not.
%          - A cell vector in which the standard arc name is listed.
%            If empty cell then the corresponding image is not an
%            arc image.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [IsArc,ArcName]=is_arc_image('blue*.fits');
%          [IsArc,ArcName]=is_arc_image('blue*.fits','Compare','bit');
% Reliable: 2
%--------------------------------------------------------------------------



ImageField  = 'Im';
HeaderField = 'Header';
%FileField   = 'ImageFileName';
%MaskField   = 'Mask';
%BackImField = 'BackIm';
%ErrImField  = 'ErrIm';

% input parameters
% read header
DefV.FieldName     = [];
DefV.Ind           = 1;
% arc identification
DefV.Instrum       = 'P200-DBSP';
DefV.Compare       = 'strcmpi';  % {'strcmpi'|'bit'}
% check images
DefV.CheckImage    = false;  % FFU - not supported
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

[ArcKey,ArcKeyValue,ArcName]=def_instrument_arc(InPar.Instrum);

% go over images
IsArc       = false(Nim,1);
ImageArcName= cell(Nim,1);

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
    
    % check candidate image header
    NewCellHead = cell_fitshead_getkey(Header,{ArcKey},'NaN');
    Val         = NewCellHead{1,2};
    if (isnan(Val)),
        % ARC keyword is not available
        error('ARC image header keyword is not available');
    end
    
    % truncate spaces
    Val = Util.string.spacedel(Val);
    switch lower(InPar.Compare)
        case 'strcmpi'
            FoundFlag     = strcmpi(Val,ArcKeyValue);
            if (any(FoundFlag)),
               ImageArcName{Iim} = ArcName{FoundFlag};
               IsArc(Iim)    = true;
            end
        case 'bit'
            %Val           = bin2dec(Val);
            BitStatus     = find(fliplr(bitget(bin2dec(Val),(1:length(Val)))));
            if (~isempty(BitStatus)),
                ImageArcName{Iim} = ArcName(BitStatus);
                IsArc(Iim)        = true;
            end
            
        otherwise
            error('Unknown Compare option');
    end
    
end


if (InPar.Verbose),
    fprintf('Arc images search include total of %d images\n',Nim);
    if (InPar.CheckImage),
        %fprintf('Found %d good bias images out of %d candidate bias images\n',length(find(IsBias & IsGoodNoise & IsGoodMean)),length(find(IsBias)));
    else
        fprintf('Found %d candidate arc images\n',length(find(IsArc)));
    end
end
    
