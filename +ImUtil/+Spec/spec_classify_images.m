function InfoS=spec_classify_images(Sim,varargin)
%--------------------------------------------------------------------------
% spec_classify_images function                                     ImSpec
% Description: Classify a list of images given their header keywords.
%              For each imagee the program retrieve or calculate the
%              value of some header keywords or some constants.
%              Moreover, for another set of keywords the retrival mode
%              may depends on the value of a specific keyword (e.g.,
%              specifying if a blue or red arm image).
% Input  : - Set of images.
%            The following inputs are possible:
%            (1) Cell array of image names in string format.
%            (2) String containing wild cards (see create_list.m for
%                option). E.g., 'lred00[15-28].fits' or 'lred001*.fits'.
%            (3) Structure array of images (SIM).
%                The image should be stored in the 'Im' field.
%                This may contains also mask image (in the 'Mask' field),
%                and an error image (in the 'ErrIm' field).
%            (4) Cell array of matrices.
%            (5) A file contains a list of image (e.g., '@list').
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'CheckImage' - Read the image {false|true}. Default is false.
%            'Key_*' - A general keyword to retrieve which value will be
%                      populated. For example if the keword is 'Key_JD'
%                      then the outout structure will contain a 'JD'
%                      field. The value of 'Key_JD' can beone of the
%                      following:
%                      (1) A scalar which will be copied to the output
%                          structure 'JD' field.
%                      (2) A string. In this case the value will be
%                          retrieved from the header keyword with this
%                          name (e.g., 'OBSJD').
%                      (3) A function handle, which return the value.
%                      (4) A cell array in which the first element is
%                          a function handle, which return the value.
%                          The second element is an header keyword name.
%                          The value of this header keyword will be 
%                          provided to the function as the first input
%                          argument (e.g., {@julday, 'UTSHUT'}).
%                          Additional elements will be send to the function
%                          as additional input arguments
%                          (e.g., {@convertdms,'RA','gH','d'}).
%            'Key_ArmIndex' - A special keyword that specify how to get
%                          the index of the arm in the 'ArmList' input.
%                          (e.g., {@Util.string.find_strcmpi,'FPA',DefV.ArmList}).
%            'ArmList' - A cella array of possible arm names.
%            'Arm_*'   - A cell array which number of elements is
%                        identical to the number of possible arms in
%                        'ArmList'. Each one will be interpreted as
%                        'Key_*' for the specific arm.
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
%            images2sim.m
% Output : - A structure array will all the requested keyword values.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: InfoS=spec_classify_images('*.fits');  % use default values
%          % get default values from instrument definition structure
%          InstCell=struct2keyvalcell(InstrumentDef(2));
%          InfoS=spec_classify_images('*.fits',InstCell{:});
% Reliable: 2
%--------------------------------------------------------------------------


%ImageField  = 'Im';
%HeaderField = 'Header';
FileField   = 'ImageFileName';
%MaskField   = 'Mask';
%BackImField = 'BackIm';
%ErrImField  = 'ErrIm';


% input arguments
DefV.CheckImage        = false;
% Arm definitions
DefV.ArmList           = {'DBSP_BLUE','DBSP_RED2'};
% keyword header which value depends on arm
DefV.Arm_DispAxis      = {@()'y'     ,@()'x'};  % e.g., functions that return strings (if string get from header...)
DefV.Arm_FlipDisp      = {true       ,false};   % flip dispersion direction so increasing wavelength
DefV.Arm_PixScale      = {0.389      ,0.293};
DefV.Arm_DispScale     = {1.54       ,0.86};
DefV.Arm_SatLevel      = {65000      ,65000};
DefV.Arm_CollapseRange = {[1001 1500],[501 1000]};
DefV.Arm_DispRange     = {[200 2600] ,[420 3720]};
DefV.Arm_SpatRange     = {[60 370]   ,[20 330]};
DefV.Arm_NrangeWC      = {4          ,2};       % number of ranges in wavelengtth calibration
DefV.Arm_WaveCalibType = {@(){'arc'} ,@(){'arc','sky'}};  % will perform both sky and arc wave calib % the first one is the primary
% spatial keyword
DefV.Key_ArmIndex      = {@Util.string.find_strcmpi,'FPA',DefV.ArmList};
% keyword header which value is arm independendt
DefV.Key_Long          = -116.8649;
DefV.Key_Lat           = 33.3563;
DefV.Key_Date          = 'UTSHUT';
DefV.Key_JD            = {@julday, 'UTSHUT'};
DefV.Key_RA            = {@convertdms,'RA','gH','d'};
DefV.Key_Dec           = {@convertdms,'DEC','gD','d'};
DefV.Key_Equinox       = 2000.0;
DefV.Key_Object        = 'OBJECT';
DefV.Key_ImType        = 'IMGTYPE';
DefV.Key_Gain          = 'GAIN';
DefV.Key_RN            = 'RON';
DefV.Key_Dichroic      = {@Util.string.spacedel, 'DICHROIC'};
DefV.Key_SlitWidth     = {@Util.string.spacedel, 'APERTURE'};
DefV.Key_Turret        = {@Util.string.spacedel, 'TURRET'};
DefV.Key_Lamps         = {@Util.string.spacedel, 'LAMPS'};
DefV.Key_Arcs          = {@spec_dbsp_lamps, 'LAMPS'};
DefV.Key_ExpTime       = 'EXPTIME';
DefV.Key_BinX          = {@spec_dbsp_bin,'CCDSUM','x'};
DefV.Key_BinY          = {@spec_dbsp_bin,'CCDSUM','y'};
DefV.Key_TrimSec       = {@ccdsec_convert,'TSEC1'};
DefV.Key_BiasSec       = {@ccdsec_convert,'BSEC1'};
DefV.Key_Detector      = 'DETNAM';
DefV.Key_ArmName       = 'FPA';
DefV.Key_AM            = {@str2num_nan,'AIRMASS'};
DefV.Key_LST           = {@convertdms,'OBSLST','gH','f'};
DefV.Key_SlitPA        = {@str2num_nan,'CASSPA'};
DefV.Key_ParAng        = {@str2num_nan,'PARALLAC'};
DefV.Key_Grating       = 'GRATING';
DefV.Key_GratingAng    = 'ANGLE';

InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

%--- dill with Arm dependent keywords ---
FN       = fieldnames(InPar);
PossKey  = regexp(FN,'Arm_','split');
Nf       = cellfun(@numel,PossKey);
ArmKeys  = FN(Nf==2).';
Narmkey  = numel(ArmKeys);

%--- Identify all keys to retrieve from header
PossKey  = regexp(FN,'Key_','split');
Nf       = cellfun(@numel,PossKey);
Keys     = FN(Nf==2).';
Nkey     = numel(Keys);


%--------------------------------------
%--- Read the headers of all images ---
%--------------------------------------

Sim = images2sim(Sim,varargin{:},'ReadHead',true,'ReadImage',InPar.CheckImage);
Nim = numel(Sim);

for Iim=1:1:Nim,
    % image name
    if (isfield(Sim(Iim),FileField)),
        InfoS(Iim).(FileField) = Sim(Iim).(FileField);
    else
        InfoS(Iim).(FileField) = [];
    end
    
    % arm independent values
    for Ikey=1:1:Nkey,
        % check the content of the input argument
        KeyToPass = InPar.(Keys{Ikey});
        Val       = spec_classify_images_populate(KeyToPass,Sim(Iim));

        InfoS(Iim).(Keys{Ikey}(5:end)) = Val;
    end
    
    % arm dependent values
    if (isfield(InfoS(Iim),'ArmIndex')),
        for Iarmkey=1:1:Narmkey,
            if (~isnan(InfoS(Iim).ArmIndex)),
                KeyToPass = InPar.(ArmKeys{Iarmkey}){InfoS(Iim).ArmIndex};
            else
                KeyToPass = [];
            end
            Val       = spec_classify_images_populate(KeyToPass,Sim(Iim));

            InfoS(Iim).(ArmKeys{Iarmkey}(5:end)) = Val;
        end
    end
    
end


%-----------------------------------------------------------
function Val=spec_classify_images_populate(KeyToPass,SimIim)
%-----------------------------------------------------------

        if (isnumeric(KeyToPass) || islogical(KeyToPass)),
            if (isempty(KeyToPass)),
                % ignore
                Val = [];
            else
               % numeric keyword value - populate with value
               Val = KeyToPass;
            end
        elseif (ischar(KeyToPass)),
            % char keyword - get keyword value from header as is
            ValKey = sim_getkeyval(SimIim,KeyToPass,'ConvNum',false);
            Val    = ValKey{1};
        elseif (isa(KeyToPass,'function_handle')),
            % get value from a function
            Val = feval(KeyToPass);
        elseif (iscell(KeyToPass)),
            % get vale from a function which first which first argument is
            % the avlue of the keyword
            if (isempty(KeyToPass{2})),
                % call function without header keyword
                Val    = feval(KeyToPass{3:end});
            else
                % call function with header keyword
                ValKey = sim_getkeyval(SimIim,KeyToPass{2},'ConvNum',false);
                if (isnan(ValKey{1})),
                    Val = NaN;
                else
                    Val    = feval(KeyToPass{1},ValKey{1}, KeyToPass{3:end});
                end
            end
        else
            error('Unknwon Key type');
        end
