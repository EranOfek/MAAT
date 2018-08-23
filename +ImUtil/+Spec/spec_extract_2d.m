function ExtractedSpec=spec_extract_2d(Image,Tr,varargin)
%--------------------------------------------------------------------------
% spec_extract_2d function                                          ImSpec
% Description: Given an image and a trace, cut a sub image along the
%              trace position in which the trace is aligned along the
%              X-axis.
% Input  : - A single image in one of the following forms:
%            (1) A structure array (SIM).
%            (2) A file name.
%            (3) A matrix.
%            (4) A cell array with a single file name string.
%            (5) A cell array with a single matrix image.
%          - Trace matrix [Disp_pos, Spat_pos] or the Trace structure
%            returned by spec_trace.m (will use [Trace.X, Trace.SmY]).
%            The trace structure may optionbal contain a bit Mask vector
%            in the 'Mask' field.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExtSemiW'   - Extraction region semi width (pixels).
%                           Default is 50.
%            'BinSpatSize'- Bin size in spatial direction (pixels).
%                           Default is 1.
%            'InterpMethod'- Interpolation method. Default is 'cubic'.
%            'ExMethod'   - One of the following extraction methods:
%                           'y'   - extract the spectrum in each dispersion
%                                   position, along the y axis.
%                           'vert'- extract the spectrum in each dispersion
%                                   position, along the direction
%                                   perpendicular to the trace (default).
%            'DispDir'    - Dispersion direction of image {'x'|'y'}.
%                           Default is 'x'. 'y' option may not work in
%                           some cases.
%            'Mask'       - An optional mask image. The mask image bits
%                           will be propogated into the extracted spectrum.
%                           Default is empty.
%                           Alternatively the mask can be provided in the
%                           structure image input ("Mask" field).
%                           This parameter overirde the content
%                           of the image Mask field.
%            'MaskTr'     - An optional trace mask vector.
%                           The trace mask vector can be provided either
%                           through the trace structure or through this
%                           parameter. This parameter will override
%                           the content of the trace structure mask.
%                           This mask corresponds to every dispersion
%                           position and it will be combined
%                           using an 'or' operator with
%                           the image mask. Default is empty.
%           'MaskType'   - Bit mask type to us, if created.
%                          Default is 'uint16'.
% Output : - A structure containing the following fields:
%            .Im      - Extracted sub image centered on trace.
%                      The width of this image is
%                      (2*ExtSemiW+1)*BinSpatSize.
%            .VecDisp - Vector of dispersion position.
%            .VecSpatR - Vector of spatial position relative to the trace
%                        spatial poistion.
%            .Mask    - An optional bit mask extracted in the position of
%                       the extracted image.
%            .BackIm  - If the input contains a background image,
%                       then the it will be extracted (around the trace).
%            .ErrIm  - If the input contains an error image,
%                       then the it will be extracted (around the trace).
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: ExtractedSpec=spec_extract_2d('lred0064.fits',Trace,'ExMethod','y');
% Reliable: 2
%--------------------------------------------------------------------------

TraceXField    = 'X';
TraceYField    = 'SmY';
TraceMaskField = 'Mask';


ImageField  = 'Im';
%HeaderField = 'Header';
%FileField   = 'ImageFileName';
MaskField   = 'Mask';
BackImField = 'BackIm';
ErrImField  = 'ErrIm';

DefV.ExtSemiW      = 50;
DefV.BinSpatSize   = 1;
DefV.InterpMethod  = 'cubic'; %'linear';
DefV.ExMethod      = 'vert';     % {'y' | 'vert'}
DefV.DispDir       = 'x';        % {'x' | 'y'}
DefV.Mask          = [];
DefV.MaskTr        = [];
DefV.MaskType      = 'uint16';
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

Sim=image2sim(Image);

% if (ischar(Image)),
%     Sim.(ImageField)  = fitsread(Image);
%     Sim.(HeaderField) = fitsinfo(Image);
%     Sim.(FileField)   = Image;
% elseif (isnumeric(Image)),
%     Sim.(ImageField)  = Image;
%     clear Image;
% elseif (isstruct(Image)),
%     Sim = Image;
%     clear Image;
% else
%     error('Unknown Image input type');
% end

if (isstruct(Tr)),
    Trace = [Tr.(TraceXField), Tr.(TraceYField)];
else
    Trace = Tr;
    Tr    = [];
end

switch lower(InPar.DispDir)
    case 'y'
        Sim.(ImageField) = [Sim.(ImageField)].';
    otherwise
        % do nothing
end

% MaskExist = false;
% if (isfield(Sim,MaskField)),
%     % .Mask field exist - use it.
%     MaskExist = true;
% else
%     if (~isempty(InPar.Mask)),
%         Sim.(MaskField) = InPar.Mask;
%         MaskExist       = true;
%     end
% end

SizeIm    = size(Sim.(ImageField));
SizeDisp  = SizeIm(2);    % dispersion axis size
SizeSpat  = SizeIm(1);    % spatial axis size
VecSpat   = (1:1:SizeSpat).';
VecDisp   = (1:1:SizeDisp).';


%----------------------------
%--- Extract the spectrum ---
%----------------------------
%Ntr                    = size(Trace,1);
ExtractedSpec.VecSpatR = (-InPar.ExtSemiW:InPar.BinSpatSize:InPar.ExtSemiW).';
Nspat                  = numel(ExtractedSpec.VecSpatR);   % number of extracted pixel
                                                         % in cross-dispersion direction
ExtractedSpec.VecDisp  = VecDisp;                                       

%ExtractedSpec.( ImageField) = zeros(Nspat,Ntr);

switch lower(InPar.ExMethod)
    case 'y'
        % extract in y-direction
        
        YI = bsxfun(@plus, Trace(:,2).', ExtractedSpec.VecSpatR );
        XI = repmat(Trace(:,1).',Nspat,1);
        
    case 'vert'
        % extract in direction perpendicular to the trace
        Alpha = atan(derivative(Trace(:,1),Trace(:,2)));
        
        YI = bsxfun(@plus, Trace(:,2).', ExtractedSpec.VecSpatR * cos(Alpha.') );
        XI = bsxfun(@minus, Trace(:,1).', ExtractedSpec.VecSpatR * sin(Alpha.') );
        
    otherwise
        error('Unknown ExMethod option');
end


% interpolate extracted region
ExtractedSpec.(ImageField) = interp2fast(VecDisp,VecSpat, Sim.(ImageField), XI, YI, InPar.InterpMethod);

% interpolate extrcated region of background image and error image
if (isfield(Sim,BackImField)),
    ExtractedSpec.(BackImField) = interp2fast(VecDisp,VecSpat, Sim.(BackImField), XI, YI, InPar.InterpMethod);
end
if (isfield(Sim,ErrImField)),
    ExtractedSpec.(ErrImField) = interp2fast(VecDisp,VecSpat, Sim.(ErrImField), XI, YI, InPar.InterpMethod);
end


%--- Image mask ---
% interpolate Mask field of original image
if (~isempty(InPar.Mask)),
    % if InPar.Mask exist then overide Sim.(MaskField)
    ExtractedSpec.(MaskField) = interp2fast(VecDisp,VecSpat, InPar.Mask, XI, YI, 'nearest');
else
    if (isfield(Sim,MaskField)),
        % interpolate
        ExtractedSpec.(MaskField) = interp2fast(VecDisp,VecSpat, Sim.(MaskField), XI, YI, 'nearest');
    else
        % init
        ExtractedSpec.(MaskField) = zeros(size(ExtractedSpec.(ImageField)),InPar.MaskType);
    end
end
switch lower(InPar.DispDir)
    case 'y'
        ExtractedSpec.(MaskField) = ExtractedSpec.(MaskField).';
    otherwise
        % do nothing
end


%--- Trace mask ---
if (~isempty(InPar.MaskTr)),
    % if InPar.MaskTr exist then overide Trace.(TraceMaskField)
    % data already in InPar.MaskTr
else
    if (~isempty(Tr)),
        if (isfield(Tr,TraceMaskField)),
            InPar.MaskTr = Tr.(TraceMaskField);
        else
            InPar.MaskTr = []; %zeros(SizeDisp,1,InPar.MaskType);
        end
    else
        InPar.MaskTr = []; %zeros(SizeDisp,1,InPar.MaskType);
    end
end

% combine ExtractedSpec.(MaskField) and InPar.MaskTr
if (~isempty(InPar.MaskTr)),
    InPar.MaskTr              = repmat(InPar.MaskTr.',Nspat,1);
    ExtractedSpec.(MaskField) = bitor(ExtractedSpec.(MaskField),InPar.MaskTr);
end


      

