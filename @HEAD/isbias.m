function [IsBias,Res]=isbias(Sim,varargin)
% Check if HEAD/SIM object is a bias image
% Package: @HEAD
% Description: Check if HEAD/SIM objects are bias images.
%              The program can look for bias images in a set of SIM or HEAD
%              objects, using header keyword or/and file name. It also
%              check if the images global std is consistent with the
%              readnoise, and if the images are similar (in difference or
%              ratio) to a template image.
% Input  : - An HEAD object or a SIM object. For HEAD objects can look
%            for bias images only based on header keywords.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'TypeKeyVal' - The value of the IMTYPE like keyword in
%                           the header in case of a bias image.
%                           Either a string or a cell array of strings
%                           (i.e., multiple options).
%                           Default is {'bias','Bias','BIAS'}.
%            'TypeKeyDic' - IMTYPE keyword names. If empty use the istype.m
%                           default which is:
%                           {'TYPE','IMTYPE','OBSTYPE','IMGTYP','IMGTYPE'}.
%                           Default is empty.
%            'FileNameStr'- Optional search for bias images based on file
%                           name. This is a substring that if contained
%                           in the file name then the image is a bias.
%                           If empty then do not use this option.
%                           Default is empty.
%            'RN'         - Readnoise of the detector. This is either
%                           a numeric value for the readnoise [e-], or
%                           a string or cell array of strings containing
%                           possible readnoise header keywords [e-].
%                           If empty then do not select bias images based
%                           on std.
%                           Default is {'READNOI','READNOIS','RON'}.
%            'DefRN'      - A default value for the readnoise in case it
%                           is not found in header. Default is 10 e-.
%            'Nrn'        - The image is declared as a possible bias image
%                           if its global std multiplied by the gain is
%                           smaller than the readnoise multiplied by this
%                           value. Default is 2.
%            'StdFun'     - Function to use for the calculation of the
%                           global std of a SIM object {@std | @rstd}.
%                           Default is @rstd (slower than @std).
%            'GAIN'       - The detector gain [e-/ADU]. Either a numeric
%                           value or a string or a cell array of strings
%                           of possible header keywords containing the 
%                           detector gain.
%            'DefGAIN'    - A default value for the gain in case it
%                           is not found in header. Default is 1.
%            'Template'   - A matrix or a SIM image containing a template
%                           image which will be compared with each input
%                           image.
%                           If empty then do not use template search.
%                           Default is empty.
%            'TemplateType'- The comparison with the template can either
%                           done by difference ('diff') or by ratio
%                           ('ratio'). Default is 'diff'.
%            'TemplateNoise'- The image is a possible bias image if the
%                           global std of the comparison with the template
%                           is smaller than this value (in the native units
%                           of the image). Default is 30.
%            'CombType'   - The function tha will be used to combine all
%                           the bias search criteria {@all|@any}.
%                           Default is @all (i.e., requires that all the
%                           criteria are fullfilled).
%                           However, only active searches are being
%                           combined. For example, if 'Template' is empty
%                           then its results (false) will not be combined.
%            'SelectMethod'- Method by which to select the best keyword
%                           value. See getkey_fromlist.m for details.
%                           Default is 'first'.
% Output : - A vector of logical flags indicating if each image is a
%            candidate bias image, based on the combined criteria.
%          - A structure array with additional information.
%            The following fields are available:
%            .IsBiasKey - IsBias based on IMTYPE header keyword
%            .IsBiasFN  - IsBias based on file name.
%            .IsBiasStd - IsBias based on std.
%            .IsBiasTempStd - IsBias based on comparison with template.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [IsBias,R]=isbias(S);
%          IsBias = isbias(S,'RN',[]);
% Reliable: 2
%--------------------------------------------------------------------------

ImageField         = 'Im';
ImageFileNameField = 'ImageFileName';


DefV.TypeKeyVal         = {'bias','Bias','BIAS'};
DefV.TypeKeyDic         = [];    % if empty use istype default.
DefV.FileNameStr        = [];    % e.g.., 'Bias' - if empty do not use file name
DefV.RN                 = {'READNOI','READNOIS','RON'};  % if empty do not use   % [e-]
DefV.DefRN              = 10; % [e-]
DefV.Nrn                = 2;
DefV.StdFun             = @rstd;
DefV.GAIN               = {'GAIN'};
DefV.DefGAIN            = 1;
DefV.Template           = [];    % either SIM or a matrix
DefV.TemplateType       = 'diff';
DefV.TemplateNoise      = 30;     % [ADU or ratio]
DefV.CombType           = @all;   % @all | @any
DefV.SelectMethod       = 'first';
if (numel(varargin)>0)
    %InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
    InPar = InArg.populate_keyval(DefV,varargin,mfilename);
else
    InPar = DefV;
end

if (~isempty(InPar.Template))
    if (isnumeric(InPar.Template))
        Template = SIM;
        Template.(ImageField) = InPar.Template;
    elseif (SIM.issim(InPar.Template))
        Template = InPar.Template;
    else
        error('Unknown Template format');
    end
end

% treat the input in case its an HEAD object
% Select bias images based on image TYPE keywords
IsBiasKey = istype(Sim,InPar.TypeKeyVal,InPar.TypeKeyDic);
IsBias    = IsBiasKey;

% treat the input in case its a SIM object
IsBiasFN      = false(numel(IsBias),1);
IsBiasStd     = false(numel(IsBias),1);
IsBiasTempStd = false(numel(IsBias),1);
if (SIM.issim(Sim))
    % select bias images based on file name
    if (~isempty(InPar.FileNameStr))
        % do not use file name
        IsBiasFN = ~Util.cell.isempty_cell(strfind({Sim.(ImageFileNameField)}.',InPar.FileNameStr));
        
        IsBias   = InPar.CombType([IsBias,IsBiasFN],2);
    end
    
    % select bias images based on std of image
    if (~isempty(InPar.RN))
        Std = InPar.StdFun(Sim);
        % get readnoise
        RN = cell2mat(getkey_fromlist(Sim,InPar.RN,InPar.SelectMethod));
        % change RN with NaN value to the default value
        RN(isnan(RN)) = InPar.DefRN;
        
        % get GAIN
        Gain = cell2mat(getkey_fromlist(Sim,InPar.GAIN,InPar.SelectMethod));
        % change RN with NaN value to the default value
        Gain(isnan(Gain)) = InPar.DefGAIN;
        
        IsBiasStd = (Std(:).*Gain)<(RN.*InPar.Nrn);
        
        IsBias   = InPar.CombType([IsBias,IsBiasStd],2);
    end
    
    
    % select images based on similarity to template
    if (~isempty(InPar.Template))
        switch lower(InPar.TemplateType)
            case 'diff'
                StdTempResid   = InPar.StdFun(Sim - Template);
            case 'ratio'
                StdTempResid   = InPar.StdFun(Sim./Template);
            otherwise
                error('Unknown TemplateType option');
        end
        IsBiasTempStd  = StdTempResid<InPar.TemplateNoise;
        
        IsBias   = InPar.CombType([IsBias,IsBiasTempStd],2); 
    end
        
end

Res.IsBiasKey     = IsBiasKey;
Res.IsBiasFN      = IsBiasFN;
Res.IsBiasStd     = IsBiasStd;
Res.IsBiasTempStd = IsBiasTempStd;





