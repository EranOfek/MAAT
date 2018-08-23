function Res=sedm_bias_subtraction(varargin)
%--------------------------------------------------------------------------
% sedm_bias_subtraction function                                      SEDM
% Description: Subtract bias from SEDM image/s.
% Input  : * Arbitrary number of ...,key,val,... input arguments.
%            The following keywords are available:
%            'ListAll' - List of all images available (or only bias
%                        images) from which to select
%                        the bias images. The bias images will be selected
%                        automatically using the header information.
%                        Note that IFU and RC images will be reduced
%                        seperatly. See create_list.m for allowed input.
%                        Default is '*.fits'.
%            'ListRed' - List of images from which to subtract bias. This
%                        can be a subset or non-overlapping set of the
%                        first input argument.
%                        If empty then use the first input argument.
%                        Default is empty.
%            'BaseIFU' - Image base name for the IFU images.
%                        Default is 'ifu'.
%            'BaseRC'  - Image base name for the RC images.
%                        Default is 'rc'.
%            'BiasKeyword' - The value of the OBJECT header keyword for
%                        bias images. Default is 'bias'.
%            'IFU_BiasImageName' - IFU output bias image name.
%                        Default is 'BiasIFU.fits'.
%            'RC_BiasImageName' - RC output bias image name.
%                        Default is 'BiasRC.fits'.
%            'IFU_BiasStdImageName' - IFU output bias std image name.
%                        Default is 'BiasStdIFU.fits'.
%            'RC_BiasStdImageName' - RC output bias std image name.
%                        Default is 'BiasStdRC.fits'.
%            'Prefix'  - The prefix of the bias subtracted images.
%                        Default is 'b_'.
%            'BiasAlgo'- Bias subtraction algorithm {'new'|'old'}.
%                        Default is 'new'.
% Output : - Structure containing the statistics of the output bias images.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Sep 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: StatBias=sedm_bias_subtraction;
% Reliable: 2
%--------------------------------------------------------------------------
import Util.cell.*

DefV.ListAll              = '*.fits';
DefV.ListRed              = [];
DefV.BaseIFU              = 'ifu';
DefV.BaseRC               = 'rc';
DefV.BiasKeyword          = 'Calib: bias'; %'bias';
DefV.IFU_BiasImageName    = 'BiasIFU.fits';
DefV.IFU_BiasStDImageName = 'BiasStdIFU.fits';
DefV.RC_BiasImageName     = 'BiasRC.fits';
DefV.RC_BiasStDImageName  = 'BiasStdRC.fits';
DefV.Prefix               = 'b_';
DefV.BiasAlgo             = 'new';
DefV.IFU_OverScanSection  = [1 2048 2045 2047];
DefV.RC_OverScanSection   = []; %[1 2048 2045 2047];
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

if (isempty(InPar.ListRed)),
   InPar.ListRed = InPar.ListAll;
end

[~,ListCell.Red] = create_list(InPar.ListRed,NaN);
[~,ListCell.All] = create_list(InPar.ListAll,NaN);

ListCell.RedIFU = ListCell.Red(find(~isempty_cell(strfind(ListCell.Red,InPar.BaseIFU))));
ListCell.AllIFU = ListCell.All(find(~isempty_cell(strfind(ListCell.All,InPar.BaseIFU))));
ListCell.RedRC  = ListCell.Red(find(~isempty_cell(strfind(ListCell.Red,InPar.BaseRC))));
ListCell.AllRC  = ListCell.All(find(~isempty_cell(strfind(ListCell.All,InPar.BaseRC))));


% select bias images from All
HeadKeys = {'OBJECT'};
%[~,KeywordS] = mget_fits_keyword(ListCell.AllIFU,HeadKeys);
[~,~,KeywordS]=fits_mget_keys(ListCell.AllIFU,HeadKeys);
if (~isempty(KeywordS)),
    ListCell.BiasIFU = ListCell.AllIFU(strcmpi({KeywordS.OBJECT},InPar.BiasKeyword));
else
    ListCell.BiasIFU = [];
end
   
[~,KeywordS] = mget_fits_keyword(ListCell.AllRC,HeadKeys);
if (~isempty(KeywordS)),
    ListCell.BiasRC = ListCell.AllRC(strcmpi({KeywordS.OBJECT},InPar.BiasKeyword));
else
    ListCell.BiasRC = [];
end


switch lower(InPar.BiasAlgo)
    case 'new'
        % new algorithm, including overscan bias subtraction

        % read IFU images
        I = 1;
        Res(I).ImageType    = 'IFU';
        if (~isempty(ListCell.AllIFU)),
            Sim = images2sim(ListCell.AllIFU);
            % subtract bias
            [Sim,BiasSim]=sim_bias(Sim,'BiasList',ListCell.BiasIFU);
            Res(I).BiasName    = sims2fits(BiasSim,'OutName',InPar.IFU_BiasImageName);
            Res(I).BiasStdName = sims2fits(BiasSim,'OutName',InPar.IFU_BiasStDImageName);
            % subtract overscan bias
            if (~isempty(InPar.IFU_OverScanSection)),
                Sim=sim_suboverscan(Sim,'BiasSec',InPar.IFU_OverScanSection);
            end
            % read bias subtracted images as FITS files
            Res(I).OutImageName = sims2fits(Sim,'OutName',ListCell.AllIFU,'Prefix',InPar.Prefix);
        end
        
        I = 2;
        Res(I).ImageType    = 'RC';
        if (~isempty(ListCell.AllRC)),
            Sim = images2sim(ListCell.AllRC);
            % subtract bias
            [Sim,BiasSim]=sim_bias(Sim,'BiasList',ListCell.BiasRC);
            Res(I).BiasName    = sims2fits(BiasSim,'OutName',InPar.RC_BiasImageName);
            Res(I).BiasStdName = sims2fits(BiasSim,'OutName',InPar.RC_BiasStDImageName);
            % subtract overscan bias
            if (~isempty(InPar.RC_OverScanSection)),
                Sim=sim_suboverscan(Sim,'BiasSec',InPar.RC_OverScanSection);
            end
            % read bias subtracted images as FITS files
            Res(I).OutImageName = sims2fits(Sim,'OutName',ListCell.AllIFU,'Prefix',InPar.Prefix);
        end
        
        
        
    case 'old'
        % construct and subtract bias
        if (~isempty(ListCell.BiasIFU)),
           bias_fits(ListCell.AllIFU,ListCell.BiasIFU,[],InPar.IFU_BiasImageName,...
                     'OutPrefix',InPar.Prefix,...
                     'StD',InPar.IFU_BiasStDImageName);
           Res.IFU = imstat_fits(InPar.IFU_BiasImageName);
        end

        if (~isempty(ListCell.BiasRC)),
           bias_fits(ListCell.AllRC,ListCell.BiasRC,[],InPar.RC_BiasImageName,...
                     'OutPrefix',InPar.Prefix,...
                     'StD',InPar.RC_BiasStDImageName);

           Res.RC  = imstat_fits(InPar.RC_BiasImageName);
        end

    otherwise
        error('Unknown BiasAlgo option');
end

