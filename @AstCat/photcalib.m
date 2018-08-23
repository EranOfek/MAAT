function [Sim,ZP]=photcalib(Sim,varargin)
% Calculate and apply photometric calibration.
% Package: @AstCat
% Description: Calculate and apply photometric calibration, against
%              a reference catalog, to an AstCat object array
%              (or a SIM images with catalogs).
%              The catalog should contains all the necesseary columns.
%              The user can specify the fitted formula.
% Input  : - SIM array with catalogs.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'PreDef'  - Use predefined set of parameters:
%                        Check out the inner function pre_defined_par
%                        for parameter definitions.
%                        Options include 'r_sdss_gr' (i.e., solve r-band
%                        observations using SDSS r and g-r), 'g_sdss_gr',
%                        'u_sdss_ug',... 'r_apass_gr',...
%                        You can overide any of the parameters by set
%                        them individually.
%            'Xcol'    - Column name, in the catalog, containing the sources
%                        X position. Default is 'XWIN_IMAGE'.
%            'Ycol'    - Column name, in the catalog, containing the sources
%                        Y position. Default is 'YWIN_IMAGE'.
%            'Flagcol' - Column name, in the catalog, containing flags.
%                        Default is 'FLAGS'.
%            'FlagVal' - Value in FLAGS indicating good object for
%                        photometric calibration.
%                        Default is 0 (i.e., all the bits are false).
%            'ColSelectRange' - Three column cell array of of column name,
%                        lower range and upper range. Only sources which
%                        column value is in this range will participate in
%                        the fit. Default is to use 'PreDef'.
%            'EquationMag' - The "observable" of the set of equations to
%                        solve. This is a string containing some functions
%                        on existing columns. Default is to use 'PreDef'.
%            'EquationZP' - A cell array of string containing the functions
%                        defining the design matrix section.
%                        E.g., {'APASS_htm_g - PASS_htm_r','AIRMASS',...
%                               'mod(XWIN_IMAGE,1)','mod(YWIN_IMAGE,1)',
%                               'XWIN_IMAGE.^2'}.
%                        Default is to use 'PreDef'.
%            'EquationApply' - A vector of logicals indicating which columns
%                        in 'EquationZP' should be used when applying the
%                        ZP to the instrumental magnitude.
%                        Default is to use 'PreDef'.
%                        These defaults mean that the color term will not
%                        be applied (this is should be done after the
%                        colors are known).
%            'EquationZPerr' - Equation defining the error in the
%                        observables.
%                        Default is to use 'PreDef'.
%            'BlockSize' - Indicating if solution should be done on
%                        entire image or several regions separatly
%                        (i.e., use this if the zero point may be position
%                        depandent).
%                        Default is [] (i.e., solve the entire image).
%                        Other examples are [1024 1024] will break the
%                        image to 1024x1024 pixels sub images,
%                        while [-2 -2] will divide the image to 2 by 2 sub
%                        images.
%            'BlockBuffer' - Extend sub images sideways to overlap with
%                        nearby sub images [pix]. Default is 32.
%            'ApplyZP' - A cell array of column names to which to apply
%                        the ZP. If empty then do not apply the ZP.
%                        Note that the ZP will be applied only for sources
%                        which are found within the image bounderies.
%                        Default is {'MAG_APER'}.
% Output : - SIM array of images in which the zero points are applied to
%            the magnitude columnds.
%          - Structure array containing the zero point solutions for each
%            image.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% Sim=images2sim('PTF_201511075122_i_p_scie_t121737_u026000155_f02_p100072_c10.fits');
% Sim=mextractor(Sim);
% Sim=simcat_addcol(Sim);
% Sim=xcat(Sim,'ExtCats','APASS_htm');
% [Sim1,ZP]=photcalib(Sim);
% Reliable: 2
%--------------------------------------------------------------------------
import Util.array.*

%ImageField       = 'Im';
CatField         = 'Cat';
CatColField      = 'Col';
%CatColCellField  = 'ColCell';


DefV.PreDef             = 'r_apass_gr';
% look for pre-defined parameters
Ipd = find(strcmpi(varargin(1:2:end-1),'predef'));
if (~isempty(Ipd))
    PreDef = varargin{Ipd+1};
    DefV   = pre_defined_par(PreDef);
else
    PreDef = DefV.PreDef;
    DefV   = pre_defined_par(PreDef);
end
DefV.Xcol               = 'XWIN_IMAGE';
DefV.Ycol               = 'YWIN_IMAGE';
DefV.Flagcol            = 'FLAGS';
DefV.FlagVal            = 0;
%DefV.ColSelectRange     = {'APASS_htm_r',14,16};
%{'wget_sdss_modelMag_r',15,17.5; 'wget_sdss_type',5.5,6.5}; %
%DefV.EquationMag        = 'APASS_htm_r - MAG_APER';
%'wget_sdss_modelMag_r - MAG_APER';% 
%DefV.EquationZP         = {'XWIN_IMAGE-1024','YWIN_IMAGE-2048',...
%                           '(XWIN_IMAGE-1024).^2','(YWIN_IMAGE-2048).^2'...
%                           '(XWIN_IMAGE-1024).*(YWIN_IMAGE-2048)',...
%DefV.EquationZP         = {'APASS_htm_g - APASS_htm_r'}; %,...
%{'wget_sdss_modelMag_g - wget_sdss_modelMag_r'};
                          % 'mod(XWIN_IMAGE,1)',...
                          % 'mod(YWIN_IMAGE,1)'};
%DefV.EquationApply      = [false]; %, true, true];
                           % 'MAG_APER',...
                           %'AIRMASS',...
                           %'mod(XWIN_IMAGE,1)',...
                           %'mod(YWIN_IMAGE,1)'}; 
%DefV.EquationZPerr      = {'sqrt(APASS_htm_rerr.^2 + APASS_htm_gerr.^2 + MAGERR_APER.^2)'};
%{'sqrt(wget_sdss_modelMagErr_g.^2 + wget_sdss_modelMagErr_r.^2 + MAG_APER.^2)'};
DefV.BlockSize          = []; %[1024 1024];  % if empty use entire image
DefV.BlockBuffer        = 32;
DefV.ApplyZP            = {}; %{'MAG_APER'};

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            
Nmag = numel(InPar.ApplyZP);

% Read images into SIM array
Sim  = images2sim(Sim,varargin{:});
Nsim = numel(Sim);

%ResY   = simcat_colarith(Sim,InPar.EquationMag);
%ResH   = simcat_colarith(Sim,InPar.EquationZP);
%ResErr = simcat_colarith(Sim,InPar.EquationZPerr);
ResY   = col_arith(Sim,InPar.EquationMag,'mat');
ResH   = col_arith(Sim,InPar.EquationZP,'mat');
ResErr = col_arith(Sim,InPar.EquationZPerr,'mat');


%SelCol = simcat_colrange(Sim,InPar.ColSelectRange,@all);
SelCol = col_range(Sim,InPar.ColSelectRange,@all);


Npar = numel(InPar.EquationZP) + 1;  % the +1 is for the ZP par that is added
InPar.EquationApply = [true, InPar.EquationApply];
% for each SIM
ZP     = struct('Par',cell(Nsim,1),'ParErr',cell(Nsim,1),...
                'RMS',cell(Nsim,1),'Chi2',cell(Nsim,1),...
                'Npar',cell(Nsim,1),'Ndof',cell(Nsim,1),...
                'ListEdge',cell(Nsim,1),'ListCenter',cell(Nsim,1));
for Isim=1:1:Nsim
    %  assume that SIM contains catalog with all the necessery columns
    if (SIM.issim(Sim(Isim)))
        ImSize=imagesize(Sim,true);
    else
        ImSize=naxis(Sim);
    end
    
    % calculate the edges and centers of blocks in which to solve for ZP
    [ListEdge,ListCenter] = ImUtil.Im.image_blocks(ImSize,InPar.BlockSize,InPar.BlockBuffer,'simple');
    ListEdgeNoBuffer      = ImUtil.Im.image_blocks(ImSize,InPar.BlockSize,0,'simple');
    Nblock = size(ListEdge,1);   % number of blocks
    
    % select good stars
    %Flag = Sim(Isim).(CatField)(:,Sim(Isim).(CatColField).(InPar.Flagcol))==InPar.FlagVal & SelCol(Isim).Flag;
    Flag = maskflag_check(Sim(Isim).(CatField)(:,Sim(Isim).(CatColField).(InPar.Flagcol)),InPar.FlagVal) & SelCol(Isim).Cat(:,1);
    
        
    
    ZP(Isim).Par    = zeros(Nblock,Npar);
    ZP(Isim).ParErr = zeros(Nblock,Npar);
    ZP(Isim).RMS    = zeros(Nblock,1);
    ZP(Isim).Chi2   = zeros(Nblock,1);
    ZP(Isim).Npar   = zeros(Nblock,1);
    ZP(Isim).Ndof   = zeros(Nblock,1);
    for Iblock=1:1:Nblock
        % for each block
        
        % select stars in block & Good
        Iobj = find(Sim(Isim).(CatField)(:,Sim(Isim).(CatColField).(InPar.Xcol))>=ListEdge(Iblock,1) & ...
                    Sim(Isim).(CatField)(:,Sim(Isim).(CatColField).(InPar.Xcol))<=ListEdge(Iblock,2) & ...
                    Sim(Isim).(CatField)(:,Sim(Isim).(CatColField).(InPar.Ycol))>=ListEdge(Iblock,3) & ...
                    Sim(Isim).(CatField)(:,Sim(Isim).(CatColField).(InPar.Ycol))<=ListEdge(Iblock,4) & ...
                    Flag);
        Nobj = numel(Iobj);
        % construct design matrix
        H   = [ones(Nobj,1), ResH(Isim).Cat(Iobj,:)];
        Y   = ResY(Isim).Cat(Iobj);
        Err = ResErr(Isim).Cat(Iobj);
        
        %Tmp = Sim(Isim).(CatField)(Iobj,Sim(Isim).(CatColField).APASS_htm_r);

        Fnn = ~isnan(sum(H,2)) & ~isnan(Y) & ~isnan(Err); % & Tmp>14 & Tmp<16;
               
        
        [Par,~]      = lscov(H(Fnn,:),Y(Fnn),1./Err(Fnn).^2);
        Resid        = Y(Fnn) - H(Fnn,:)*Par;
        
        [~,~,Fresid] = clip_resid(Resid,'Method','StdP','Clip',[2 2]);
        %hist(Resid(Fresid))
        [Par,ParErr] = lscov(H(Fnn & Fresid,:),Y(Fnn & Fresid),1./Err(Fnn & Fresid).^2);
        Resid        = Y(Fnn & Fresid) - H(Fnn & Fresid,:)*Par;
        
        %Tmp = Sim(Isim).(CatField)(:,Sim(Isim).(CatColField).wget_sdss_modelMag_r);
        
        %plot(Tmp(Iobj(Fnn)),Resid,'.')
        %jj
        
        ZP(Isim).ListEdge   = ListEdge;
        ZP(Isim).ListCenter = ListCenter;
        
        ZP(Isim).Par(Iblock,:)    = Par.';
        ZP(Isim).ParErr(Iblock,:) = ParErr.';
        ZP(Isim).RMS(Iblock)      = std(Resid);
        ZP(Isim).Chi2(Iblock)     = sum(Resid.^2./Err(Fnn & Fresid).^2);
        ZP(Isim).Npar(Iblock)     = Npar;
        ZP(Isim).Ndof(Iblock)     = numel(Resid) - Npar;
        %ZP(Isim).AM(Iblock) = nanmedian(Sim(Isim).Cat(:,14));
      
        if (~isempty(InPar.ApplyZP))
            % Apply ZP to specific columns
            Iobj = find(Sim(Isim).(CatField)(:,Sim(Isim).(CatColField).(InPar.Xcol))>=ListEdgeNoBuffer(Iblock,1) & ...
                        Sim(Isim).(CatField)(:,Sim(Isim).(CatColField).(InPar.Xcol))<ListEdgeNoBuffer(Iblock,2) & ...
                        Sim(Isim).(CatField)(:,Sim(Isim).(CatColField).(InPar.Ycol))>=ListEdgeNoBuffer(Iblock,3) & ...
                        Sim(Isim).(CatField)(:,Sim(Isim).(CatColField).(InPar.Ycol))<ListEdgeNoBuffer(Iblock,4) );
            Nobj = numel(Iobj);
            H   = [ones(Nobj,1), ResH(Isim).Cat(Iobj,:)];
            H   = H(:,InPar.EquationApply);
            MagZP   = H*Par(InPar.EquationApply);
            
            % apply the ZP to the magnitude columns
            for Imag=1:1:Nmag
                Sim(Isim).(CatField)(Iobj,Sim(Isim).(CatColField).(InPar.ApplyZP{Imag})) = MagZP + Sim(Isim).(CatField)(Iobj, Sim(Isim).(CatColField).(InPar.ApplyZP{Imag}));
            end
          
        end
        
    end
    
end
    



function DefV=pre_defined_par(PreDef)

% set pre-defined solutions
if (~isempty(PreDef))
    switch lower(PreDef)
        case 'r_sdss_gr'
            % SDSS
            % r-band observations vs. SDSS r and g-r
            DefV.ColSelectRange     = {'wget_sdss_modelMag_r',15,17.5; 'wget_sdss_type',5.5,6.5};
            DefV.EquationMag        = 'wget_sdss_modelMag_r - MAG_APER';
            DefV.EquationZP         = {'wget_sdss_modelMag_g - wget_sdss_modelMag_r'};
            DefV.EquationApply      = (false);
            DefV.EquationZPerr      = {'sqrt(wget_sdss_modelMagErr_g.^2 + wget_sdss_modelMagErr_r.^2 + MAG_APER.^2)'};
        case 'g_sdss_gr'
            % SDSS
            % g-band observations vs. SDSS g and g-r
            DefV.ColSelectRange     = {'wget_sdss_modelMag_g',15,17.5; 'wget_sdss_type',5.5,6.5};
            DefV.EquationMag        = 'wget_sdss_modelMag_g - MAG_APER';
            DefV.EquationZP         = {'wget_sdss_modelMag_g - wget_sdss_modelMag_r'};
            DefV.EquationApply      = (false);
            DefV.EquationZPerr      = {'sqrt(wget_sdss_modelMagErr_g.^2 + wget_sdss_modelMagErr_r.^2 + MAG_APER.^2)'};
        case 'i_sdss_ri'
            % SDSS
            % i-band observations vs. SDSS i and r-i
            DefV.ColSelectRange     = {'wget_sdss_modelMag_i',15,17.5; 'wget_sdss_type',5.5,6.5};
            DefV.EquationMag        = 'wget_sdss_modelMag_i - MAG_APER';
            DefV.EquationZP         = {'wget_sdss_modelMag_r - wget_sdss_modelMag_i'};
            DefV.EquationApply      = (false);
            DefV.EquationZPerr      = {'sqrt(wget_sdss_modelMagErr_i.^2 + wget_sdss_modelMagErr_r.^2 + MAG_APER.^2)'};
        case 'u_sdss_ug'
            % SDSS
            % u-band observations vs. SDSS u and u-g
            DefV.ColSelectRange     = {'wget_sdss_modelMag_u',14,17.5; 'wget_sdss_type',5.5,6.5};
            DefV.EquationMag        = 'wget_sdss_modelMag_u - MAG_APER';
            DefV.EquationZP         = {'wget_sdss_modelMag_u - wget_sdss_modelMag_g'};
            DefV.EquationApply      = (false);
            DefV.EquationZPerr      = {'sqrt(wget_sdss_modelMagErr_u.^2 + wget_sdss_modelMagErr_g.^2 + MAG_APER.^2)'};
        case 'z_sdss_iz'
            % SDSS
            % z-band observations vs. SDSS z and i-z
            DefV.ColSelectRange     = {'wget_sdss_modelMag_z',14,17.5; 'wget_sdss_type',5.5,6.5};
            DefV.EquationMag        = 'wget_sdss_modelMag_z - MAG_APER';
            DefV.EquationZP         = {'wget_sdss_modelMag_i - wget_sdss_modelMag_z'};
            DefV.EquationApply      = (false);
            DefV.EquationZPerr      = {'sqrt(wget_sdss_modelMagErr_z.^2 + wget_sdss_modelMagErr_i.^2 + MAG_APER.^2)'};
        case 'r_apass_gr'
            % APASS
            % r-band observations vs. APASS r and g-r
            DefV.ColSelectRange     = {'APASS_htm_r',14,17.5};
            DefV.EquationMag        = 'APASS_htm_r - MAG_APER';
            DefV.EquationZP         = {'APASS_htm_g - APASS_htm_r'};
            DefV.EquationApply      = (false);
            DefV.EquationZPerr      = {'sqrt(APASS_htm_g.^2 + APASS_htm_r.^2 + MAG_APER.^2)'};
        case 'g_apass_gr'
            % APASS
            % g-band observations vs. APASS g and g-r
            DefV.ColSelectRange     = {'APASS_htm_g',14,17.5};
            DefV.EquationMag        = 'APASS_htm_g - MAG_APER';
            DefV.EquationZP         = {'APASS_htm_g - APASS_htm_r'};
            DefV.EquationApply      = (false);
            DefV.EquationZPerr      = {'sqrt(APASS_htm_g.^2 + APASS_htm_r.^2 + MAG_APER.^2)'};
        otherwise
            error('Unknown PreDef option');
    end
end
            