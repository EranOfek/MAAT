function wget_all_skymapper(varargin)
% SHORT DESCRIPTION HERE
% Package: VO
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Apr 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: VO.prep.wget_all_skymapper
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

CatField = AstCat.CatField;

DefV.CatNameBase          = 'SkyMapper';
DefV.HTM_Level            = 8;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Col = {'raj2000','dej2000','e_raj2000','e_dej2000','mean_epoch','rms_epoch',...
            'flags','class_star','radius_petro','a','b','pa',...
            'u_psf','e_u_psf','u_petro',...
            'v_psf','e_v_psf','v_petro',...
            'g_psf','e_g_psf','g_petro',...
            'r_psf','e_r_psf','r_petro',...
            'i_psf','e_i_psf','i_petro',...
            'z_psf','e_z_psf','z_petro'};
        
ColCell = {'RA','Dec','ErrRA','ErrDec','MeanEpoch','MeanEpochRMS',...
           'Flags','ClassStar','radPetro','A','B','PA',...
           'uPSF','uErrPSF','uPetro',...
           'vPSF','vErrPSF','vPetro',...
           'gPSF','gErrPSF','gPetro',...
           'rPSF','rErrPSF','rPetro',...
           'iPSF','iErrPSF','iPetro',...
           'zPSF','zErrPSF','zPetro'};
ColUnits = {'rad','rad','mas','mas','MJD','day','','','arcsec','arcsec','arcsec','deg',...
            'mag','mag','mag',...
            'mag','mag','mag',...
            'mag','mag','mag',...
            'mag','mag','mag',...
            'mag','mag','mag',...
            'mag','mag','mag'};


[HTM,LevList] = celestial.htm.htm_build(InPar.HTM_Level);

Ptr  = LevList(end).ptr;
Side = LevList(end).side;
Nhtm = numel(Ptr);

for Ihtm=1:1:Nhtm
    Ihtm
    IndHTM = Ptr(Ihtm);
    RA  = mean(HTM(Ptr(Ihtm)).coo(:,1));
    Dec = mean(HTM(Ptr(Ihtm)).coo(:,2));
    
    if (Dec<(15./RAD))
        pause(10);
        
        tic;
        [Out,Link] = VO.SkyMapper.skymapper_cat_search(RA,Dec,'Radius',Side.*RAD.*3600);
        toc
        
        Nsrc = size(Out.(CatField),1);
        if (Nsrc>0)
        
            SelOut = Out.col_select(Col);
            SelOut.Cat = table2array(SelOut.(CatField));
            SelOut.Cat(:,1:2) = SelOut.(CatField)(:,1:2)./RAD;

            Flag = celestial.htm.in_polysphere(SelOut.Cat(:,1:2),HTM(Ptr(Ihtm)).coo);

            if (sum(Flag)>0)
            
                [FileName,DataName] = catsHTM.get_file_var_from_htmid(InPar.CatNameBase,IndHTM,100);

                catsHTM.save_cat(FileName,DataName,SelOut.(CatField)(Flag,:),2,30);
            end
        end
    end
end
