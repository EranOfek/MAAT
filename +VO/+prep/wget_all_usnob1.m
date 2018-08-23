function wget_all_usnob1(varargin)
% Retrieve USNO-B1 catalog from VizieR and format into HDF5/HTM (catsHTM)
% Package: VO.prep
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: VO.prep.wget_all_usnob1
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

DefV.Nsrc                 = [];
DefV.LevelHTM             = 9;
DefV.NfilesInHDF          = 100;
DefV.FileBase             = 'USNOB1';
DefV.ColDec               = 2;
DefV.StepInd              = 30;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

[HTM,LevList]=celestial.htm.htm_build(InPar.LevelHTM);

Nh = numel(LevList(end).ptr);
for Ih=1:1:Nh
  
    
    IndHTM  = LevList(end).ptr(Ih);
    
    if (isempty( find(IndHTM==InPar.Nsrc(:,1)) ))
        
        tic;
        
        MeanRA  = mean(HTM(IndHTM).coo(:,1));
        MeanDec = mean(HTM(IndHTM).coo(:,2));

        D = celestial.coo.sphere_dist_fast(MeanRA,MeanDec,HTM(IndHTM).coo(:,1),HTM(IndHTM).coo(:,2));
        Radius = max(D);

        % search
        pause(2);
        try
            Cat = VO.VizieR.cds_astcat_search('findusnob1',MeanRA,MeanDec,'OutType','astcat','RadiusUnits','rad','Radius',Radius,'CooUnits','rad');
            if (size(Cat.Cat,1)==0)
                pause(60);
                Cat = VO.VizieR.cds_astcat_search('findusnob1',MeanRA,MeanDec,'OutType','astcat','RadiusUnits','rad','Radius',Radius,'CooUnits','rad');
            end
            % convert coo to radians
            Cat.Cat(:,1:2) = Cat.Cat(:,1:2)./RAD;
            % select in HTM
            Flag = celestial.htm.in_polysphere(Cat.Cat(:,1:2),HTM(IndHTM).coo,2);
            Cat.Cat = Cat.Cat(Flag,:);
            Cat.Cat = sortrows(Cat.Cat,InPar.ColDec);
            Ncat    = size(Cat.Cat,1);

        catch
            % failed
            load Failed.mat
            Failed = [Failed; Ih];
            save Failed.mat Failed

            Ncat = 0;
        end


        if (Ncat>0)
            [FileName,DataName]=catsHTM.get_file_var_from_htmid(InPar.FileBase,IndHTM,InPar.NfilesInHDF);

            catsHTM.save_cat(FileName,DataName,Cat.Cat,InPar.ColDec,InPar.StepInd);
        end
        [Ih, Nh, IndHTM, Ncat, toc]
    end
end

    
    