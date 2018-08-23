function download_galex

RAD = 180./pi;

Query{1} = {'ra', 'dec', 'nuv_mag', 'nuv_magerr', 'fuv_mag', 'fuv_magerr', 'NUV_FWHM_WORLD', 'nuv_weight', 'fuv_weight'};  %'PhotoExtract.nobs_dat','PhotoExtract.fobs_dat'};
Query{2} = 'FROM photoobjall'; % left join PhotoExtract on photoobjall.photoExtractID=PhotoExtract.photoExtractID';
Query{3} = '';
%[Out,~,~,RS] = VO.MAST.query_casjobs_recur(Query,'boxcoo',[100 101 0 1].*pi./180,'StrRA','ra','StrDec','dec','Table','GALEX_GR6Plus7');

Level = 8;
[HTM, LevelHTM] = celestial.htm.htm_build(Level);
IndHTM = LevelHTM(Level).ptr;
Nhtm   = numel(IndHTM);

NfilesInHDF = 100;
FileBase    = 'GALEX';

%Update = zeros(0,4);
load Update.mat

for Ihtm=1:1:Nhtm
    pause(2);
    tic;
    Ih = IndHTM(Ihtm);
    
    MinRA  = min(HTM(Ih).coo(:,1));
    MaxRA  = max(HTM(Ih).coo(:,1));
    MinDec = min(HTM(Ih).coo(:,2));
    MaxDec = max(HTM(Ih).coo(:,2));
    
    try
        [Out,ColCell,Status,ResultOrig,QueryS] = VO.MAST.query_casjobs_recur(Query,'boxcoo',[MinRA MaxRA MinDec MaxDec],'StrRA','ra','StrDec','dec','Table','GALEX_GR6Plus7');
    
        % convert to radians
        Out.Cat(:,1:2) = Out.Cat(:,1:2)./RAD;

        % select in HTM
        Flag = celestial.htm.in_polysphere(Out.Cat(:,1:2),HTM(Ih).coo);
        Ns   = sum(Flag);

        Update = [Update; [Ihtm, Ih, Ns, Status]];

        if (Ns>0)
            [FileName,DataName]=HDF5.get_file_var_from_htmid(FileBase,Ih,NfilesInHDF);
            HDF5.save_cat(FileName,DataName,Out.Cat(Flag,:),2,30);
        end
    catch
        Ns     = NaN;
        Status = -1;
    end
    [Ihtm, Ih, toc, Ns, Status]
    
    if (Ihtm./100 == floor(Ihtm./100))
        save Update.mat Update    
    end
end

    