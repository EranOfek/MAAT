function []=add_meta_data2ps1(varargin)
% SHORT DESCRIPTION HERE
% Package: VO.PS1
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Dec 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------



DefV.FilterList          = {'PAN-STARRS/PS1','g';'PAN-STARRS/PS1','r';'PAN-STARRS/PS1','i';'PAN-STARRS/PS1','z';'PAN-STARRS/PS1','y'};
DefV.MagSystem           = 'AB';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


zVec = (0:0.01:1)';
Nz   = numel(zVec);

% prep synphot for all spectra

% galaxy/qso spectra
SpecName = AstSpec.get_galspec;
Nfilt = size(InPar.FilterList,1);
Nsn = numel(SpecName);
Ind = 0;
for Isn=1:1:Nsn
    Ind  = Ind + 1;
    AstS = AstSpec.get_galspec(SpecName{Isn});
    for Ifilt=1:1:Nfilt
        [Data(Ifilt,Ind).Mag,Data(Ind,Ifilt).Frac] = synphot(AstS,InPar.FilterList{Ifilt,1},InPar.FilterList{Ifilt,2},InPar.MagSystem);
        Data(Ifilt,Ind).SpecName = SpecName{Isn};
    end
end

% GAIA stellar spectra
W = Util.IO.load2('GAIA_Wave1A.mat');
Dir1 = dir('T*G*P*V000*.mat');
Dir2 = dir('T*G*M*V000*.mat');
Dir  = [Dir1;Dir2];
Nspec = numel(Dir);
for Ispec=1:1:Nspec
    Ind = Ind + 1;
    Spec = Util.IO.load2(Dir(Ispec).name);
    for Ifilt=1:1:Nfilt
        [Data(Ifilt,Ind).Mag,Data(Ind,Ifilt).Frac] = AstroUtil.spec.synphot([W, Spec],InPar.FilterList{Ifilt,1},InPar.FilterList{Ifilt,2},InPar.MagSystem);
        Data(Ifilt,Ind).SpecName = Dir(Ispec).name;
    end
end

% gal spec
Dir = dir('*_spec.dat');
Nspec = numel(Dir);
for Ispec=1:1:Nspec
    Ind = Ind + 1;
    Spec = Util.IO.load2(Dir(Ispec).name);
    for Iz=1:1:Nz
        Spec1 = Spec(:,1:2);
        Spec1(:,1) = Spec(:,1).*(1+zVec(Iz));
        Spec1(:,2) = Spec(:,2).*(1+zVec(Iz));  ????
        for Ifilt=1:1:Nfilt
            [Data(Ifilt,Ind).Mag,Data(Ind,Ifilt).Frac] = AstroUtil.spec.synphot(Spec1,InPar.FilterList{Ifilt,1},InPar.FilterList{Ifilt,2},InPar.MagSystem);
            Data(Ifilt,Ind).SpecName = Dir(Ispec).name;
        end
    end
end

% need to apply redshifts


% format Data into a matrix
Ndata  = numel(Data);
SynMag = zeros(Nfilt,Ndata);
for Ifilt=1:1:Nfilt
    SynMag(Ifilt,:) = [Data(Ifilt,:)];
end


% extinction steps 0.03, Rv step 0.1 -> in order to have 10mmag resolution


% load PS1 HTM file
Cat = 
Col.Mag    = [];  % g r i z y
Col.MagErr = [];  % g r i z y


% select bright objects with good photometry (about ~1000)
% phot err should be small [we wll use rms minimization due to this...]
FlagGP = 
IndGP  = find(FlagGP);
NGP    = numel(IndGP);


VecR   = (1:0.1:8)';
VecEbv = (0.000:0.003:0.1)';
Nr     = numel(VecR);
Nebv   = numel(VecEbv);


Ext = cell(Nfilt,1);
for Ifilt=1:1:Nfilt
    Ext{Ifilt} = zeros(Nebv,Nr);
    for Ir=1:1:Nr
        Ext{Ifilt}(:,Ir) = AstroUtil.spec.extinction(VecEbv,InPar.FilterList{Ifilt,1},InPar.FilterList{Ifilt,1},VecR(Ir));
    end
end

    
MedianMinRMS = zeros(Nebv,Nr);
for Ir=1:1:Nr
    for Iebv=1:1:Nebv
        
        Resid = zeros(Nfilt,NGP,Ndata);
        for Ifilt=1:1:Nfilt
            Resid(Ifilt,:,:) = Cat(FlagGP,Col.Mag) - (SynMag(Ifilt,:)+Ext{Ifilt}(Iebv,Ir));
        end

        % rms over all filters
        RMS = std(Resid,[],1);   % star, spec

        % for each star select minimum rms
        [MinRMS,MinI] = min(RMS,[],2);  % the minimum RMS for each star

        % remove MinRMS which are above some threshold
        % but what about the number of sources?
        
        
        % median of minimum rms
        MedianMinRMS(Iebv,Ir) = median(MinRMS);
    end
end

surface(VecEbv,VecR,MedianMinRMS);









