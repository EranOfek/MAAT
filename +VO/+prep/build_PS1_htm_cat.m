function build_PS1_htm_cat
% build PS1 HDF5/HTM catalog

RAD = 180./pi;

Col.RA  = 1;
Col.Dec = 2;

Level = 9;
NfilesInHDF = 1000;
[HTM,LevList]=celestial.htm.htm_build(Level);

Ind = find([HTM.level] == (Level - 1));
Nind = numel(Ind);

% addpath('/raid/eran/Catalogue/PS1/PS_Tile/Z2');
% addpath('/raid/eran/Catalogue/PS1/PS_Tile/Z3');
% addpath('/raid/eran/Catalogue/PS1/PS_Tile/Z4');
% addpath('/raid/eran/Catalogue/PS1/PS_Tile/Z5');
% addpath('/raid/eran/Catalogue/PS1/PS_Tile/Z6');
% addpath('/raid/eran/Catalogue/PS1/PS_Tile/Z7');

load('SummaryMaster.txt')  % SummaryMaster [Tile #src]
SummaryMaster=sortrows(SummaryMaster,1);
load('TileListPS1.mat');   % TileList [CenterRA,CenterDec,MinRA,MaxRA,MinDec,MaxDec]  [deg]
TileSize = 0.2;  % [deg]

% missing fields
% III=TileList(II,2)>-15;      
% >> II(III)
% [II(III), TileList(II(III),:)]
% 
% ans =
% 
%       160097
%       160107
%       268070
%       268071
%       268072
%       269869
%       269870
%       269871
%       269872
%       271669
%       271670
%       271671
%       271672
%       273469
%       273470
%       273471
%       273472
%       275269
%       431423
%       708674
%       711157
%       711167

NfilesInHDF = 100;
FileBase    = 'PS1';
Update      = zeros(0,3);
for Iind=442369:1:Nind
    Iind
    tic;
    Ihtm = Ind(Iind);
    
    MeanRA  = mean(HTM(Ihtm).coo(:,1));   % [rad]
    MeanDec = mean(HTM(Ihtm).coo(:,2));   % [rad]
    
    D = celestial.coo.sphere_dist_fast(HTM(Ihtm).coo(:,1),HTM(Ihtm).coo(:,2),MeanRA,MeanDec);
    SearchRadius = max(D).*1.5 + TileSize.*0.5.*1.5./RAD;
    
    % search for all relevant tiles near HTM
    Dt = celestial.coo.sphere_dist_fast(MeanRA,MeanDec,TileList(:,1)./RAD,TileList(:,2)./RAD);
    
    % list of tiles to load
    Itile = find(Dt<(SearchRadius) & SummaryMaster(:,2)>0);
    Ntile = numel(Itile);
    Cat = zeros(0,41);
    for It=1:1:Ntile
        %It
        %[Itile(It), SummaryMaster(Itile(It),:)]
        TileFileName = sprintf('PS_Tile_%06d.hdf5',Itile(It));
        Cat=[Cat; HDF5.load(TileFileName,'V')];
        
    end   
    
    % convert RA/Dec to radians
    Cat(:,[Col.RA, Col.Dec]) = Cat(:,[Col.RA, Col.Dec])./RAD;
    
    % search for stars within HTM
    MaxErr = 0.1;
    Flag = celestial.htm.in_polysphere(Cat(:,[Col.RA, Col.Dec]),HTM(Ihtm).coo,2) & ...
           (abs(Cat(:,8))<MaxErr | abs(Cat(:,15))<MaxErr | abs(Cat(:,22))<MaxErr | abs(Cat(:,29))<MaxErr | abs(Cat(:,36))<MaxErr);
    Cat = Cat(Flag,:);
    
    [FileName,DataName]=HDF5.get_file_var_from_htmid(FileBase,Ihtm,NfilesInHDF);
   
    if (~isempty(Cat))
        HDF5.save_cat(FileName,DataName,Cat,2,30);
    end
    
    Ns = size(Cat,1);
    Update = [Update; [Iind, Ihtm, Ns]];
    if (Iind./100==floor(Iind./100))
        save Update.mat Update
    end
    
    toc
    
end

