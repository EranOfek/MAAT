function []=prep_decals_htmcat(varargin)
% SHORT DESCRIPTION HERE
% Package: VO
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

CatField = AstCat.CatField;

DefV.DataURL              = 'ftp://archive.noao.edu/public/hlsp/ls/dr5/sweep/5.0/';
DefV.Download             = true;
DefV.MakeCat              = true;
DefV.LevelHTM             = 9;
DefV.BricksCatDR          = 'cats.DECaLS.DECaLS_bricks_DR5';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% download files
if (InPar.Download)
    % download sweep files
    [List,FileNames] = www.find_urls_ftp(InPar.DataURL);
end

if (InPar.MakeCat)
    % make HDF5/HTM from sweep files
    
    [HTM,LevelH] = celestial.htm.htm_build(InPar.LevelHTM);
    LevelL = LevelH(InPar.LevelHTM);
    Nhtm   = numel(LevelL.ptr);
    
    % load bricks catalog
    DECaLS_bricks_DR  = eval(InPar.BricksCatDR);
    DECaLS_bricks     = cats.DECaLS.DECaLS_bricks;
    
    Nbrick = numel(DECaLS_bricks_DR.(CatField),1);
    
   
    for Ibrick=1:1:Nbrick
        % load all adjuscent bricks
        BrickRA  = DECaLS_bricks_DR.(CatField)(Ibrick,1);
        BrickDec = DECaLS_bricks_DR.(CatField)(Ibrick,2);
        D = celestial.coo.sphere_dist(BrickRA,BrickDec,...
                                           DECaLS_bricks.(CatField)(:,1),...
                                           DECaLS_bricks.(CatField)(:,2),'deg');
        Ib = find(D<0.001);  % brick index
        BrickCol = DECaLS_bricks.(CatField)(Ib,10);
        BrickRow = DECaLS_bricks.(CatField)(Ib,9);
        
        % Bricks in the same or nearby declination zones
        I1 = find(DECaLS_bricks.(CatField)(:,9)==BrickRow | ...
                  DECaLS_bricks.(CatField)(:,9)==(BrickRow-1) | ...
                  DECaLS_bricks.(CatField)(:,9)==(BrickRow+1));
              
        BrickRAwidth = DECaLS_bricks.(CatField)(Ib,4) - DECaLS_bricks.(CatField)(Ib,3);
        BrickWidth   = BrickRAwidth.*cosd(DECaLS_bricks.(CatField)(Ib,2));
        D = celestial.coo.sphere_dist(BrickRA,BrickDec,...
                                      DECaLS_bricks.(CatField)(I1,1),...
                                      DECaLS_bricks.(CatField)(I1,2),'deg');
         
        Iadj = I1(D<(BrickWidth.*1.5./RAD));
        
        % load all Iadj bricks
        Nadj = numel(Iadj);
        for Ia=1:1:Nadj
            I = Iadj(Ia);
            % prep brick file name
            BrickMin = DECaLS_bricks.(CatField)(I,
            BrickMax = 
            sprintf('sweep-%03d%c%03d-%03d%c%03d.fits
        
        % find HTMs which centers within Ib brick
        
        
        
        
    
end

    
