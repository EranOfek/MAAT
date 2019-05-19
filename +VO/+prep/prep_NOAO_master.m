function []=prep_NOAO_master(varargin)
% SHORT DESCRIPTION HERE
% Package: VO.prep
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;

DefV.DirNum               = [1];
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Ndir = numel(InPar.DurNum);
for Idir=1:1:Ndir
    DirStr = sprintf('%d',InPar.DirNum(Idir));
    cd(DirStr);
    
    Files  = dir('*.fits');
    Nfiles = numel(Files);
    for Ifiles=1:1:Nfiles
        T = FITS.read_table(Files(Ifiles).name);
        %Cols = [T.Col.RA, T.Col.DEC];
        RA   = T.Cat(:,T.Col.RA)./RAD;
        Dec  = T.Cat(:,T.Col.DEC)./RAD;
        Mjd  = T.Cat(:,T.Col.MJD);
        FWHM = T.Cat(:,T.Col.FWHM);
        CCDNUM = T.Cat(:,T.Col.CCDNUM);
        FILTER = T.Cat(:,T.Col.FILTER);
        X      = T.Cat(:,T.Col.X);
        Y      = T.Cat(:,T.Col.Y);
        
        [X,Y,Z] = celestial.coo.coo2cosined(RA,Dec);
        [MidRA,MidDec] = celestial.coo.cosined2coo(median(X),median(Y),median(Z));
        Dist    = celestial.coo.sphere_dist_fast(MidRA,MidDec,RA,Dec);
        MaxRad  = max(Dist);
        
        
        
        