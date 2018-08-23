function Coo=coo2xy(Sim,varargin)
%--------------------------------------------------------------------------
% coo2xy function                                               class/@SIM
% Description: Given a SIM image with a WCS in header, convert RA,Dec
%              coordinates [radians] to X,Y pixel coordinates.
%              Using xy2sky.m.
% Input  : - SIM class object.
%          - Vector of RA positions [rad].
%          - Vector of Dec poistions [rad].
%          * Additional arguments to pass to sky2xy.m.
% Output : - AstCat object in which each catalog has two columns [X,Y].
%            Coordinates are given in pixels.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:  S.coo2xy(1,1)
% Reliable: 3
%--------------------------------------------------------------------------

Nsim = numel(Sim);
%Coo  = struct('RA',cell(Nsim,1),'Dec',cell(Nsim,1));
Coo  = astcatdef(size(Sim));
for Isim=1:1:Nsim,
    [Coo(Isim).Cat(:,1),Coo(Isim).Cat(:,2)] = sky2xy(Sim(Isim),varargin{:});
    Coo(Isim).ColCell  = {'X','Y'};
    Coo(Isim).ColUnits = {'pix','pix'};
    Coo(Isim)          = colcell2col(Coo);
end