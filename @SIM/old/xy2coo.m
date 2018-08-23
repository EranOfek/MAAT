function Coo=xy2coo(Sim,varargin)
%--------------------------------------------------------------------------
% xy2coo function                                               class/@SIM
% Description: Given a SIM image with a WCS in header, convert X,Y pixel
%              coordinates to RA/Dec. Using xy2sky.m.
% Input  : - SIM class object.
%          - Vector of X positions.
%          - Vector of Y poistions.
%          * Additional arguments to pass to xy2sky.m.
% Output : - AstCat object in which each catalog has two columns [RA,Dec].
%            Coordinates are given in radians.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:  S.xy2coo(1,1)
% Reliable: 2
%--------------------------------------------------------------------------

Nsim = numel(Sim);
%Coo  = struct('RA',cell(Nsim,1),'Dec',cell(Nsim,1));
Coo  = astcatdef(size(Sim));
for Isim=1:1:Nsim,
    [Coo(Isim).Cat(:,1),Coo(Isim).Cat(:,2)] = xy2sky(Sim(Isim),varargin{:});
    Coo(Isim).ColCell  = {'RA','Dec'};
    Coo(Isim).ColUnits = {'rad','rad'};
    Coo(Isim)          = colcell2col(Coo);
end