function []=prep_binary_asteroid(varargin)
% Create a table of bknown binary asteroids
% Package: VO.prep
% Description: Create a table of bknown binary asteroids
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------



%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

URL = 'http://www.johnstonsarchive.net/astro/astmoontable.html';

Table = urlreadtable(URL);

Table = Table{1};
ColCell  = {'Name','Type','a_sun','P_sun','e','r1','rot1','a','p','r2','disc_method','disc_date','annon_date'};
ColUnits = {''    ,''    ,'au'   ,'yr'   ,'' ,'km','hr'  ,'km','day','km','','',''};
Table.Properties.VariableNames= ColCell;
Table.Properties.VariableUnits= ColUnits;
Table = Table(2:end,:);
Table.a_sun = str2double(Table.a_sun);
Table.P_sun = str2double(Table.P_sun);
Table.e     = str2double(Table.e);
Table.r1    = str2double(Table.r1);
Table.rot1  = str2double(Table.rot1);
Table.a     = str2double(Table.a);
Table.p     = str2double(Table.p);
Table.r2    = str2double(Table.r2);

BinaryMP = AstCat;
BinaryMP.Cat = Table;
BinaryMP.ColCell  = ColCell;
BinaryMP.ColUnits = ColUnits;
BinaryMP          = colcell2col(BinaryMP);
BinaryMP.Source   = URL;
 
Nast = size(BinaryMP.Cat,1);
AN=nan(N,1);
for I=1:1:Nast
    SSS=regexp(BinaryMP.Cat.Name{I},'\((?<num>\d+)\)','names');
    if (~isempty(SSS))
        AN(I) = str2double(SSS(1).num);
    end
end

