function []=read_table(File,varargin)
% SHORT DESCRIPTION HERE
% Package: VO.TNS
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Apr 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


DefV.Seperator            = ',';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

FID     = fopen(File,'r');
Line    = fgetl(FID);
ColCell = regexp(Line,InPar.Seperator,'split');
fclose(FID);
ColCell = regexprep(ColCell,'"','');
