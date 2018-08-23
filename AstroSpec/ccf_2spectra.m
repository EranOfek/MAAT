function []=ccf_2spectra(Spec,Template1,Template2,varargin)
%--------------------------------------------------------------------------
% ccf_2spectra function                                          AstroSpec
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jun 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


DefV.VecZ1             = (0:0.01:2).';
DefV.VecZ2             = (2:0.01:4).';
DefV.Ebv               = 0;

InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

