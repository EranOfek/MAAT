function []=loop_ndgrid(t,Obs,Fun,VarNames,varargin)
% SHORT DESCRIPTION HERE
% Package: Util.fit
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------



Ngrid = numel(varargin);

for Igrid=1:1:Ngrid
    N = numel(varargin{Igrid});
    
    for In=1:1:N
        if (isempty(VarNames))
        Y = Fun(t,vara
    