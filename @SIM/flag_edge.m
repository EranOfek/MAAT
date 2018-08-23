function Sim=flag_edge(Sim,varargin)
% Set a bit mask for pixels near image edge.
% Package: @SIM
% Description: Set a bit mask for pixels near image edge.
% Input  : - A SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'BufferWidth' - pixels located up to this distance from image
%                       edge will be flagged. Default is 12.
%            'BitEdgeName' - Edge bit mask name. Default is 'Bit_Edge'.
% Output : - A SIM object.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=flag_edge(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

ImageField = SIM.ImageField;

DefV.BufferWidth          = 12;
DefV.BitEdgeName          = 'Bit_Edge';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Nsim = numel(Sim);
for Isim=1:1:Nsim
    SizeIm = size(Sim(Isim).(ImageField));
    [MatX,MatY] = meshgrid((1:1:SizeIm(2)),(1:1:SizeIm(1)));
    FlagIm = false(SizeIm);
    FlagIm(MatX<InPar.BufferWidth | MatY<InPar.BufferWidth | MatX>(SizeIm(2)-InPar.BufferWidth) | MatY>(SizeIm(1)-InPar.BufferWidth)) = true;
    
    Sim(Isim) = bitmask_set(Sim(Isim),FlagIm,InPar.BitEdgeName);
end
