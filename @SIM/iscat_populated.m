function IsPop=iscat_populated(Sim)
%--------------------------------------------------------------------------
% iscat_populated function                                      class/@SIM
% Description: Check if the catalog entry (.Cat) in a SIM object is
%              populated.
% Input  : - A SIM object.
% Output : - An array of logicals indicating if each element in the SIM
%            array has a populated catalog.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: IsPop=iscat_populated(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

CatField = 'Cat';

Nsim  = numel(Sim);
IsPop = false(size(Sim));
for Isim=1:1:Nsim,
    IsPop(Isim) = ~isempty(Sim(Isim).(CatField));
end