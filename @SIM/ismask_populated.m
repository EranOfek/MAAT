function IsPop=ismask_populated(Sim)
% Check if the mask entry (.Mask) in a SIM object is populated.
% Package: @SIM
% Description: Check if the mask entry (.Mask) in a SIM object is
%              populated.
% Input  : - A SIM object.
% Output : - An array of logicals indicating if each element in the SIM
%            array has a populated mask.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: IsPop=ismask_populated(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

MaskField = 'Mask';

Nsim  = numel(Sim);
IsPop = false(size(Sim));
for Isim=1:1:Nsim
    IsPop(Isim) = ~isempty(Sim(Isim).(MaskField));
end