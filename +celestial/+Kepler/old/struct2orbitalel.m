function OrbEl=struct2orbitalel(St)
%--------------------------------------------------------------------------
% struct2orbitalel function                                          ephem
% Description: Convert a structure array into an OrbitalEl class object.
% Input  : - A structure array with the appropriate fields.
% Output : - An OrbitalEl object.
% See also: OrbitalEl.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: OrbEl=struct2orbitalel(Data);
% Reliable: 2
%--------------------------------------------------------------------------


OrbEl = OrbitalEl;
Nst   = numel(St);
FN    = fieldnames(St);
Nfn   = numel(FN);

for Ist=1:1:Nst,
    OrbEl(Ist) = OrbitalEl;
    for Ifn=1:1:Nfn,
        OrbEl(Ist).(FN{Ifn}) = St(Ist).(FN{Ifn});
    end
end

        