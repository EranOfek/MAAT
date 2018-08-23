function Flag=isfield_populated(Sim,Field)
% Check if a given field in a SIM object is populated.
% Package: @SIM
% Description: Check if a given field in a SIM object is populated.
% Input  : - A SIM object.
%          - A field name (e.g., 'Im').
% Output : - Logical flags indicating if the field, for each SIM element,
%            is not empty.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Flag=isfield_populated(Sim,'BackIm');
% Reliable: 2
%--------------------------------------------------------------------------


Nsim = numel(Sim);
Flag = false(size(Sim));
for Isim=1:1:Nsim
    %Flag(Isim) = ~isempty(Sim(Isim).(Field));
    Flag(Isim) = any(size(Sim(Isim).(Field))>0);
end
