function Cell=sim2cell(Sim,Field)
% Convert the images in SIM into a cell array.
% Package: @SIM
% Description: Convert the images in SIM into a cell array.
% Input  : - A SIM object.
%          - A field name which content to convert into a cell array.
%            Default is 'Im'.
% Output : - A cell array.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jul 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cell=sim2cell(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==1)
    Field = SIM.ImageField;
end

Nsim = numel(Sim);
Cell = cell(size(Sim));
for Isim=1:1:Nsim
    Cell{Isim} = Sim.(Field);
end