function Sim=astcat2sim(AstC,Sim)
%--------------------------------------------------------------------------
% astcat2sim function                                        class/@AstCat
% Description: Copy an AstCat object into an existing or new SIM object.
% Input  : - An AstCat object.
%          - An optional SIM object (with the same size as the AstCat
%            object) to which to copy the AstCat object.
% Output : - A SIM object with the AstCat catalog.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=astcat2sim(AstC);
% Reliable: 2
%--------------------------------------------------------------------------

CopyFields = {'Cat','Col','ColCell','ColUnits','SortedBy','SortedByCol'};

if (nargin==1),
    Sim = SIM(size(AstC));
end

Ncat = numel(AstC);
Nf   = numel(CopyFields);

for Icat=1:1:Ncat,
    % copy each AstCat object into the SIM
    for If=1:1:Nf,
        % copy each field
        Sim(Icat).(CopyFields{If}) = AstC(Icat).(CopyFields{If});
    end
end
    