function ImSize=naxis(Head)
%--------------------------------------------------------------------------
% naxis function                                               class/@HEAD
% Description: Get the value of the NAXIS1, NAXIS2,... keywords from
%              the header.
% Input  : - An HEAD object, or e.g., a SIM object.
% Output : - A matrix of NAXIS value keywords. each line corresponds to
%            on HEAD element (or e.g., a SIM image). The number of columns
%            is equal to the maximum NAXIS value, and the columns
%            corresponds to NAXIS1, NAXIS2, etc.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: ImSize=naxis(S)
% Reliable: 2
%--------------------------------------------------------------------------

Nh       = numel(Head);
Naxis    = Head.getkey('NAXIS');
if (any(isnan(cell2mat(Naxis))))
    % Naxis is not available for all images - return empty
    ImSize = [];
else
    MaxNaxis = max([Naxis{:}]);
    ImSize   = nan(Nh,MaxNaxis);

    for Ih=1:1:Nh
        for Inaxis=1:1:Naxis{Ih}
            Tmp = Head(Ih).getkey(sprintf('NAXIS%d',Inaxis));
            ImSize(Ih,Inaxis) = Tmp{1};
        end
    end
end
    