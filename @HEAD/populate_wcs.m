function Head=populate_wcs(Head,WCS)
%--------------------------------------------------------------------------
% populate_wcs function                                        class/@HEAD
% Description: Populate the WorldCooSys object in an HEAD object.
%              This can be used to populate the WCS field
%              in an HEAD object by copying it from an input structure or
%              by extracting the WCS information from the header.
% Input  : - An HEAD object.
%          - An optional WCS object or structure. If empty, will attempt
%            to extract the WCS information from the HEAD object
%            Header field. Default is empty,
% Output : - An HEAD object in which the WCS field is populated.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Head=populate_wcs(S);
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<2)
    WCS = [];
end


if (~isempty(WCS))
    % WCS structure is provided as input
    Nw = numel(WCS);
    Nh = numel(Head);
    
    N = max(Nw,Nh);
    for I=1:1:N
        Iw = min(N,Nw);
        Ih = min(N,Nh);
        
        Head(Ih).WCS      = WCS(Iw).WCS;
        Head(Ih).UserData = WCS(Iw).UserData;
    end
else
    % WCS was not provided
    % attmpt to get it from header
    WCS = get_wcs(Head);
    Nw = numel(WCS);
    for Iw=1:1:Nw
        Head(Iw).WCS      = WCS(Iw).WCS;
        Head(Iw).UserData = WCS(Iw).UserData;
    end
end

    