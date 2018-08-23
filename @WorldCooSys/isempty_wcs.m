function IsEmpty=isempty_wcs(WCS)
% Check if the WCS fields in a WorldCooSys array object are empty.
% Package: @WorldCooSys
% Description: Check if the WCS fields in a WorldCooSys array object are 
%              empty.
% Input  : - An WorldCooSys object.
% Output : - A flag (for each elelemt in the input WorldCooSys object),
%            indictaing if its WCS field is empty.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: IsEmpty=isempty_wcs(WCS)
% Reliable: 
%--------------------------------------------------------------------------


Nw = numel(WCS);
IsEmpty = true(size(WCS));
for Iw=1:1:Nw
    IsEmpty(Iw) = isempty(WCS(Iw).WCS);
end