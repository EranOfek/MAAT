function WCS=get_wcs(Sim,OutType)
%--------------------------------------------------------------------------
% get_wcs function                                             class/@HEAD
% Description: Get the World Coordinate System (WCS) from an HEAD object
%              (e.g., an header of a FITS image).
% Input  : - An HEAD object, or e.g., SIM images.
%          - Output type: 'struct'|'WorldCooSys'. Default is 'WorldCooSys'.
% Output : - A structure array, or an WorldCooSys object
%            containing the WCS keywords.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: WCS=get_wcs(S);
% Reliable: 2
%--------------------------------------------------------------------------

HeaderField   = 'Header';

if (nargin<2)
    OutType = 'WorldCooSys';
end

KeywordCell{1}  = 'RADECSYS';
KeywordCell{2}  = 'CTYPE1';
KeywordCell{3}  = 'CTYPE2';
KeywordCell{4}  = 'CUNIT1';
KeywordCell{5}  = 'CUNIT2';
KeywordCell{6}  = 'CRPIX1';
KeywordCell{7}  = 'CRPIX2';
KeywordCell{8}  = 'CRVAL1';
KeywordCell{9}  = 'CRVAL2';
KeywordCell{10} = 'CD1_1';
KeywordCell{11} = 'CD1_2';
KeywordCell{12} = 'CD2_1';
KeywordCell{13} = 'CD2_2';
KeywordCell{14} = 'CDELT1';
KeywordCell{15} = 'CDELT2';
KeywordCell{16} = 'PC1_1';
KeywordCell{17} = 'PC1_2';
KeywordCell{18} = 'PC2_1';
KeywordCell{19} = 'PC2_2';
KeywordCell{20} = 'LONPOLE';
KeywordCell{21} = 'LATPOLE';
KeywordCell{22} = 'EQUINOX';

Nim = numel(Sim);

% keyval from an Header object
if (any(isempty_wcs(Sim)))
    WCS = keyval2struct(Sim,KeywordCell,{},true);
    switch lower(OutType)
        case 'worldcoosys'
            WCS = WorldCooSys.struct2wcs(WCS);
    end  
else
    WCS = WorldCooSys(Nim,1);
    for Iim=1:1:Nim
        WCS(Iim).WCS = Sim(Iim).WCS;
    end
end



for Iim=1:1:Nim
    %--- SIP/PV ---
    % Get the sip distortion information 
    if ((~isempty(strfind(WCS(Iim).WCS.CTYPE1,'SIP')))||(~isempty(strfind(WCS(Iim).WCS.CTYPE2,'SIP')))) 
        WCS(Iim).WCS.sip = FITS.get_sip(Sim(Iim).(HeaderField));
    end
    % Get the tpv distortion information 
    if ((~isempty(strfind(WCS(Iim).WCS.CTYPE1,'TPV')))||(~isempty(strfind(WCS(Iim).WCS.CTYPE2,'TPV')))) 
        WCS(Iim).WCS.tpv = FITS.get_tpv(Sim(Iim).(HeaderField));
    end
end

%-----------------------
%--- fix WCS key/val ---
%-----------------------
for Iim=1:1:Nim
    
    if (isnan(WCS(Iim).WCS.CUNIT1))
       WCS(Iim).WCS.CUNIT1 = 'deg';
    end
    if (isnan(WCS(Iim).WCS.CUNIT2))
       WCS(Iim).WCS.CUNIT2 = 'deg';
    end

    WCS(Iim).WCS.CUNIT1 = Util.string.spacedel(WCS(Iim).WCS.CUNIT1);
    WCS(Iim).WCS.CUNIT2 = Util.string.spacedel(WCS(Iim).WCS.CUNIT2);
    
    if (~isnan(WCS(Iim).WCS.PC1_1))
        % use PC matrix instead of CD matrix
        % see Eq. 1 in Calabertta & Greisen (2002)
        if (isnan(WCS(Iim).WCS.CD1_1))
            WCS(Iim).WCS.CD1_1 = WCS(Iim).WCS.CDELT1.*WCS(Iim).PC1_1;
            WCS(Iim).WCS.CD1_2 = WCS(Iim).WCS.CDELT1.*WCS(Iim).WCS.PC1_2;
            WCS(Iim).WCS.CD2_1 = WCS(Iim).WCS.CDELT2.*WCS(Iim).WCS.PC2_1;
            WCS(Iim).WCS.CD2_2 = WCS(Iim).WCS.CDELT2.*WCS(Iim).WCS.PC2_2;
        else
            % ignore PC matrix sibcd CD matrix is given too
        end
    end


    if (isnan(WCS(Iim).WCS.CD1_1))
       % try to use CDELT1/2
       WCS(Iim).WCS.CD1_1 = WCS(Iim).WCS.CDELT1; %CellWCS{14};
       WCS(Iim).WCS.CD2_2 = WCS(Iim).WCS.CDELT2; %CellWCS{15};
       WCS(Iim).WCS.CD1_2 = 0;
       WCS(Iim).WCS.CD2_1 = 0;
    end

    WCS(Iim).WCS.CD = [WCS(Iim).WCS.CD1_1, WCS(Iim).WCS.CD1_2; WCS(Iim).WCS.CD2_1, WCS(Iim).WCS.CD2_2];
end


      
