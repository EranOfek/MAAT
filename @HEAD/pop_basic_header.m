function Head=pop_basic_header(Head,varargin)
%--------------------------------------------------------------------------
% pop_basic_header function                                    class/@HEAD
% Description: Populate basic header into an HEAD object or a SIM image.
% Input  : - An HEAD object (or a SIM image).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ReplaceKey'- A three column cell array of {key,val,comment},
%                          or an HEAD object which keywords to replace
%                          or add to the SIM header.
%                          Default is {}.
%            'AddKey'    - Like 'ReplaceKey', but adding keywords,
%                          without replacment.
%                          Default is {}.
% Output : - An HEAD object with populated basic header.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: H=pop_basic_header(H);
% Reliable: 2
%--------------------------------------------------------------------------

Def.Origin    = 'MATLAB Astro';
ImageField    = 'Im';

DefV.ReplaceKey         = {};
DefV.AddKey             = {};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% basic header definition
Nh = numel(Head);
for Ih=1:1:Nh
    % for each header/image
    Ik = 0;
    if (SIM.issim(Head(Ih)))
        BitClass = class(Head(Ih).(ImageField));
        switch lower(BitClass)
            case {'double','uint64','int64'}
                BitPix = 64;
            case {'single','uint32','int32'}
                BitPix = 32;
            case {'uint16','int16'}
                BitPix = 16;
            case {'uint8','int8'}
                BitPix = 8;
            case {'logical'}
                BitPix = 1;
            otherwise
                error('Unknown data type');
        end
        Size     = size(Head(Ih).(ImageField));
        Naxis    = numel(Size);
    
        
        Ik = Ik + 1;
        KeyH(Ik,:) = {'SIMPLE', true,      'STANDARD FITS'};
        Ik = Ik + 1;
        KeyH(Ik,:) = {'BITPIX', BitPix,    'BITS/PIXEL'};
        Ik = Ik + 1;
        KeyH(Ik,:) = {'NAXIS',  Naxis,     'NUMBER OF AXES'};
        for Iaxis=1:1:Naxis
            if (Iaxis==1)
                Ind = 2;
            end
            if (Iaxis==2)
                Ind = 1;
            end
            Ik = Ik + 1;
            KeyH(Ik,:) = {sprintf('NAXIS%d',Iaxis), Size(Ind), sprintf('SIZE OF DIMENSION %d',Ind)};
        end
        Ik = Ik + 1;
        KeyH(Ik,:) = {'BSCALE',1,'PHYSICAL = BSCALE * DATA + BZERO'};
        Ik = Ik + 1;
        KeyH(Ik,:) = {'BZERO'  0,'BZERO'};
    end
    Ik = Ik + 1;
    KeyH(Ik,:) = {'ORIGIN',Def.Origin,'ORIGIN'};
    
    % populate the header
    Head(Ih) = replace_key(Head(Ih),KeyH);
end


%--- Update header ---
% replace
if (~isempty(InPar.ReplaceKey))
    Head = replace_key(Head,InPar.ReplaceKey);
end
% add
if (~isempty(InPar.AddKey))
    Head = add_key(Head,InPar.AddKey);
end
    