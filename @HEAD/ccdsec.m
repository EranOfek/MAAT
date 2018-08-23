function OutCCDSEC=ccdsec(Head,CCDSEC,UseNaN,GetKeyMethod)
% Get CCDSEC keyword value from an HEAD/SIM object.
% Package: @HEAD
% Description: Get CCDSEC keyword value from an HEAD object (or e.g.,
%              a SIM image).
% Input  : - An HEAD object (or e.g., a SIM images).
%          - Either a CCDSEC header keyword name from which to exctract
%            the CCDSEC (a string or a cell array of strings),
%            or a CCDSEC vector that will override the
%            header keyword. If empty, or header keyword is not available
%            then use the naxis command if input is HEAD and use the
%            the imagesize if input is SIM.
%            If this is a cell array of strings use getkey_fromlist.m
%            to extract the CCDSEC from several possible keywords.
%          - A flag indicating what to do if CCDSEC value is NaN.
%            If true then return NaN, if false then get CCDSEC using
%            naxis or imagesize (of input is SIM).
%            Default is false.
%          - Method by which to select the best keyword value.
%            See getkey_fromlist.m for options.
%            Default is 'first'.
% Output : - A 4 column Matrix with CCDSEC [Xmin, Xmax, Ymin, Ymax].
%            Row per header.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Oct 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: CCDSEC=ccdsec(H,'CCDSEC');
%          CCDSEC=ccdsec(H,{'CCDSEC','IMSEC'});
% Reliable: 2
%--------------------------------------------------------------------------


ImageField  = 'Im';
%HeaderField = 'Header';
%FileField   = 'ImageFileName';
%MaskField   = 'Mask';
%BackImField = 'BackIm';
%ErrImField  = 'ErrIm';

Def.CCDSEC = [];
Def.UseNaN = false;
Def.GetKeyMethod = 'first';
if (nargin==1)
    CCDSEC = Def.CCDSEC;
    UseNaN = Def.UseNaN;
    GetKeyMethod = Def.GetKeyMethod;
elseif (nargin==2)
    UseNaN = Def.UseNaN;
    GetKeyMethod = Def.GetKeyMethod;
elseif (nargin==3)
    GetKeyMethod = Def.GetKeyMethod;
elseif (nargin==4)
    % do nothing
else
    error('Illegal number of input arguments');
end

Nim = numel(Head);
OutCCDSEC = zeros(Nim,4);
for Iim=1:1:Nim
    if (isempty(CCDSEC))
        if (SIM.issim(Head))
            % get CCDSEC from image size
            Size = size(Head(Iim).(ImageField));
            OutCCDSEC(Iim,:) = [1 Size(2) 1 Size(1)];
        else
            % get CCDSEC using the naxis command
            Size = naxis(Head(Iim));
            OutCCDSEC(Iim,:) = [1 Size(1) 1 Size(2)];
        end
    else
        if (ischar(CCDSEC) || iscellstr(CCDSEC))
            % get CCDSEC from header
            [NewCellHead] = getkey_fromlist(Head(Iim),CCDSEC,GetKeyMethod);
            CCDSEC_Val = NewCellHead{1};

            if (isnan(CCDSEC_Val))
                if (UseNaN)
                    OutCCDSEC(Iim,:) = [NaN NaN NaN NaN];
                else
                    if (SIM.issim(Head))
                        % get CCDSEC from image size
                        Size = size(Head(Iim).(ImageField));
                        OutCCDSEC(Iim,:) = [1 Size(2) 1 Size(1)];
                    else
                        % get CCDSEC using the naxis command
                        Size = naxis(Head(Iim));
                        OutCCDSEC(Iim,:) = [1 Size(1) 1 Size(2)];
                    end
                  
                end
            else
                % CCDSEC is not NaN
                Splitted   = regexp(CCDSEC_Val,'[:\[\],]','split');
                OutCCDSEC(Iim,:)  = [str2double(Splitted{2}), str2double(Splitted{3}), str2double(Splitted{4}), str2double(Splitted{5})];
            end
        elseif (isnumeric(CCDSEC))
            % Use CCDSEC as is
            OutCCDSEC(Iim,:) = CCDSEC;
        else
            error('Unknown CCDSEC type');
        end
    end
end

    
        
        
    

