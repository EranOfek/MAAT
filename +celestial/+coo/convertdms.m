function [OutData]=convertdms(InData,InType,OutType)
% Convert between various representations of coordinates and time 
% Package: celestial.coo
% Description: Convert between various representations of coordinates and
%              time as sexagesimal coordinates, degrees and radians.
% Input  : - Input data in one of the following formats:
%            degrees, radians, [H M S.S], [Sign D M S.S],
%            [H M.M], [Sign D M.M], 
%            'HH:MM:SS.S', '+DD:MM:SS.S'.
%          - Input type:
%            'r'  - radians.
%            'd'  - degrees.
%            'D'  - [Sign D M S.S] format.
%            'DM' - [Sign D M.M] format.
%            'h'  - hours.
%            'H'  - [H M S.S] format.
%            'HM' - [H M.M] format.
%            'f'  - fraction in range [0,1].
%            'SH' - String of sexagesimal hours 'HH:MM:SS.SSS'
%                   In this case the input can be strings or cell array
%                   of strings.
%            'SHb'- like 'SH', but with blanks instead of ":".
%            'SHh'- String of sexagesimal hours 'HHhMMmSS.SSSs'
%                   In this case the input can be strings or cell array
%                   of strings.
%            'SD' - String of sexagesimal degrees '+DD:MM:SS.SS'
%                   In this case the input can be strings or cell array
%                   of strings.
%            'SDb'- like 'SH', but with blanks instead of ":".
%            'gH' - Automatic identification of hours,
%                   [rad], [H M S], [sexagesimal string].
%            'gD' - Automatic identification of degrees,
%                   [rad], [Sign D M S], [sexagesimal string].
%          - Output type:
%            'r'  - radians, in range [0,2*pi] (default).
%            'R'  - radians, in range [-pi,pi].
%            'd'  - degrees.
%            'h'  - hours.
%            'f'  - fraction in range [0,1].
%            'D'  - [Sign D M S.S] format.
%            'H'  - [H M S.S] format.
%            'DM' - [Sign D M.M] format.
%            'HM' - [H M.M] format.
%            'SH' - String of sexagesimal hours 'HH:MM:SS.SSS'
%            'SD' - String of sexagesimal degrees '+DD:MM:SS.SS'
%            'SHn'- String of sexagesimal hours 'HHMMSS.SSS'
%            'SHh'- String of sexagesimal hours 'HHhMMmSS.SSSs'
%            'SDn'- String of sexagesimal degrees '+DDMMSS.SS'
%            'SHb'- like 'SH', but with blanks instead of ":".
%            'SDb'- like 'SH', but with blanks instead of ":".
% Output : - Requested output.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Jun 2000
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: celestial.coo.convertdms([1;2;3],'gH','d');
%          celestial.coo.convertdms([-1;0.1;1],'gD','SD');
% Reliable: 2
%--------------------------------------------------------------------------
import celestial.coo.*

RAD     = 180./pi;
Epsilon = 1e-8;
IsCell  = 0;
SignVec = ['-';'+'];



%InSize=size(InData);
%if (InSize(2)==1),
%   InType = InType;
%elseif (InSize(2)==3),
%   InType = 'H';
%elseif (InSize(2)==4),
%   InType = 'D';
%else
%   error('InData has illegal size');
%end

%if (nargin==1),
%   %InType  = InType;
%   OutType = 'r';   
if (nargin==2)
   OutType = 'r';
elseif (nargin==3)
   % no default
else
   error('Illegal number of input arguments');
end

%--- automatic type recognition ---
switch InType
 case 'gH'
    if (ischar(InData)==1 || iscellstr(InData)==1)
       InType = 'SH';
    else
       if (size(InData,2)==1)
          InType = 'r';
       elseif (size(InData,2)==2)
          InType = 'HM';
       elseif (size(InData,2)==3)
          InType = 'H';
       else
          error('Unknown InData format');
       end
    end
 case 'gD'
    if (ischar(InData)==1 || iscell(InData)==1)
       InType = 'SD';
    else
       if (size(InData,2)==1)
          InType = 'R';
       elseif (size(InData,2)==3)
          InType = 'DM';
       elseif (size(InData,2)==4)
          InType = 'D';
       else
          error('Unknown InData format');
       end
    end
 otherwise
    % do nothing
end



% convert InData to radians
switch InType
 case {'D'}
    RadData = (InData(:,2)+InData(:,3)./60+InData(:,4)./3600).*InData(:,1)./RAD;
 case {'DM'}
    RadData = (InData(:,2)+InData(:,3)./60).*InData(:,1)./RAD;
 case {'H'}
    RadData = (InData(:,1)+InData(:,2)./60+InData(:,3)./3600).*15./RAD;
 case {'HM'}
    RadData = (InData(:,1)+InData(:,2)./60).*15./RAD;
 case {'R','r'}
    RadData = InData;
    % do nothing
 case {'d'}
    RadData = InData./RAD;
 case {'h'}
    RadData = InData.*15./RAD;
 case {'f'}
    RadData = InData.*2.*pi;
 case {'SH'}
    IsCell  = iscell(InData);
    InData  = char(InData);     % convert cell array to verticaly con-strings
    N       = size(InData,1);
    RadData = zeros(N,1);
    for I=1:1:N
       [H,M,S] = strread(InData(I,:),'%2d:%2d:%f');
       RadData(I) = (H + M./60 + S./3600).*15./RAD;
    end
 case {'SHn'}
    IsCell  = iscell(InData);
    InData  = char(InData);     % convert cell array to verticaly con-strings
    N       = size(InData,1);
    RadData = zeros(N,1);
    for I=1:1:N
       [H,M,S] = strread(InData(I,:),'%2d%2d%f');
       RadData(I) = (H + M./60 + S./3600).*15./RAD;
    end
 case {'SHh'}
    IsCell  = iscell(InData);
    InData  = char(InData);     % convert cell array to verticaly con-strings
    N       = size(InData,1);
    RadData = zeros(N,1);
    for I=1:1:N
       [H,M,S] = strread(InData(I,:),'%2dh%2dm%fs');
       RadData(I) = (H + M./60 + S./3600).*15./RAD;
    end

 case {'SHb'}
    IsCell  = iscell(InData);
    InData  = char(InData);     % convert cell array to verticaly con-strings
    N       = size(InData,1);
    RadData = zeros(N,1);
    for I=1:1:N
       [H,M,S] = strread(InData(I,:),'%2d %2d %f');
       RadData(I) = (H + M./60 + S./3600).*15./RAD;
    end
 case {'SD'}
    IsCell  = iscell(InData);
    InData  = char(InData);     % convert cell array to verticaly con-strings
    N       = size(InData,1);
    RadData = zeros(N,1);
    for I=1:1:N
       [Sign,D,M,S] = strread(InData(I,:),'%c%2d:%2d:%f');
       switch Sign
        case '+'
           DecSign = +1;
        case '-'
           DecSign = -1;
        otherwise
           error('Unknown declination sign');
       end
       RadData(I) = (D + M./60 + S./3600).*DecSign./RAD;
    end
 case {'SDn'}
    IsCell  = iscell(InData);
    InData  = char(InData);     % convert cell array to verticaly con-strings
    N       = size(InData,1);
    RadData = zeros(N,1);
    for I=1:1:N
       [Sign,D,M,S] = strread(InData(I,:),'%c%2d%2d%f');
       switch Sign
        case '+'
           DecSign = +1;
        case '-'
           DecSign = -1;
        otherwise
          error('Unknown declination sign');
       end
       RadData(I) = (D + M./60 + S./3600).*DecSign./RAD;
    end
    
case {'SDh'}
    IsCell  = iscell(InData);
    InData  = char(InData);     % convert cell array to verticaly con-strings
    N       = size(InData,1);
    RadData = zeros(N,1);
    for I=1:1:N
       [Sign,D,M,S] = strread(InData(I,:),'%c%2dd%2dm%fs');
       switch Sign
        case '+'
           DecSign = +1;
        case '-'
           DecSign = -1;
        otherwise
              error('Unknown declination sign');
       end
       RadData(I) = DecSign.*(D + M./60 + S./3600)./RAD;
    end
 

 case {'SDb'}
    IsCell  = iscell(InData);
    InData  = char(InData);     % convert cell array to verticaly con-strings
    N       = size(InData,1);
    RadData = zeros(N,1);
    for I=1:1:N

       % This line is in order to correct for a bug in Horizons
       % replacing ** by 00
       InData(I,:) = strrep(InData(I,:),'**','00');

       [Sign,D,M,S] = strread(InData(I,:),'%c%2d %2d %f');
       switch Sign
        case '+'
           DecSign = +1;
        case '-'
           DecSign = -1;
        otherwise
           error('Unknown declination sign');
       end
       RadData(I) = (D + M./60 + S./3600).*DecSign./RAD;
    end
 otherwise
    error('Illegal InType');
end

% convert RadData to InType


%if (~strcmp(OutType,'D') || ~strcmp(InType,'D')),
%   RadData = 2.*pi.*(RadData./(2.*pi) - floor(RadData./(2.*pi)));
%end


switch OutType
 case {'D'}
    Sign    = sign(RadData);
    RadData = abs(RadData*RAD);
    D       = floor(RadData);
    M       = (RadData - D).*60;
    S       = (M - floor(M)).*60;
    M       = floor(M);
    I    = find(S>=60);  %abs(S-60)<Epsilon);
    S(I) = 0;
    M(I) = M(I) + 1;
    I    = find(M>=60);  %abs(M-60)<Epsilon);
    M(I) = 0;
    D(I) = D(I) + 1;
    OutData = [Sign, D, M, S];
 case {'H'}
    H       = RadData.*RAD./15;
    M       = (H - floor(H)).*60;
    S       = (M - floor(M)).*60;
    H       = floor(H);
    M       = floor(M);
    I    = find(S>=60);    %abs(S-60)<Epsilon);
    S(I) = 0;
    M(I) = M(I) + 1;
    I    = find(M>=60);   %abs(M-60)<Epsilon);
    M(I) = 0;
    H(I) = H(I) + 1;
    OutData = [H M S];
 case {'DM'}
    Sign    = sign(RadData);
    RadData = abs(RadData*RAD);
    D       = floor(RadData);
    M       = (RadData - D).*60;
    I    = find(M>=60);  %abs(M-60)<Epsilon);
    M(I) = 0;
    D(I) = D(I) + 1;
    OutData = [Sign, D, M];
 case {'HM'}
    H       = RadData.*RAD./15;
    M       = (H - floor(H)).*60;
    H       = floor(H);
    I    = find(M>=60);   %abs(M-60)<Epsilon);
    M(I) = 0;
    H(I) = H(I) + 1;
    OutData = [H M];
 case {'r'}
    % do nothing, allready in range [0,2*pi]
    OutData = RadData;
 case {'R'}
    RadData = 2.*pi.*(RadData./(2.*pi) - floor(RadData./(2.*pi)));
    I = find(RadData>pi);
    OutData = RadData;
    OutData(I) = RadData(I) - 2.*pi;
 case {'f'}
    OutData = RadData./(2.*pi);
 case {'d'}
    OutData = RadData.*RAD;
 case {'h'}
    OutData = RadData.*RAD./15;
 case {'SH','SHn','SHb','SH7'}
    RA = convertdms(RadData,'r','H');
    N  = size(RA,1);

    %StrRA = Util.cell.sprintf2cell('%02d:%02d:%06.3f',RA);  % faster to use...
    for I=1:1:N
       switch OutType
        case 'SH'
           StrRA = sprintf('%02d:%02d:%06.3f',RA(I,:));
        case 'SH7'
           StrRA = sprintf('%02d:%02d:%011.7f',RA(I,:));
        case 'SHn'
           StrRA = sprintf('%02d%02d%06.3f',RA(I,:));
        case 'SHb'
           StrRA = sprintf('%02d %02d %06.3f',RA(I,:));
        otherwise
           error('Unknown  OutType Option');
       end
       if (IsCell==1)
          OutData{I}   = StrRA;
       else 
          OutData(I,:) = StrRA;
       end
    end
 case {'SD','SDn','SDb','SD6'}
    Dec = convertdms(RadData,'r','D');
    DecSign = SignVec(floor(0.5.*Dec(:,1)+1.5+eps));
    N  = size(Dec,1);
    for I=1:1:N
       switch OutType
        case 'SD'
           StrDec = sprintf('%s%02d:%02d:%05.2f',DecSign(I),Dec(I,2:4));
        case 'SD6'
           StrDec = sprintf('%s%02d:%02d:%010.6f',DecSign(I),Dec(I,2:4));   
        case 'SDn'
           StrDec = sprintf('%s%02d%02d%05.2f',DecSign(I),Dec(I,2:4));
        case 'SDb'
           StrDec = sprintf('%s%02d %02d %05.2f',DecSign(I),Dec(I,2:4));
        otherwise
           error('Unknown  OutType Option');
       end
       if (IsCell==1)
          OutData{I}   = StrDec;
       else 
          OutData(I,:) = StrDec;
       end
    end

 otherwise
    error('Illegal OutType');
end

    
 
