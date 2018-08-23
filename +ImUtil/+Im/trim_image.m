function [Mat,Shift,Section,SectionF]=trim_image(Mat,Section,Method,PadVal)
% Trim/crop a section from a 2D matrix.
% Package: ImUtil.Im
% Description: Trim a section from a 2D matrix (image). If the trim region
%              is outside the image bounderies, then the trimmed
%              image can be optionally padded to have the size of the
%              requested trim region.
% Input  : - A 2D matrix (image).
%          - A trim section: [Xmin, Xmax, Ymin, Ymax]
%            or [Xcenter, Ycenter, Xhalf_width, Yhalf_width].
%          - Trim section representation method:
%            'section' - Trim section format is: [Xmin, Xmax, Ymin, Ymax].
%            'center'  - Format: [Xcenter, Ycenter, Xhalf_width, Yhalf_width].
%            Default is 'section'.
%          - Pad value for out of bounderies regions. If empty then
%            do not pad. Default is NaN.
% Output : - A trimmed image matrix.
%          - Shift (new-old) [X,Y] of new trimed image relative to
%            original image.
%          - The section in 'section' format.
%          - The truncted section.
% See also: SIM/trim_image.m, SIM/stamp_coo.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Mat=ImUtil.Im.trim_image(rand(3,3),[1 2 2 3])
%          ImUtil.Im.trim_image([1 2 3;4 5 6; 7 8 9],[0 4 3 4]);
%          ImUtil.Im.trim_image([1 2 3;4 5 6; 7 8 9],[0 4 3 4],'section',0);
% Reliable: 2
%--------------------------------------------------------------------------


Def.Method = 'section';
Def.PadVal = NaN;
if (nargin==2),
    Method = Def.Method;
    PadVal = Def.PadVal;
elseif (nargin==3),
    PadVal = Def.PadVal;
elseif (nargin==4),
    % do nothing
else
    errId = 'trim_image.m:TooManyInputArguments';
    errMsg = 'Mat, Section, [Method, Boundry]';
    error(errId, errMsg);
end
     
switch lower(Method)
    case 'section'
        % do nothing - already in section format [Xmin Xmax Ymin Ymax]
    case 'center'
        % convert [Xcenter, Ycenter, Xhalfwidth, Yhalfwidth] to [Xmin Xmax Ymin Ymax]
        Section = [Section(1)-Section(3), Section(1)+Section(3), Section(2)-Section(4), Section(2)+Section(4)];
    otherwise
        error('Unknown Method option');
end

[SizeY, SizeX] = size(Mat);
if (Section(2)<1 || Section(4)<1 || Section(1)>SizeX || Section(3)>SizeY),
    error('Requested trim position is completly outside image boundries');
end

% if Boundry is empty then truncate the Section to be within
% the image bounderies.
SectionF(1) = max(1,Section(1));     % Xmin
SectionF(3) = max(1,Section(3));     % Ymin
SectionF(2) = min(SizeX,Section(2)); % Xmax
SectionF(4) = min(SizeY,Section(4)); % Ymax

    
% trim image
Mat = Mat(SectionF(3):SectionF(4),SectionF(1):SectionF(2));

% pad image if needed
if (~isempty(PadVal)),
    % Boundry may be out of image - pad bounderies
    PadXlow    = max(0,1-Section(1));
    PadXhigh   = max(0,Section(2)-SizeX);
    PadYlow    = max(0,1-Section(3));
    PadYhigh   = max(0,Section(4)-SizeY);
    
    Mat = padarray(Mat,[PadYlow, PadXlow],  PadVal,'pre');
    Mat = padarray(Mat,[PadYhigh, PadXhigh],PadVal,'post');

    % shift: new-old
    Shift    = zeros(1,2);
    Shift(1) = Section(1) - 1;
    Shift(2) = Section(3) - 1;
    
else
    % Mat is already correct
    
    % shift: new-old
    if (nargout>1),
        Shift    = zeros(1,2);
        Shift(1) = max(1,Section(1) - 1);
        Shift(2) = max(1,Section(3) - 1);
    end
end



    
