function ZI=interp2fast(VecX,VecY,Z,XI,YI,varargin)
% Faster version of interp2
% Package: Util.interp
% Description: A faster version of the interp2.m built in function.
%              This function is faster than interp2.m when the XI and YI
%              vectors span over a small range in the VecX and VecY
%              space.
% Input  : - Vector defining the X coordinates in the interpolated matrix.
%          - Vector defining the Y coordinates in the interpolated matrix.
%          - Matrix to interpolate.
%          - X values in which to interoplate the matrix.
%          - Y values in which to interpolate the matrix.
%          * Arbirtrary number of input arguments to pass to interp2.m
%            (e.g., interpolation method). Otherwise will use the interp2
%            defaults.
% Output : - The interoplated vector.
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Aug 2013
%    UR: : http://weizmann.ac.il/home/eofek/matlab/
% Example: InterpTrace = Util.interp.interp2fast(VecX,VecY,ScienceImage,XI,YI,'linear');
% Reliable: 2
%-------------------------------------------------------------------------
Pad = 2;

MinXI = min(XI(:)); % check boundries of XI, YI
MaxXI = max(XI(:));
MinYI = min(YI(:));
MaxYI = max(YI(:));


StartX = VecX(find((VecX-MinXI) <0,1,'last'))-Pad;
if (isempty(StartX))
    StartX = 1;
end
if (StartX<1)
    StartX = 1;
end

EndX   = VecX(find((VecX-MaxXI) >0,1,'first'))+Pad;
if (isempty(EndX))
    EndX = max(VecX);
end
if (EndX>max(VecX))
    EndX = max(VecX);
end

StartY = VecY(find((VecY-MinYI) <0,1,'last'))-Pad;
if (isempty(StartY))
    StartY = 1;
end
if (StartY<1)
    StartY = 1;
end

EndY   = VecY(find((VecY-MaxYI) >0,1,'first'))+Pad;
if (isempty(EndY))
    EndY = max(VecY);
end
if (EndY>max(VecY))
    EndY = max(VecY);
end

% regular (slow version) example
%tic;
%ZI = interp2(VecX,VecY,Z,XI,YI,varargin{:});
%toc

ZI = interp2(VecX(StartX:EndX),VecY(StartY:EndY),Z(StartY:EndY,StartX:EndX),XI,YI,varargin{:});
