function image2avi(List,OutputName,FPS,Map)
%------------------------------------------------------------------------------
% image2avi function                                                   General
% Description: Create avi file from list of images
% Input  : - File name containing list of images (e.g., in jpg).
%            Alternatively, this can be a cube in which the third dimesion
%            corresponds to the image index (time).
%          - Output avi file name.
%          - frame per sec, default is 25.
%          - Color map parameter to pass to imshow.m.
%            Default is [0 1].
% Example : image2avi('image.list','out.avi',5);
% Tested : Matlab 6.1
%     By : Eran O. Ofek                       May 2002
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Reliable: 2
%------------------------------------------------------------------------------

if (nargin==2),
   FPS = 25;
   Map = [0 1];
elseif (nargin==3),
   Map = [0 1];
elseif (nargin==4),
   % do nothing
else
   error('Illigal number of input arguments');
end



if (ischar(List)),
   FID = fopen(List,'r');
end

%fig=figure;
%set(fig,'DoubleBuffer','on');
%MovH = avifile(OutputName,'compression','none','fps',FPS);
MovH = VideoWriter(OutputName);
open(MovH);

%N = wc(List,'l');
N = size(List,3);
for I=1:1:130,
%while (feof(FID)~=1),

    if (ischar(List)),
       ImageName = fgetl(FID);
       Image     = imread(ImageName);
    else
       Image     = squeeze(List(:,:,I));
    end
%    H = imshow(uint8(Image),Map);
    H = imshow(Image,Map);

    CurrFrame = getframe;
    writeVideo(MovH,CurrFrame);
    %%set(H,'EraseMode','xor');
    %F = getframe(gca);
    %MovH = addframe(MovH,F);
end
close(MovH);


if (ischar(List)),
   fclose(FID);
end

