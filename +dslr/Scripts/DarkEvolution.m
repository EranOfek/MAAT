%This script will measure the evolution of Dark Current with time over long
%exposure periods due to heating of the sensor. Exposure time, number of
%shots, and waitTime must be specified within the script as variables.
%After every 10 exposures, standard deviation and mean are calculated pixel-by-pixel
%through the stack. Only the median value of these for each channel is
%stored.

C = CameraController;           %Initialize digiCam interface, sets fixed exposure params
C.session.deletefileaftertransfer = 1;
C.camera.fnumber = 22;
C.camera.isonumber = '200';


%-----SPECIFY THESE VALUES --------------
exptime = 10;       %exposure time in seconds
shots = 2000;        %number of exposures
waitTime = 5;       %time to wait between exposures. For D700, recommended is 5
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


C.camera.shutterspeed = strcat(num2str(exptime),'s');   %sets camera shutter speed

%once again this method won't be ideal in terms of memory allocation. To be
%fixed in later iterations

rawStack = [];
elapsedTime = [];
medianSTD = [];
medianMean = [];

tic;

for i = 1:shots 
    C.Capture;                             %capture image
    pause(exptime+waitTime);               %wait long enough for files to transfer
    
    lastfilename = C.lastfile;             %retrieve file name of image just captured   
     
    count = 0;
    while exist(lastfilename,'file')==0    %will wait for file to appear in search time for 90s
            pause(1);
            count = count+1;
        
            if count == 90                 %max time allowed for file to appear
            disp(strcat(lastfilename,' not found.'));
            break
        end
    end     

    rawImage = dslr.getImageRGBG(lastfilename);         %load raw image with four channels in RGBG order               
    rawStack = cat(4,rawStack,rawImage);                %stack images along the fourth dimension (instead of dim 3, this preserves channel individuality
                    
    elapsedTime = cat(2,elapsedTime,toc);               %logs current elapsed time
    
    if rem(i,10) == 0                                   %every 10 images, std and mean images will be calculated from the 10 stacked images
        stdIm = std(double(rawStack),0,4);
        meanIm = mean(rawStack,4);
        rawStack = [];                                  %last 10 stacked images are purged to free up memory
        
        medianSTD = cat(2,medianSTD,median(median(stdIm,1),2));     %only median of STD and mean images are stored
        medianMean = cat(2,medianMean,median(median(meanIm,1),2));         
    end
end

medianSTD = squeeze(medianSTD);     %singleton dimensions are removed from median STD and mean arrays
medianMean = squeeze(medianMean);   %resulting arrays are of dimension shots/10 x 4

Inc = 10;                               %number of images per increment - leave at 10
elapsedTimeStack = [];                  
nInc = size(elapsedTime,2)/Inc;         %calculate number of increments

for a = 1:nInc                          %compresses elapsedTime vector for number of entries to correspond to those in medianmean and median STD
    elapsedTimeTemp = mean(elapsedTime(1,(a-1)*Inc+1:a*Inc));       
    elapsedTimeStack = cat(2,elapsedTimeStack,elapsedTimeTemp);
end


clear rawStack              %clean up
clear i
clear stdIm
clear meanIm
clear rawImage
clear count
clear lastfilename
clear shots
clear waitTime
clear expTime
clear elapsedTimeTemp
clear elapsedTime