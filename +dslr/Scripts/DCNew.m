%Revised method of plotting Dark Current vs. exposure time (using method
%from DarkEvolution). Takes a number of images at fixed exposure time and
%calculates STD and mean through stacked image, whose median is then
%stored. As before, exposure time, number of exposures per fixed time,
%wait time between exposures, and cooling time must be specified in script.

%-----SPECIFY THESE VALUES-------------------------------------------------
exptimes = [1 5 10 15 20 25 30];        %vector of exposure times in seconds. Be mindful to specify exposure times accepted by camera    
cycles = 8;                                %number of exposures at fixed exposure time
waitTime = 5;                               %time to wait between exposures. For D700, recommended is 5
coolTime = 120;                              %time to cool between cycles
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

C = CameraController;                    %Initialize digiCam interface, sets fixed exposure params
C.session.deletefileaftertransfer = 1;
C.camera.fnumber = 22;
C.camera.isonumber = '200';

%once again, not ideal
rawStack = [];
medianSTDStore = [];
medianMeanStore = [];

tic;

for m = exptimes
    C.camera.shutterspeed = strcat(num2str(m),'s');         %sets camera shutter speed
    
    for k = 1:cycles
        C.Capture;                          %capture image
        pause(m+waitTime);                  %wait long enough for files to transfer
    
        lastfilename = C.lastfile;          %retrieve file name of image just captured
     
        count = 0;
        while exist(lastfilename,'file')==0 %will wait for file to appear in search time for 90s
            pause(1);
            count = count+1;
            
            if count == 90                  %max time allowed for file to appear
            disp(strcat(lastfilename,' not found.'));
            break
            end
        end     

        rawImage = dslr.getImageRGBG(lastfilename);     %load raw image with four channels in RGBG order
        rawStack = cat(4,rawStack,rawImage);            %stack images along the fourth dimension (instead of dim 3, this preserves channel individuality        
                      
    end
    
    stdIm = std(double(rawStack),0,4);              %STD and mean images are calculated through the fixed exposure stack
    meanIm = mean(rawStack,4);
    rawStack = [];                                  %purge stack of images

    medianSTD = median(median(stdIm,1),2);         %store only median values of mean and STD so std/meanIm can be reused 
    medianMean = median(median(meanIm,1),2);      

    medianSTDStore = cat(1,medianSTDStore,medianSTD);      %store median STD and mean values     
    medianMeanStore = cat(1,medianMeanStore,medianMean);     

    pause(coolTime);            %pause for coolTime to allow camera to cool
end

toc                                 %output elapsed time

clear rawStack                      %clean up
clear waitTime
clear coolTime
clear m
clear k
clear medianSTD
clear medianMean
clear stdIm
clear meanIm
clear lastfilename
clear rawImage
clear cycles
clear count
