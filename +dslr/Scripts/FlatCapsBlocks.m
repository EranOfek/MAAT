%This script automatically captures flats at different exposure times at fixed
%illumination. It constructs median and variance images from stacked images
%at fixed exposure time. Then, the median and variance images are cropped
%then divided into sections whose sizes are determined by X/YchopDim. The
%median of the sections is calulated and stored. The script outputs
%medianStore and varStore arrays of dimensions m x n x 4 x num of exp times
%
%NOTE: at the moment, crop dimensions/final image dimention selection are hard-coded. REMEMBER TO CHANGE
%THEM FOR DIFFERENT IMAGE SIZES


C = CameraController;                        %Initialize digiCam interface, sets fixed exposure params
C.session.deletefileaftertransfer = 1;
C.camera.fnumber = 22;
C.camera.isonumber = '800';

%------SPECIFY THESE VALUES-------------------------
exptimes = [1 2 3 4 5];     %vector of exposure times in seconds. Be mindful to specify exposure times accepted by camera    
waitTime = 5;               %specify time to waitbetween exposures (Minimum for D700: 5)
cycles = 4;                 %number of images at fixed exposure times to be stacked
XchopDim = 20;              %x,y dims of sections divided to be averaged over in median and variance image
YchopDim = 20;              
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

disp('Begin Capture')

%This is far from optimal - ideally, create array of proper size to begin with.
%Doing that led me to some strange bugs, and this is a somewhat temporary
%fix for them
              
rawStack = [];       
medianStore = [];
varStore = [];

for i = exptimes    
    C.camera.shutterspeed = strcat(num2str(i),'s');        %set camera exposure time
    
    for k = 1:cycles
        C.Capture;                          %capture image
        pause(i+waitTime);                  %wait long enough for files to transfer
    end

    %NOTE: this script retrieves file after cycle is over. This changes in
    %later scripts
    
    lastfilename = C.lastfile;                              %retreive file name of image just captured
    lastfilenumber = sscanf(lastfilename(5:8),'%f');        %retrieve file number from string
    firstfilenumber = lastfilenumber - cycles +1;           %obtain first file number


    for filenumber = firstfilenumber:lastfilenumber                     %ranges over file numbers
        filename = strcat('DSC_',num2str(filenumber,'%04.f'),'.nef');    %construct individual filenames from filenumbers

        count = 0;
        while exist(filename,'file')==0                 %will wait for file to appear in search time for 90s
            pause(1);
            count = count+1;
            
            if count == 60                              %max time allowed for file to appear
            disp(strcat(filename,' not found.'));
            break
            end
        end     

        rawImage = dslr.getImageRGBG(filename);         %load raw image with four channels in RGBG order
        rawStack = cat(4,rawStack,rawImage);            %stack images along the fourth dimension (instead of dim 3, this preserves channel individuality
    end
    
    medianIm = median(double(rawStack),4);              %calculate median and variance images through stack
    varIm = var(double(rawStack),0,4);                  

    rawStack = [];                                      %purge rawStack

    medianImCrop = medianIm(11:1410,2:2141,:);      %Crop images arbitrarily FUTURE: automate cropping
    varImCrop = varIm(11:1410,2:2141,:);            %**CHANGE HERE**

    medianChop = [];            
    varChop = [];

    for Xchop = 1:2140/XchopDim                 %**CHANGE HERE** FUTURE: automate total image dimension selection
        XchopVec = ((Xchop-1)*XchopDim+1):(Xchop*XchopDim); %create a vector that gives Xchop coordinates for indexing later

        for Ychop = 1:1400/YchopDim
            YchopVec = ((Ychop-1)*YchopDim+1):(Ychop*YchopDim); %similarly as above for Ychop

            for CH = 1:4
                medianChop(Ychop,Xchop,CH) = median(median(medianImCrop(YchopVec,XchopVec,CH)));    %create new median and variance images with "chopped" squares
                varChop(Ychop,Xchop,CH) = median(median(varImCrop(YchopVec,XchopVec,CH)));          %this is done chanel-by-channel     
            end                                                  
        end
    end

    medianStore = cat(4,medianStore,medianChop);    %store median and variance chopped images as array where dim4 corresponds to exposure time number
    varStore = cat(4,varStore,varChop); 
end


disp('Capture and Analysis Finished.')

clear Xchop             %clean up
clear Ychop
clear nYchop
clear nXchop
clear waitTime
clear k
clear n
clear i
clear count
clear CH
clear varCH
clear meanCH
