%This script runs FlatCapsBlocks to extimate gain, DC, and read noise
%using designMatrixEstimation. Outputs Gstore, RNstore, and DCstore of
%dimensions number of runs x 4. Gain for each channel is calculated
%separately


%------SPECIFY VALUE------
runs = 5;
%^^^^^^^^^^^^^^^^^^^^^^^^^


Gstore = [];
RNstore = [];
DCstore = [];

for h = 1:runs
    
    FlatCapsBlocks;           %call FlatCapsBlocks (make sure parameters are properly set within script
    
    designMatrixEstimation;     %run design matrix 
    
    Gstore = cat(1,Gstore,G);       %store gain, read noise, and DC
    
    RNact = sqrt(RN);
    
    RNstore = cat(1,RNstore,RNact);
    
    DCstore = cat(1,DCstore,DC);
 
    clear H                 %clean up from designMatrixEstmation
    clear Y
    clear P
end

clear runs          %clean up
clear h


