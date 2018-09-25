%Script accepts inputs from FlatCapsBlocks to estimate channel-by-channel
%gain =, read noise, and dark current


%P = [fG fG^2 dcG^2 rnG^2];

G = [];

for CH = 1:4
    uH = [];
    lH = [];
    medianVecStore = [];
    varVecStore = [];
    
    for m = 1:size(exptimes,2)

        medianVec = reshape(medianStore(:,:,CH,m),[],1);        %reshape individual channel and exposure time median and var
        varVec = reshape(varStore(:,:,CH,m),[],1);              %into column vector

        nMedVar = size(medianVec,1);                            %number of median and var measurements

        exptime = exptimes(m);                          

        tHolder = ones(nMedVar,1)*exptime;        %creates column vector of constant exposure times of height equal to number of median/var measurements for exposure time
        zHolder = zeros(nMedVar,1);                 %creates column vector of zeroes of same dim as above
        oneHolder = ones(nMedVar,1);              % creates column vector of ones of same dim as "tHolder"

        uHpc = [];               
        uHpc = cat(2,tHolder,zHolder,zHolder,zHolder);      %"uHpc" is matrix corresponding to median values at fixed t
        uH = cat(1,uH,uHpc);                        %"uH" is matrix corresponding to all median values

        lHpc = [];
        lHpc = cat(2,zHolder,tHolder,tHolder,oneHolder);    %similarly, "lHpc" is matrix corresponding to var values at fixed t
        lH = cat(1,lH,lHpc);                        %similarly, "lH" is matrix corresponding to all var values

        medianVecStore = cat(1,medianVecStore,medianVec);       %concatenates successive median vectors to form Y
        varVecStore = cat(1,varVecStore,varVec);                %concatentaes succesive var vectors to form Y
    end

    H = cat(1,uH,lH);                           %form full H matrix
    Y = cat(1,medianVecStore,varVecStore);      %form full Y vector
 
    P = H\Y;                %find LSS of HP = Y

    G = cat(2,G,P(2)/P(1));  %vector of gain for each channel, if gain = fGsquare/fG
    DC = P(3)./((P(2)/P(1))^2);          
    RN = P(4)./((P(2)/P(1))^2);        %NOTE: this is actually RN^2

end