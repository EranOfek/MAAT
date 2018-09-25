%diagnostic tool for checking mean value of each channel to ensure that
%channels are ordered as expected. Note: values are dependent on white
%balance parameters set in open_raw.cpp

function [mean1,mean2,mean3,mean4] = meanCheck(filename)
    Im = dslr.getImageRGBG(filename);
    mean1 = mean(mean(Im(:,:,1)));
    mean2 = mean(mean(Im(:,:,2)));
    mean3 = mean(mean(Im(:,:,3)));
    mean4 = mean(mean(Im(:,:,4)));
end