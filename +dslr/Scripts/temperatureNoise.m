%selects reference images for cold and hot DC at different ISO values
%use with the following syntax:
%
%temperatureNoise(coldLoFile,hotLoFile,coldHiFile,hotHiFile,fact) 
%
%where cold/hot correspond to camera temperature, lo/hi correspond to hi and
%low ISO, and fact corresponds to a multiplicative factor for easier image
%viewing


function temperatureNoise(coldLoFile,hotLoFile,coldHiFile,hotHiFile,fact)

coldLo = dslr.getImageRGB(coldLoFile);
hotLo = dslr.getImageRGB(hotLoFile);
coldHi = dslr.getImageRGB(coldHiFile);
hotHi = dslr.getImageRGB(hotHiFile);

coldHiCrop = coldHi(501:800,801:1100,:)*fact;   %crop all images
coldLoCrop = coldLo(501:800,801:1100,:)*fact;
hotLoCrop = hotLo(501:800,801:1100,:)*fact;
hotHiCrop = hotHi(501:800,801:1100,:)*fact;

%display bright cropped images side-by-side
figure
subplot(2,2,4)
image(hotHiCrop)
title('High Temp. ISO 800')
subplot(2,2,2)
image(hotLoCrop)
title('High Temp. ISO 200')
subplot(2,2,3)
image(coldHiCrop)
title('Low Temp. ISO 800')
subplot(2,2,1)
image(coldLoCrop)
title('Low Temp. ISO 200')


clear fact      %clean up