This package uses LibRaw and digiCamControl (to be added) to control DSLR cameras and to open the RAW files. 

RAW files
-------------------
This is based on a windows 64 bit pre-compiled LibRaw instance. For other OS compile your own. 
The mex file open_raw.mexw64 is for windows 64. 
The function getImage takes a filename and returns an image (MxNx3 16 bit). 
The resulting image has correct RGB order but might be very dark and green tinted. 
Multiply the image by ~16 to convert from 14 bit to 16 bit brightness. 
The green channel is averaged between the two green pixel in the Bayer array. 

Camera control
------------------
to be added. 