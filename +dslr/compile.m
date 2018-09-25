currentFold = cd;

cd +dslr

mex CXXFLAGS="$CXXFLAGS -std=c++11 -static" -Llibraw-0.19.0/lib -lraw -Ilibraw-0.19.0/libraw open_raw.cpp

cd(currentFold);

clear currentFold
clear ans