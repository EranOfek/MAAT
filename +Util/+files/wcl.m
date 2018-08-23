function Counter=wcl(File)
% Count the number of lines in a file.
% Package: Util.files
% Description: Count the number of lines in a file.
% Input  : - String containing file name.
% Output : - Number of lines in file.
% Tested : Matlab 7.8
%     By : Eran O. Ofek                    Nov 2009
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Counter=wcl('wcl.m');
% Reliable: 2
%--------------------------------------------------------------------------
FID = fopen(File,'r');
Counter = 0;
while (feof(FID)==0),
   fgetl(FID);
   Counter = Counter + 1;
end
fclose(FID);

