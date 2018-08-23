function [Struct,List]=sedm_prep_nightlog(Images,OutputFile)
%--------------------------------------------------------------------------
% sedm_prep_nightlog function                                         SEDM
% Description: Given a list of SEDM images, prepare a log file of images
%              with basic information derived from the image headers.
% Input  : - List of images (see create_list.m for options).
%            Default is '*.fits'.
%          - String containing log output file name. Default is 'log'.
% Output : - Structure array containing the requested keyword values per
%            image.
%          - Cell array of image names.
% Tested : Matlab R2013a
%     By : Eran O. Ofek                    Aug 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: sedm_prep_nightlog;
% Reliable: 2
%--------------------------------------------------------------------------

Def.Images     = '*.fits';
Def.OutputFile = 'log';
if (nargin==0),
    Images     = Def.Images;
    OutputFile = Def.OutputFile;
elseif (nargin==1),
    OutputFile = Def.OutputFile;
else
    % do nothing
end
    
[~,List] = create_list(Images,NaN);

KeysToGet = {'OBJECT','EXPTIME','JD','RA','DEC','AIRMASS','UTC'};
[~,~,Struct,List]=fits_mget_keys(List,KeysToGet);

FID = fopen(OutputFile,'w');
fprintf(FID,'%%%19s %35s %7s %12s %14s %14s %7s %20s\n','Image',KeysToGet{:});
for Iim=1:1:numel(Struct),
   fprintf(FID,'%20s %35s %7.1f %12.4f %14s %14s %7.4f %20s\n',List{Iim},...
                                                         Struct(Iim).(KeysToGet{1}),...
                                                         Struct(Iim).(KeysToGet{2}),...
                                                         Struct(Iim).(KeysToGet{3}),...
                                                         Struct(Iim).(KeysToGet{4}),...
                                                         Struct(Iim).(KeysToGet{5}),...
                                                         Struct(Iim).(KeysToGet{6}),...
                                                         Struct(Iim).(KeysToGet{7}));
                                                     
end
fclose(FID);


