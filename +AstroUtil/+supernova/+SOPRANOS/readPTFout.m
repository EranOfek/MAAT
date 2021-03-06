function [tbl, RA, Dec, redshift, discJD] = readPTFout(SNname, startRow, endRow)
% read PTF out file  
% Package: AstroUtil.supernove.SOPRANOS
% Description: read PTF out file and read the matching GALEX data of the
%              GALEX/PTF experiment. 
% Input  : - the supernova name (for example 'PTF12gnt')
% Output : - Table with file contents
%               
% See also: AstroUtil.supernova.SOPRANOS.calcGrid
% Tested : Matlab 9.5
%     By : Noam Ganot                      Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% AstroUtil.supernova.SOPRANOS.readPTFout('PTF12gnt');
% Reliable: 2
%--------------------------------------------------------------------------

% Auto-generated by MATLAB on 2019/06/29 10:28:55

filename = sprintf('%s.out_PTF48R',SNname);

%% Initialize variables.
if nargin<=2
    startRow = 16;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: categorical (%C)
%	column6: categorical (%C)
%   column7: categorical (%C)
%	column8: categorical (%C)
%   column9: categorical (%C)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%11f%11f%10f%7f%5C%3C%4C%7C%C%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r','n','UTF-8');
% Skip the BOM (Byte Order Mark).
fseek(fileID, 3, 'bof');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

% %% Create output variable
% tbl = table(dataArray{1:end-1}, 'VariableNames', {'MJD','counts','dcounts','zp', 'dzp','magsys','Telescope','Detector','PhotSource'});

%% Add filter instr fields
mAB = zeros(size(dataArray{2}));
mAB(dataArray{2}>0)= -2.5 * log10(dataArray{2}(dataArray{2}>0)) + dataArray{4}(dataArray{2}>0);
mAB(dataArray{2}<=0) = 99;
dataArray{10}= mAB;
mErr = zeros(size(dataArray{2}));
mErr(dataArray{2}>0)= 2.5/log(10) * dataArray{3}(dataArray{2}>0)./dataArray{2}(dataArray{2}>0);
mErr(dataArray{2}<=0) = 99;
dataArray{11}= mErr;
% convert ZP to flux to get PTF flux factor:
FF = convert.flux(dataArray{4},'AB','cgs/A',AstFilter.get('PTF','r').pivot_wl_photon,'A');
dataArray{12}= dataArray{2}.*FF;
dataArray{13}= dataArray{3}.*FF;
dataArray{14}=categorical(ones(size(dataArray{1})),1,{'P48+PTF'});
dataArray{15}=categorical(ones(size(dataArray{1})),1,{'r_p48'});

%% read GALEX data from InfoAll.mat
[ ~, MJD, OBS, ERR, f, fERR, ~, ~, RA, Dec, DiscMJD, redshift, mInd, mAB,  mERR ] = AstroUtil.supernova.SOPRANOS.load_GALEX( SNname );
MagAB    = zeros(size(OBS));
MagErrAB = zeros(size(OBS));
MagAB(mInd) = mAB; MagAB(~mInd) = nan;
MagErrAB(mInd) = mERR; MagErrAB(~mInd) = nan;

% add GALEX entries to dataArray
dataArray{1}=[dataArray{1};MJD];
dataArray{2}=[dataArray{2};OBS];
dataArray{3}=[dataArray{3};ERR];
dataArray{4}=[dataArray{4};ones(size(MJD))*AstFilter.get('GALEX','NUV').pivot_wl_photon];
dataArray{5}=[dataArray{5};categorical(ones(size(MJD)),1,{'None'})];
dataArray{6}=[dataArray{6};categorical(ones(size(MJD)),1,{'AB'})];
dataArray{7}=[dataArray{7};categorical(ones(size(MJD)),1,{'GALEX'})];
dataArray{8}=[dataArray{8};categorical(ones(size(MJD)),1,{'NUV'})];
dataArray{9}=[dataArray{9};categorical(ones(size(MJD)),1,{'AperPhot'})];
dataArray{10}=[dataArray{10}; MagAB];
dataArray{11}=[dataArray{11}; MagErrAB];
dataArray{12}=[dataArray{12}; f];
dataArray{13}=[dataArray{13}; fERR];
dataArray{14}=[dataArray{14};categorical(ones(size(MJD)),1,{'GALEX'})];
dataArray{15}=[dataArray{15};categorical(ones(size(MJD)),1,{'NUV'})];


%% Create unified output variable
tbl = table(dataArray{:}, 'VariableNames', {'MJD','counts','dcounts','zp','dzp','magsys','Telescope','Detector','PhotSource','mag','magerr','flux','fluxerr','instr','filter'});



