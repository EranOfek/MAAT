function List=wget_lat_weekly_data(WeekNumber,Get,pwgetPar)
% Retrieve the Fermi/LAT weekly photon data in FITS format
% Package: VO.Fermi
% Description: Retrieve the last version of the Fermi/LAT weekly files
%              from the HEASARC FTP site in FITS format.
%              The function can get the list of files and their URLs,
%              retrieve all files or retrieve specific files.
% Input  : - A vector of files to retrieve.
%            If Inf then get all files. Default is Inf.
%          - Get files: true|false. Default is true.
%          - A cell array of additional arguments to pass to www.pwget.
%            Default is {}.
% Output : - A cell array of the links to retrive.
%     By : Eran O. Ofek                    Jan 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: List=VO.Fermi.wget_lat_weekly_data(Inf,false);
% Reliable: 2

if (nargin<3)
    pwgetPar = {};
    if (nargin<2)
        Get = true;
        if (nargin<1)
            WeekNumber = Inf;
        end
    end
end

BaseURL = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly/photon/';

% Get all FITS file links from Fermi HEASARC archive:
List    = www.find_urls(BaseURL,'match','https://heasarc.gsfc.nasa.gov/.+lat_photon.+fits');
Nfile   = numel(List);

% break file names by "_" and find the week number of each file
RegList = regexp(List,'_','split');

% week number is stored in field 4
WeekNumberInd = 4;

WeekList = zeros(Nfile,1);
for Ifile=1:1:Nfile
    WeekList(Ifile) = str2double(RegList{Ifile}{WeekNumberInd}(2:end));
end

% search for requested week numbers
if isinf(WeekNumber)
    % get all weeks
    Flag = true(Nfile,1);
else
    Flag = Util.array.findmany(WeekList,WeekNumber);
end

List = List(Flag);

if (Get)
    www.pwget(List,pwgetPar{:});
end
