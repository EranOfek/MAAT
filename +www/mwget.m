function [Status]=mwget(URL,varargin)
% A wrapper around the wget command
% Package: www
% Description: A wrapper around the wget command. Retrieve a URL using
%              the wget command. 
%              OBSOLETE: Use www.pwget instead.
% Input  : - A string containing a URL.
%            Alternativly this can be a cell array of URLs.
%            However, if this is a cell array than 'rftp' will be set to 0.
%          * Arbitary number of ...,keyword,value,... pairs.
%            Where allowed keywords are:
%            'user' - user
%            'pass' - password
%            'clob' - Clobber {'y'|'n'}. Default is 'n'.
%                     No clobber means that if the file exist it will
%                     be overwritten.
%            'par'  - string of additional parameters to pass to wget.
%            'out'  - Output file name. If the first inpur parameter
%                     is a cell array this must be a cell array too.
% Output : - Vector of flags [0 | 1] indicating if the program was
%            able to retrive the requested file (0) or failed (1).
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    May 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: www.pwget.m, www.rftpget.m
% Example: [Res]=www.mwget('http://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/2008_05/00037562001/xrt/event/sw00037562001xpcw2po_cl.evt.gz');
% Reliable: 2
%----------------------------------------------------------------------------
IsFTP = ~isempty(strfind(URL,'ftp://'));



DefV.Par     = '';
DefV.Clob    = 'n';
DefV.User    = NaN;
DefV.Pass    = NaN;
DefV.Out     = NaN;

InPar = set_varargin_keyval(DefV,'y','def',varargin{:});

if (~iscell(URL)),
   URL = {URL};
end

AddArg = '';
if (IsFTP),
   AddFTP = 'ftp-';
else
   AddFTP = '';
end

if (~isnan(InPar.User)),
   AddArg = sprintf('%s--%suser=%s ',AddArg,AddFTP,InPar.User);
   AddArg = sprintf('%s--%supassword=%s ',AddArg,AddFTP,InPar.Pass);
end

switch lower(InPar.Clob)
 case 'n'
    AddArg = sprintf('%s -nc ',AddArg);
 otherwise
    % do nothing
end

if (~isempty(InPar.Par)),
   AddArg = sprintf('%s %s ',AddArg,InPar.Par);
end
if (~isnan(InPar.Out)),
   AddArg = sprintf('%s--output-document=%s ',AddArg,InPar.Out);
end



Nurl = length(URL);
for Iurl=1:1:Nurl,   
   Command = sprintf('wget %s%s',AddArg,URL{Iurl});
   Status(Iurl) = system(Command);
end

