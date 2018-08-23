function [Cat,GoodTimes]=read_chandra_acis(ObsID,varargin)
% Read Chandra ACIS event files associated with a Chandra ObsID.
% Package: AstroX
% Description: Read Chandra ACIS event files associated with a Chandra
%              ObsID and add columns containing the RA and Dec for
%              each event. If needed cd to osbid directory.
%              OBSOLETE.
% Input  : - Observation ID (numeric) or directory (string) in which to
%            look for event file. If empty, assume the evt file is in
%            the current directory. Default is empty.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'WgetPar' - Cell array of additional parameters to pass to
%                        VO.Chandra.wget_obsid. Default is {}.
% Output : - AstCat object of ACIS events.
%          - Matrix of [start end] of good times.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Cat=AstroX.read_chandra_acis(ObsID);
%          Cat=AstroX.read_chandra_acis('/raid/eran/projects/Algorithms/PoissonMatchedFilter/Chandra/ObsID366/');
% Reliable: 2
%--------------------------------------------------------------------------



if (nargin==0)
    ObsID = [];
end

DefV.WgetPar             = {};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

PWD = pwd;
if (ischar(ObsID))
    % cd to data directory
    cd(ObsID);
else
    % get data
    if (isempty(ObsID))
        % assume data is in current directory
    else
        VO.Chandra.wget_obsid(ObsID,InPar.WgetPar{:});
    end
end

if (isdir('./primary'))
    cd('primary');
end
% Look for evt file
FileEvt = dir('*_evt2.fits*');
Nf      = numel(FileEvt);
if (Nf==0)
    %error('More than one evt2 file');
    Cat = [];
else
        
    % ungzip if needed
    if (any(strcmp(FileEvt(1).name(end-2:end),'.gz')))
        system_list('gzip -d %s',{FileEvt.name});
    end
    FileEvt = dir('*_evt2.fits');
    Nf      = numel(FileEvt);

    for If=1:1:Nf

        % read evt2 file
        Head  = FITS.get_head(FileEvt(If).name,1);
        InstVal = getkey(Head,'INSTRUME');
        switch lower(InstVal{1})
            case 'acis'
                
                Cat=FITS.read_table(FileEvt(If).name);

                [RA,Dec]=AstroX.xy2coo(FileEvt(If).name,Cat.Cat(:,Cat.Col.x),Cat.Cat(:,Cat.Col.y),'chandra');
                Cat.Cat = [Cat.Cat, RA, Dec];
                Cat.ColCell(end+1:end+2) = {'RA','Dec'};
                Cat     = colcell2col(Cat);
                
                TableGT     = fitsread(FileEvt(If).name,'BinTable',2);
                GoodTimes   = [TableGT{1}, TableGT{2}];
                
            otherwise
                Cat = [];
        end
    end
end

% cd to original directory
cd(PWD);
