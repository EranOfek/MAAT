function [List]=find_urls_ftp(URL,varargin)
% Find files in a FTP link
% Package: www
% Description: Make a list of files in an FTP link.
% Input  : - FTP URL.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Match' - Regular expression for file matching. If empty, then
%                      select all files. Default is [].
%            'Verbose' - Verbose. Default is true.
% Output : - Cell array of links names.
%          - Cell array of file names.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jan 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [List]=www.find_urls_ftp('ftp://archive.noao.edu/public/hlsp/nscdr1/instcal/','Match','fits$');
% Reliable: 2
%--------------------------------------------------------------------------

DefV.Match               = [];
DefV.Verbose             = true;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% Str = urlread(URL);
% C   = textscan(Str,'%s %d %d %d %d %s %d %s %s\n');
% 
% FileNames = C{end};
% Nf        = numel(FileNames);
% List      = cell(Nf,1);
% for If=1:1:Nf
%     List{If} = sprintf('%s%s',URL,FileNames{If});
% end



try
    Str = urlread(URL);
catch
    fprintf('urlread failed - try again in 10min');
    pause(600);
    try
        Str = urlread(URL);
    catch
        fprintf('urlread failed - try again in 60min');
        fprintf('save List.mat');
        save -v7.3 List.mat List
        pause(3600);
        
        Str = urlread(URL);
    end
    
end
    
C   = textscan(Str,'%s %d %d %d %d %s %d %s %s\n');
FileNames = C{end};
DataType = C{1};
Nf        = numel(FileNames);

List = Util.struct.struct_def({'FileName','Link','FullLink'},0,1);

if (InPar.Verbose)
    fprintf('URL = %s\n',URL)
end

Ind = 0;
for If=1:1:Nf
    if (InPar.Verbose)
        fprintf('File number %d out of %d\n',If,Nf);
    end
    
    switch lower(DataType{If}(1))
        case 'd'
            % directory
            % recursive call
            
            RecURL = sprintf('%s/%s/',URL,FileNames{If});
            [RecList] = www.find_urls_ftp(RecURL,'Match',InPar.Match,'Verbose',InPar.Verbose);
            
            List = [List; RecList];
            
        otherwise
            % a file
            
            Found = true;
            if (~isempty(InPar.Match))
                if (isempty(regexp(FileNames{If},InPar.Match,'match')))
                    Found = false;
                end
            end
            
            if (Found)
                Ind = numel(List) + 1;
                List(Ind).FileName   = FileNames{If};
                List(Ind).Link       = URL;
                List(Ind).FullLink =  sprintf('%s%s', List(Ind).Link, List(Ind).FileName);
            end
    end
end
