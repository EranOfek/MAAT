function AS=spec_read_ascii(File,varargin)
%--------------------------------------------------------------------------
% spec_read_ascii function                                       AstroSpec
% Description: Read a spectrum from an ascii/text file.
% Input  : - File name, file name with wild cards, or a cell array of
%            file names. Each file is a text file in the current directory
%            containing spectra in mtarix, cell or AstSpec formats.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'SpecOutType' - Output type - options are:
%                       'AstSpec' - AstSpec class object (default).
%                       'cell'    - Cell array of spectra.
%                       'mat'     - Matrix of a single spectrum.
%            'Delimiter' - Delimitar. Default is ' '.
%            'Comment'     - Comment style (e.g., '#'). If empty, then
%                       will attempt to guess the commnet type by reading
%                       the first character in the file. If it is
%                       '#','%' or '/' then set this to be the comment
%                       type, otherwise set the comment type to '%'.
%                       Default is empty.
%            'HeaderLines' - Number of header lines to skip. Default is 0.
%            'Ncol'    - Number of columns to read. Default is 2.
%            'ColNames'- Column names (should be consistent with the
%                        AstSpec class definition).
%                        Default is {'Wave','Int','Err','Back','Mask'}.
%            'ColCorr' - Multiplication correction factor to apply to each
%                        column. Default is [1     ,1    ,1    ,1     ,1 ].
%            'WaveUnits' - Wavelength units. Default is 'Ang'.
%            'IntUnits'  - Intensity units. Default is
%                          'erg*cm^{-2}*s^{-1}*Ang^{-1}'.
% Output : - The output spectra.
% See also: 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AS=spec_read_ascii('*.dat');
% Reliable: 2
%--------------------------------------------------------------------------

DefV.SpecOutType      = 'AstSpec';   % 'AstSpec','cell','mat'
DefV.Delimiter        = ' ';
DefV.Comment          = [];
DefV.HeaderLines      = 0;
DefV.Ncol             = 2;
DefV.ColNames         = {'Wave','Int','Err','Back','Mask'};
DefV.ColCorr          = [1     ,1    ,1    ,1     ,1     ];
DefV.WaveUnits        = 'Ang';
DefV.IntUnits         = 'erg*cm^{-2}*s^{-1}*Ang^{-1}';

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

[~,List] = Util.files.create_list(File);
Nl       = numel(List);

ColNames = InPar.ColNames(1:1:InPar.Ncol);
Ncol     = numel(ColNames);
Format   = sprintf('%s %%*[^\\n]',str_duplicate('%f ',Ncol));

% initialization
switch lower(InPar.SpecOutType)
    case 'astspec'
        AS = AstSpec(Nl,1);   % define AstSpec class
    case 'cell'
        AS = cell(Nl,1);
    case 'mat'
        if (Nl>1)
            error('mat output is possible only when reading a single file');
        end
    otherwise
        error('Unknown SpecOutType option');
end


for Il=1:1:Nl
    % for each file
    
    if (isempty(InPar.Comment))
        % try to guess comment from first character in file
        FID   = fopen(List{Il},'r');
        Char1 = fgetl(FID);
        fclose(FID);
        switch lower(Char1)
            case '#'
                InPar.Comment = '#';
            case '%'
                InPar.Comment = '%';
            case '/'
                InPar.Comment = '/';
            otherwise
                % use default of matlab
                InPar.Comment = '%';
        end
    end
                
    % read ascii
    FID = fopen(List{Il},'r');
    C   = textscan(FID,Format,...
                   'Delimiter',InPar.Delimiter,...
                   'CommentStyle',InPar.Comment,...
                   'Headerlines',InPar.HeaderLines);
    fclose(FID);

    switch lower(InPar.SpecOutType)
        case 'astspec'
            for Icol=1:1:Ncol,
                % for each column
                AS(Il).(ColNames{Icol}) = C{Icol}.*InPar.ColCorr(Icol);
            end
            
            % populate the AstSpec class
            AS(Il).WaveUnits = InPar.WaveUnits;
            AS(Il).IntUnits  = InPar.IntUnits;
            AS(Il).source    = 'spec_read_ascii.m';
            AS(Il).FileName  = List{Il};
        case 'cell'
            AS{Il} = cell2mat(C);
        case 'mat'
            AS     = cell2mat(C);
        otherwise
            error('Unknown SpecOutType option');
    end    
    
end
