function AS=spec_read(File,varargin)
%--------------------------------------------------------------------------
% spec_read function                                             AstroSpec
% Description: Read astronomical spectra from ASCII, FITS, HDF5, MAT files,
%              or matrix into a AstSpec class object.
% Input  : - The spectra to read, in one of the following formats:
%            1. A two or three column matrix [Wave, Int, [Error]].
%            2. A cell array of matrices (like in option 1).
%            3. String containing a file name, or wildcards for multiple
%               file names, and in which each file is either an
%               MAT file, ASCII file, FITS file, FITS table or HDF5 file.
%            4. An AstSpec class object (in this case the output will be
%               identical to the input).
%            All the files should have the same format.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'InType' - Input file type option. Options are:
%                       []        - try to identify file type automatically
%                                   (default).
%                       'astspec' - AstSpec class object.
%                       'ascii'   - text files containing spectra.
%                                   See spec_read_ascii.m for options.
%                       'mat'     - mat files containing spectra
%                                   See spec_read_mat.m for options.

%>>> NOT IMPLEMENTED
%                       'fits'    - FITS files containing spectra
%                                   See spec_read_fits.m for options.
%                       'hdf5'    - HDF5 files containing spectra
%                                   See spec_read_hdf5.m for options.
%>>>
%                       'numeric' - A matrix containing spectrum.
%                                   See AstSpec.mat2spec for options.
%                       'numeric-cell'-A cell array of matrices containing 
%                                   spectra.
%                                   See AstSpec.mat2spec for options.
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
%            spec_read_ascii.m, spec_read_mat.m
% Output : - An AstSpec class object that contains all the spectra.
%            Use astspec2mat.m to convert to matrices.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AS=spec_read('*.mat');
%          AS=spec_read(AS(1:10));
% Reliable: 2
%--------------------------------------------------------------------------


DefV.InType              = []; % empty -> automatic identification

if (~isempty(varargin)),
    InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
else
    InPar = DefV;
end


if (isempty(InPar.InType)),
    % attempt to identify input type automatically
    if (isnumeric(File)),
        InType = 'numeric';
    elseif (iscell(File)),
        % File is cell
        if (iscellstr(File)),
            % cell array of file names
            % try to guess file format:
            InType = identify_file_type(File{1});
            
        else
            % assume - cell array of numeric
            InType = 'numeric-cell';
            
        end
        
    elseif (ischar(File)),
        % identify file type based on extension
        InType = identify_file_type(File);
        
    elseif (isastspec(File)),
        % input is AstSpec
        InType = 'astspec';
    else
        error('Illegal input File type');
    end
    
else
    InType = InPar.InType;
end


% read files
switch lower(InType)
    case 'astspec'
        % Input is AstSpec object - do nothing
        AS = File;
    case 'ascii'
        % Input is an ascii file
        AS = spec_read_ascii(File,varargin{:});
    case 'mat'
        % input is a mat file
        AS = spec_read_mat(File,varargin{:});
    case {'fits','fits-wcs','fits-mat','fits-bintab','fits-asciitab'}
        % several fits types...
        switch lower(InType)
            case 'fits'
                % attempt to identify FITS type
                
            otherwise
                % do nothing
        end
        
        
        
        
    case {'hdf5','hd5','h5'}
        
    case 'numeric'
        % input is numeric
        AS = AstSpec;  % define an AstSpec class object
        AS = mat2spec(AS,File);
        
    case 'numeric-cell'
        % input is cell of numerics
        Nc = numel(File);
        for Ic=1:1:Nc,
            AS(Ic) = mat2spec(AS,File{Ic});
        end
        
    otherwise
        error('unknown InType option');
end
        

function InType=identify_file_type(File)
    % Description: identify file name by its extension.
    Splitted = regexp(File,'\.','split');
    switch lower(Splitted{end}),
        case {'fits','fit'}
            warning('Assume file type is FITS');
            InType = 'fits';
        case {'hdf5','hd5','h5'}
            warning('Assume file type is HDF5');
            InType = 'hdf5';
        case {'txt','ascii','flm'}
            warning('Assume file type is ASCII');
            InType = 'ascii';
        case 'mat'
            warning('Assume file type is MAT');
            InType = 'mat';
        otherwise
            error('Can not identify file type automatically');
    end
end
   
% end of spec_read function
end           
            
         