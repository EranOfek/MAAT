function AstC=get_astcat(CatName,RA,Dec,Radius,varargin)
%--------------------------------------------------------------------------
% get_astcat function                                            Catalogue
% Description: Get astronomical catalog from local disk or WWW into an
%              AstCat object.
% Input  : - Catalog file name, or catalog function handle.
%            This can be either the mat file name, or HDF5 file name
%            containing the catalog.
%            HDF5 catalogs must have an .h5, .hd5, or .hdf5 extension.
%            Alternatively, this can be a function handle:
%            F(RA,Dec,Radius,Shape)
%          - J2000.0 R.A. (scalar) in radians, or [H M S], or sexagesimal
%            string.
%          - J2000.0 Dec. (scalar) in radians, or [sign H M S], or
%            sexagesimal string.
%          - Radius. Radius units are specified in the 'RadUnits' argument.
%            Default units is radians.
%            Default search radius is 1./(60.*RAD) (i.e., 1 arcmin).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Shape'   - 'circ'|'box'
%            'OutType' - 'astcat'|'mat'|'table'|'struct'
%            'RadUnits'- 'mas'|'arcsec'|'arcmin'|'deg'|'rad'
%            'UseLoadCheck'
% Output : - The catalog.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

CatField  = 'Cat';

Def.Radius = 1./(RAD.*60);
if (nargin==3),
    Radius  = Def.Radius;
end

DefV.Shape              = 'circ';
DefV.OutType            = 'astcat';
DefV.RadUnits           = 'rad';
DefV.UseLoadCheck       = true;
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

% Convert Radius to radians
Radius = convert.angular(InPar.RadUnits,'rad',Radius);  % [rad]



if (ischar(CatName)),
    % CatName is a string - try to load the file
    if (strcmp(CatName(end-2:end),'h5') || strcmp(CatName(end-3:end),'hd5') || strcmp(CatName(end-4:end,'hdf5')),
        % Assume catalog is an HDF5 file
        error('hdf5 catalogs are not supported yet - need to define');
    else
        % Assume catalog is a MAT file
        if (InPar.UseLoadCheck),
            Cat = load_check(CatName);
        else
            Cat = load2(CatName);
        end
        % Cat is not necesserly an AstCat object
        % Convert various formats into an AstCat object
        if (isastcat(Cat)),
            % Cat is already an AstCat object
            AstC = Cat;
        else
            if (isstruct(Cat)),
                % Cat is a structure 
                % assumes that Cat is an AstCat like object
                warning('OutType was changed into a struct');
                AstC = Cat;
            elseif (isnumeric(Cat)),
                % construct an AstCat object
                AstC = AstCat;
                AstC.(CatField)   = Cat;
                % no column information
            else
                error('Unsupported catalog format');
            end
        end
    end
    % search the catalog
    % check if catalog is an HTM index catalog
    
    
    
    
    
elseif (isa(CatName,'function_handle')),
    % CatName is a function handle
    
    AstC = CatName(RA,Dec,Radius,InPar.Shape,...
        
else
    error('Unknown CatName type');
end
                
            
        
    
