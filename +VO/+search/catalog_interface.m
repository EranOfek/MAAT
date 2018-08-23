function varargout=catalog_interface(Catalog,Path,varargin)
% An interface auxilary to the catalogs in the data directory
% Package: VO.search
% Description: An auxilary function to interface the local catalogs in the
%              data directory as well as external catalogs.
%              This function is responsible for operations like:
%              datacats.binaries.SB9
% Input  : - Catalog file name, or function handle.
%          - Catalog path.
%          * Either ('q',QueryString) or (RA,Dec,Radius).
% Output : - AstCat object
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% Cat=VO.search.catalog_interface('SNR.mat','/home/eran/matlab/data/+datacats/+SN')
% Reliable: 2
%--------------------------------------------------------------------------

MatFileLoader = @Util.IO.load_check;  % other option is: Util.IO.load2

Nvar = numel(varargin);

if (strfind(Catalog,'.mat'))
    %--------------------------
    %--- File is a mat file ---
    %--------------------------
    
    % load catalog (unless requested foe help)
    OnlyHelp = false;
    if (~isempty(varargin))
        if (strcmp(varargin{1}(1),'h'))
            OnlyHelp = true;
        end
    end
    if ~OnlyHelp
        varargout{1} = MatFileLoader(sprintf('%s%s%s',Path,filesep,Catalog));
    end
    
    if (Nvar==0)
        % varargin is empty
        % only load catalog
    else
        
        
        if (Nvar>0)
            IsCoo   = true;
            IsQuery = false;
            if (ischar(varargin{1}))
                IsCoo = false;
                switch lower(varargin{1}(1))
                    case 'h'
                        % show help
                        % call the show help (private) function
                        catalog_interface_show_help;

                    case 'q'
                        % query catalog
                        IsQuery = true;
                    case 'c'
                        % show catalog columns
                        catalog_interface_show_columns(varargout{1});

                    case 'i'
                        % show catalog information
                        catalog_interface_show_info(varargout{1});

                    otherwise
                        % maybe first argument is sexagesimal
                        % coordinates
                        IsCoo = true;
                end
            end
        end

            
        if (IsCoo && Nvar>1)
            % coordinates cone search

            % Default values
            RAD = 180./pi;
            Def.Radius  = 60./(3600.*RAD);   % 60 arcsec

            RA  = varargin{1};
            Dec = varargin{2};
            if (Nvar>2)
                Radius = varargin{3};
            else
                Radius = Def.Radius;
            end
            Radius = convert.angular('arcsec','rad',Radius);
            % convert RA/Dec to radians
            if (ischar(RA) || iscellstr(RA))
                RA = celestial.coo.convertdms(RA,'SH','r');
            else
                RA = convert.angular('deg','rad',RA);
            end
            if (ischar(Dec) || iscellstr(Dec))
                Dec = celestial.coo.convertdms(Dec,'SD','R');
            else
                Dec = convert.angular('deg','rad',Dec);
            end

            if (AstCat.isastcat(varargout{1}) )
                % Catalog is stored as an AstCat object sorted by
                % declination
                
                ColNameRA  = 'RA';
                ColNameDec = 'Dec';
                CatField   = AstCat.CatField;
                ColRA  = colname2ind(varargout{1},ColNameRA);
                ColDec = colname2ind(varargout{1},ColNameDec);

                Ncoo=numel(RA);
                AllCat = varargout{1};
                if (istable(AllCat.(CatField)))
                    Mat = [AllCat.(CatField).(ColNameRA), AllCat.(CatField).(ColNameDec)];
                    for Icoo=1:1:Ncoo
                        Ind = VO.search.search_sortedlat(Mat,RA(Icoo),Dec(Icoo),Radius);
                        varargout{1}(Icoo) = varargout{1}(1);  % copy info
                        varargout{1}(Icoo).(CatField) = AllCat.(CatField)(Ind,:);
                    end
                else
                    for Icoo=1:1:Ncoo
                        Ind = VO.search.search_sortedlat(AllCat.(CatField)(:,[ColRA, ColDec]),RA(Icoo),Dec(Icoo),Radius);
                        varargout{1}(Icoo) = varargout{1}(1);  % copy info
                        varargout{1}(Icoo).(CatField) = AllCat.(CatField)(Ind,:);
                    end
                end
                
            elseif (isnumeric(varargout{1}))
                % numeric arrays
                % only possible to load
                fprintf('   Catalog is stored as an array - only load operation is possible\n');

            end
        end
        
        if (IsQuery && Nvar>1)
            % Execute a query string
            varargout{1}=query(varargout{1},varargin{2});
            
        end
    end
    
end

        
end

%--- Functions ---

function catalog_interface_show_help
    % catalog interface show help
    
    fprintf('\n')
    fprintf('--- Helpe for catalog interface ---\n');
    fprintf('h - show this help\n');
    fprintf('c - show catalog columns\n');
    fprintf('q - execute a query on catalog columns (e.g., RA>pi & Dec<0)\n');
    fprintf('\n');

end

function catalog_interface_show_columns(AstC)
    % catalog interface show catalog columns
    
    ColCellField  = AstCat.ColCellField;
    ColUnitsField = AstCat.ColUnitsField; 
    fprintf('\n');
    Ncol = numel(AstC.(ColCellField));
    for Icol=1:1:Ncol
        if isempty(AstC.(ColUnitsField))
            Units = '';
        else
            Units = AstC.(ColUnitsField){Icol};
        end
        fprintf('%20s - %-20s\n',AstC.(ColCellField){Icol},Units);
    end
    fprintf('\n');
    
end

function catalog_interface_show_info(AstC)
    % catalaog interface show information
    
    fprintf('\n');
    fprintf('Catalog name    : %s\n',AstC.Name);
    fprintf('Catalog source  : %s\n',AstC.Source);
    fprintf('Catalog version : %s\n',AstC.Version);
    fprintf('\n');
    
end
