function Cat=build_catalog_kdtree(CatName)
%--------------------------------------------------------------------------
% build_catalog_kdtree function                                  Catalogue
% Description: 
% Input  : - 
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------



CatPathName = which(CatName);
Pos         = strfind(CatPathName,filesep);
CatPath     = CatPathName(1:Pos(end));

PWD = pwd;
cd(CatPath);

if (strcmp(CatName(end-3:end),'.mat')),
    % this is a mat file
    % assume this is a matlab astronomical catalog
    CatBaseName = CatName(1:end-4);
    Cat = load2(CatName);

    % make sure catalog is sorted by Declination
    Cat.Cat = sortrows(Cat.Cat,Cat.Col.Dec);
    Cat.SortedBy    = 'Dec';
    Cat.SortedByCol = Cat.Col.Dec;
    
    % check fields
    Cat=check_fields(Cat);
    
    Cat=add_kdtree(Cat);
    
    %Cat=trans2table(Cat);
    
    % save catalog
    eval(sprintf('%s=Cat',CatBaseName));
    eval(sprintf('save %s %s',CatName,CatBaseName));

elseif (strcmp(CatName(end-4:end),'.fits')),
    error('not working')
    % assume FITS binary table 
    % build catalog from scratch
   
    CatBaseName = CatName(1:end-5);

    [~,~,ColCell,Col,Table]=get_fitstable_col(CatName);
    Ncol = length(Table);
    T = table(Table{:});
    Cat.Cat = table2array(T);
    
    
    % check for RA/Dec columns:
    CheckFields = {'RA','Dec'};
    Ncf = length(CheckFields);
    for Icf=1:1:Ncf,
        fprintf('Check field: %s\n',CheckFields{Icf});
        if (isempty(find(strcmp(ColCell,CheckFields{Icf})))),
            fprintf('Cannot find %s field\n',CheckFields{Icf});
            fprintf('List of fields in FITS binary table\n');
            ColCell
            String = input('Input name of field name as appear in the list : ','s');
            
            Icol = strcmp(ColCell,String);
            ColCell{Icol} = CheckFields{Icf};
        end
    end
    Col = cell2struct(num2cell(1:1:Ncol),ColCell,2);
    Cat.Col = Col;
    Cat.ColCell = ColCell;
    
    % make sure catalog is sorted by Declination
    Cat.Cat = sortrows(Cat.Cat,Cat.Col.Dec);
    Cat.SortedBy    = 'Dec';
    Cat.SortedByCol = Cat.Col.Dec;
    
    Cat.ColUnits    = {};
    Cat.ColComments = {};
    Cat.Source      = {};
    Cat.Reference   = {};
    Cat.Version     = {};
    
    Cat=add_kdtree(Cat);
    
end

% return to original directory
cd(PWD);




function Cat=check_fields(Cat)
    % check fields
    FN = fieldnames(Cat);
    Fields = {'Cat','Col','ColCell','ColUnits','SortedBy','SortedByCol','ColComments','Source','Reference','Version','KDTree'};
    for If=1:1:length(Fields),
        if (isempty(find(strcmp(FN,Fields{If}), 1))),
            fprintf('Field: %s - doesnt exist in Cat\n',Fields{If});
            switch Fields{If}
                case {'SortedBy','Source','Reference'}
                    String = input('Enter the missing field: ','s');
                    Cat.(Fields{If}) = String;
                otherwise
                    % do nothing
            end
            
        end
    end


        
function Cat=add_kdtree(Cat)
    % cosine directions
    [CosX,CosY,CosZ] = coo2cosined(Cat.Cat(:,Cat.Col.RA),Cat.Cat(:,Cat.Col.Dec));
    KDTree = kdtree([CosX,CosY,CosZ]);
    Cat.KDTree = KDTree;
    
    Ncol = size(Cat.Cat,2);
    Cat.Cat(:,Ncol+[1 2 3]) = [CosX, CosY, CosZ];
    Cat.Col.CosX = Ncol+1;
    Cat.Col.CosY = Ncol+2;
    Cat.Col.CosZ = Ncol+3;
    Cat.ColCell{Ncol+1} = 'CosX';
    Cat.ColCell{Ncol+2} = 'CosY';
    Cat.ColCell{Ncol+3} = 'CosZ';
    Cat.ColUnits{Ncol+1} = '';
    Cat.ColUnits{Ncol+2} = '';
    Cat.ColUnits{Ncol+3} = '';




function Cat=trans2table(Cat)
    % transform to table
    if (~istable(Cat.Cat)),
        Cat.Cat = array2table(Cat.Cat);
        Cat.Cat.Properties.VariableNames = Cat.ColCell;
        Cat.Cat.Properties.VariableUnits = Cat.ColUnits;
        String = input('Enter catalog description : ','s');
        Cat.Cat.Properties.Description   = String;
        
    end
