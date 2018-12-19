%--------------------------------------------------------------------------
% AstCat class                                                       class
% Description: A class of structure array of catalogs (AstCat).
%              Note that this class is a subset of the SIM class.
% Input  : null
% Output : null
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------


classdef AstCat < HEAD
    properties (SetAccess = public)
        Cat
        Col
        ColCell
        ColUnits
        SortedBy
    end
    properties (Hidden = true)
        SortedByCol
        Name
        Source
        Reference
        Version
        %UserData
    end
  
    
    %-------------------
    %--- Constractor ---
    %-------------------
    methods
        
        function AstC=AstCat(N,M)
            % AstCat class constructor
            % Package: @AstCat
            % Description: AstCat class constructor.
            % Input  : - Number of rows, or [row, columns] Default is 1.
            %          - Number of columns. Default is 1.
            % Output : - An AstCat object of the requested size.


            CatField = 'Cat';
            
            if (nargin==0)
                N = 1;
                M = 1;
            elseif (nargin==1)
                if (numel(N)>1)
                    M = N(2);
                else
                    M = 1;
                end
            else
                % do nothing
            end

            for I=1:1:N
                for J=1:1:M
                    AstC(I,J).(CatField) = [];
                end
            end
            
        end
        

        
        function AstC=clear_cat(AstC)
            % Clear the catalog fron an AstCat object
            % Package: @AstCat
            % Description: Clear the catalog fron an AstCat object
            % Input  : - An AstCat object
            % Output : - An AstCat object

            CatField     = AstCat.CatField;
            ColField     = AstCat.ColField;
            ColCellField = AstCat.ColCellField;

            Ncat = numel(AstC);
            for Icat=1:1:Ncat
                AstC(Icat).(CatField)     = [];
                AstC(Icat).(ColField)     = [];
                AstC(Icat).(ColCellField) = [];
            end
        end

                
    end
    
    
    %----------------------
    %--- Static methods ---
    %----------------------
    % get field names
    methods (Static)
        
        function AstC=struct2astcat(St)
            % Convert a structure array to AstCat array
            % Package: @AstCat
            % Input  : - A structure array
            % Output : - An AstCat array
            
            Nst=numel(St);
            AstC=AstCat(size(St));
            Fields = fieldnames(AstC);
            Nf     = numel(Fields);
            for Ist=1:1:Nst
                for If=1:1:Nf
                    if (isfield(St,Fields{If}))
                        AstC(Ist).(Fields{If}) = St(Ist).(Fields{If});
                    end
                end
            end
                    
            
        end
        
        function Name=CatField
            % Return the Cat field name in AstCat
            % Package: @AstCat
            % Description: Return the Cat field name in AstCat
            % Input  : null
            % Output : The cat field name
            Name = 'Cat';
        end
        function Name=ColField
            % Return the Col field name in AstCat
            % Package: @AstCat
            % Description: Return the Col field name in AstCat
            % Input  : null
            % Output : The Col field name
            Name = 'Col';
        end
        function Name=ColCellField
            % Return the ColCell field name in AstCat
            % Package: @AstCat
            % Description: Return the ColCell field name in AstCat
            % Input  : null
            % Output : The ColCell field name
            Name = 'ColCell';
        end
        function Name=ColUnitsField
            % Return the ColUnits field name in AstCat
            % Package: @AstCat
            % Description: Return the ColUnits field name in AstCat
            % Input  : null
            % Output : The ColUnits field name
            Name = 'ColUnits';
        end
    end
    
    
    methods (Static)
        function AstC=sim2astcat(Sim)
            % Convert A SIM or AstCat object to an AstCat object
            % Package: @AstCat
            % Description: Copy the elements of a SIM object into an AstCat object.
            % Input  : - A SIM or AstCat class object.
            % Output : - A AstCat class object.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Mar 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: AstCat.sim2astcat(SIM);
            % Reliable: 2


            AstC       = AstCat(size(Sim));
            FieldNames = fieldnames(AstC(1));
            Nf         = numel(FieldNames);
            Nsim       = numel(Sim);
            for Ic=1:1:Nsim
                for If=1:1:Nf
                    AstC(Ic).(FieldNames{If}) = Sim(Ic).(FieldNames{If});
                end
            end
        end
        
    end % static
    
    % load hdf5 into AstCat
    methods (Static)
        
        % load HDF5 (generated by AstCat.saveh) to AstCat object
        function AstC=loadh2astcat(FileName,CatIndex)
            % Description: Load an HDF5 file (created by AstCat.saveh) into
            %              an AstCat object.
            % Input  : - HDF5 file name.
            %          - Catalog index or name within HDF5 file. If empty then load
            %            all catalogs. Default is empty.
            % Output : - An AstCat object with the loaded catalogs.
            % Example: AstC=AstCat.loadh2astcat('SDSS_DR9_Fields_All_PolySort.hdf5');
            % Reliable: 2
            
            CatField = AstCat.CatField;
            ColField = AstCat.ColField;
            
            if (nargin==1)
                CatIndex = [];
            end
            if (ischar(CatIndex))
                CatName = CatIndex;
            else
                if (~isempty(CatIndex))
                    CatName = sprintf('Cat%d',CatIndex);
                else
                    CatName = CatIndex;
                end
            end
            
            [Data,Att]=Util.IO.loadh(FileName,CatName,'s');
            
            FN  = fieldnames(Data);
            Nfn = numel(FN);
            AstC = AstCat(Nfn,1);
            for Ifn=1:1:Nfn
                AstC(Ifn).(CatField) = Data.(FN{Ifn});
                AstC(Ifn).(ColField) = Att(Ifn);
                AstC(Ifn)  = col2colcell(AstC(Ifn));
            end
                
            
        end
           
    end % static
    
    % arrays
    methods (Static)
        function [AstC]=array2astcat(Array,ColCell,varargin)
            % Convert an array to AstCat object   
            % Description: Convert an array or table into an AstCat object.
            % Input  : - Matrix or table, or a cell array of matrix or tables.
            %            If cell array, then each element will populate an Astcat
            %            element.
            %          - Cell array of column names.
            %            If empty will attempt to read from table. Default is empty.
            %          * Arbitrary number of pairs of arguments: ...,field,value,...
            %            where field is any AstCat field. E.g.,
            %            ...'ColUnits',{'rad','rad'},...
            % Output : - An AstCat object.
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Mar 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [AstC]=AstCat.array2astcat(rand(100,3),{'A','B','C'},'Source','This source')
            % Reliable: 2
            %--------------------------------------------------------------------------

            if (nargin<2)
                ColCell = {};
            end

            if (~iscell(Array))
                Array = {Array};
            end
            Ncell = numel(Array);

            AstC = AstCat(size(Array));
            Nvar = numel(varargin);

            for Icell=1:1:Ncell

                AstC(Icell).Cat = Array{Icell};
                if (isempty(ColCell))
                    if (istable(Array{Icell}))
                        AstC(Icell).ColCell = Array{Icell}.Properties.VariableNames;
                    end
                else
                    AstC(Icell).ColCell = ColCell;
                end

                AstC(Icell) = colcell2col(AstC(Icell));

                for I=1:2:Nvar-1
                    AstC(Icell).(varargin{I}) = varargin{I+1};
                end
            end
        end
    end
    
    % static methods
    methods (Static)
        % isastcat
        function Ans=isastcat(Obj)
            % Return true if AstCat object
            % Description: Check if an object is an AstCat class object.
            % Input  : - An object.
            % Output : - True if object is AstCat, otherwise false.
            %     By : Eran O. Ofek                    Mar 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: Ans=isastcat(AstC);
            % Reliable: 2

            Ans = isa(Obj,'AstCat');
        end
        
        
        
    end % static methods
    
    % structure related functions
    methods

       

        %--- Structre functions ---
        function obj=isfield(Sim,Field)
            % isfield 
            obj = any(strcmp(fieldnames(Sim),Field));
        end

        function obj=isstruct(Sim)
            obj = true;  %isstruct(Sim) || isa(Sim,'SIM');
        end

    
     
    end
end

            
