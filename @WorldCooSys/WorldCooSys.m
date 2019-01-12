%--------------------------------------------------------------------------
% WorldCooSys class                                                  class
% Description: A class of for World Coordinate System (WCS).
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef WorldCooSys
    properties (SetAccess = public)
        WCS
        UserData
    end
    
   
    methods

        %-------------------------
        %--- Class constructor ---
        %-------------------------
        function WCS=WorldCooSys(N,M)
            % Description: WorldCooSys constructor method
            
            Field = 'WCS';
            
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
                    WCS(I,J).(Field) = [];
                end
            end
        end
        
        
        
       
        
        
        %--------------------------
        %--- Structre functions ---
        %--------------------------
        function obj=isfield(WorldC,Field)
            % isfield 
            obj = any(strcmp(fieldnames(WorldC),Field));
        end

        function obj=isstruct(Head)
            obj = true;  %isstruct(Sim) || isa(Sim,'SIM');
        end

    end
    
    
    %----------------------
    %--- Static methods ---
    %----------------------
    methods (Static)
        function Out=read_ctype(String,Fields)
            % Parse WCS CTYPE header keyword
            % Description: Given a FITS header CTYPE keyword value (e.g., 'RA---TAN')
            %              return the coordinate type (e.g., RA) and transformation
            %              type (e.g., 'TAN').
            % Input  : - A string containing a CTYPE header keyword value.
            %            Alternatively, this can be a structure array which contains
            %            a CTYPE1 or a CTYPE2 fields (or other specified fields).
            %          - In case the first input argument is a structure then this
            %            is a cell array of field names to extract.
            %            Default is {'CTYPE1','CTYPE2'}.
            % Output : - If the first input is a tring then this is a structure
            %            containing the fields .Coo and .Tran. and optionally .Dist
            %            If the first input is a structure array then this is
            %            a structure array containing field names (e.g., 'CTYPE1')
            %            which bythemselfs containing the fields .Coo and .Tran.
            %            and optionally .Dist
            % Tested : Matlab R2011b
            %     By : Eran O. Ofek                    Nov 2013
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: Out=WorldCooSys.read_ctype('RA---TAN');
            % Reliable: 2
            %--------------------------------------------------------------------------

            Def.Fields = {'CTYPE1','CTYPE2'};
            if (nargin==1)
                Fields  = Def.Fields;
            end
            if (~iscell(Fields))
                Fields = {Fields};
            end

            if (ischar(String))
                Split    = regexp(String,'-','split');
                Pair     = Split(~Util.cell.isempty_cell(Split));
                Out.Coo  = Pair{1};
                Out.Tran = Pair{2};
                if (numel(Pair)>2)
                    Out.Dist = Pair{3};
                end
            elseif (isstruct(String))
                N = numel(String);
                for I=1:1:N
                    for If=1:1:numel(Fields)
                        Split    = regexp(String.(Fields{If}),'-','split');
                        Pair     = Split(~Util.cell.isempty_cell(Split));
                        Out(I).(Fields{If}).Coo  = Pair{1};
                        Out(I).(Fields{If}).Tran = Pair{2};
                        if (numel(Pair)>2)
                            Out(I).(Fields{If}).Dist = Pair{3};
                        end
                    end
                end
            else
                error('Unknown String type');
            end

        end
        
        function Ans=isWorldCooSys(Obj)
            % Return true if object is WorldCooSys
            % Description: Check if object is of WorldCooSys class.
            % Input  : - Object
            % Output : - {true|false}.
            %     By : Eran O. Ofek                    Oct 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: S=1; WorldCooSys.isWorldCooSys(S);
            % Reliable: 2
            Ans = isa(Obj,'WorldCooSys');

        end
        
        function WCS=struct2wcs(St)
            % Convert a structure array to a WorldCooSys object
            % Input  : - A structure array.
            % Output : - A WorldCooSys object.
            % License: GNU general public license version 3
            %     By : Eran O. Ofek                    Mar 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: WCS=WorldCooSys.struct2wcs(St);
            % Reliable: 2

            WCS = WorldCooSys(size(St));
            Nw  = numel(WCS);
            for Iw=1:1:Nw
                WCS(Iw).WCS = St(Iw);
            end
        end
        
      
    end
    
    % transformations
    methods (Static)
        
        function CD=read_cd_mtarix(WCS)
            
           %if ('WCSAXES'
            
        end
        
        % not working
        function X=pix2intermediate(P,Ref,RotMat)
            %
            % Input  : - A matrix with pixel coordinates.
            %            Each column represent one axis.
            %          - Reference coordinates of the transformation.
            %            I.e., CRPIX
            %          - Rotation matrix. Either:
            %            CD matrix or a structure containing either,
            %            CD, {PC and SCALE}, {CDELT, CROT} parameters.
            %            
            
            Naxes = size(P,2);
            
            FN = fieldnames(WCS);
            if (~isempty(strfind(FN,'CD1_1')))
                Str = 'CD';
            elseif (~isempty(strfind(FN,'PC1_1')))
                Str = 'PC';
            else
                error('Unsupported WCS rotation matrix');
            end
            for Iaxes1=1:1:Naxes
                for Iaxes2=1:1:Naxes
                    KeyName = sprintf('%s%d_%d',Str,Iaxes1,Iaxes2);
                    CD(Iaxes1,Iaxes2) = WCS.(KeyName);
                    switch lower(Str)
                        case 'PC'
                            KeyScale = sprintf('CDELT%d',Iaxes1);
                            CD(Iaxes1,Iaxes2) = CD(Iaxes1,Iaxes2).*WCS.(KeyScale);
                    end
                end
            end
            
            
            
            if (~isempty(strfind(FN,'CD1_1')))
                % CD matrix is available
                CD = [WCS.CD1_1, WCS.CD1_2; WCS.CD2_1, WCS.CD2_2];
            else
                if (~isempty(strfind(FN,'PC1_1')))
                    CD = [WCS.PC1_1, WCS.PC1_2; WCS.PC2_1, WCS.PC2_2];
                    CD = CD.*[WCS.CDELT1; WCS.CDELT2];
                else
                    
                end
            end
            
            
            for Iaxes=1:1:Naxes
                
             
                
            end
            
            X = RotMat*(P - Ref(:).');
            
            
        end
            
        % not working  
        function pix2native(P,Proj,PV) 
            %
            
            switch lower(Proj)
                case 'tan'
                    
                otherwise
                    error('unknown transformation type');
            end
                    
        end
        
        % not working
        function distortion_sip_poly2keys(PolyOrder,PolyCoef)
            
        end
        
        
        
    end
    
    
end

            
