%--------------------------------------------------------------------------
% ClassWCS class                                                     class
% Description: A class of for World Coordinate System (WCS).
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef ClassWCS
    properties (SetAccess = public)
        WCS
        UserData
    end
    
    % class creation
    methods

        %-------------------------
        %--- Class constructor ---
        %-------------------------
        function W=ClassWCS(N,M)
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
                    W(I,J).(Field) = [];
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
    
    % class
    methods (Static)
        function Res=isClassWCS(W)
            % Return true if a ClassWCS object
            % Package: @ClassWCS
            % Input  : - A variable.
            % Output : - A true fror ClassWCS object, false otherwise
            % Example: ClassWCS.isClassWCS(1)
            
            if (isa(W,'ClassWCS'))
                Res = true;
            else
                Res = false;
            end
        end
    end
    
    % fields
    methods (Static)
        function F=WCSField
            % Return the WCS field name (i.e., 'WCS')
            % Package: @ClassWCS
            F = 'WCS';
        end
        
        
    end
    
    % populate and get data
    methods (Static) 
        % populate ClassWCS
        function W=populate(Header)
            % Populate the ClassWCS object with key/val from header
            % Package: @ClassWCS
            % Description: Populate the ClassWCS object with key/val from header
            % Input  : - HEAD class object or a 3 column cell array of
            %            header keywords.
            % Output : - A populated ClassWCS object.
            
            Default.CUNIT = 'deg';
            Default.CTYPE1 = 'RA---TAN';
            Default.CTYPE2 = 'DEC--TAN';
            
            HeaderField = HEAD.HeaderField;
            WCSField    = 'WCS';
            
            if (iscell(Header))
                H = HEAD;
                H.(HeaderField) = Header;
            else
                H = Header;
            end
            Nh = numel(H);
            W  = ClassWCS(Nh,1);
            
            KeysSingle = {'RADECSYS','LONPOLE','LATPOLE','EQUINOX'};
            Nsin       = numel(KeysSingle);
            KeysN      = {'CTYPE','CUNIT','CRPIX','CRVAL','CDELT'};
            Nn         = numel(KeysN);
            
            
            for Ih=1:1:Nh
                % for each header
                % Read number of axes
                % if WCSAXES is not available use NAXES as default
                Naxes = getkey(H(Ih),'WCSAXES');
                if (isempty(Naxes))
                    Naxes = getkey(H(Ih),'NAXIS');
                else
                    if (isnan(Naxes{1}))
                        Naxes = getkey(H(Ih),'NAXIS');
                    end
                end
                Naxes = Naxes{1};
                
                % read keywords from the KeysSingle list
                ValSingle = mgetkey(H(Ih),KeysSingle);
              
                TmpCtype = getkey(H(Ih),'CTYPE1');
                
                % concat
                KeyNames = {'WCSAXES', KeysSingle{:}};
                KeyVal   = {Naxes, ValSingle{:}};
                
                W(Ih).(WCSField) = cell2struct(KeyVal,KeyNames,2);
                
              
%<<<<<<< HEAD
                if isnan(Naxes) || isnan(TmpCtype{1}(1))
                    % deal with missing WCS keywords
                    W(Ih).(WCSField).Status = false;
                    W(Ih).(WCSField).CD = nan(2,2);
                    W(Ih).(WCSField).CRPIX = nan(1,2);
                    W(Ih).(WCSField).CRVAL = nan(1,2);
                    W(Ih).(WCSField).CDELT = nan(1,2);
                    W(Ih).(WCSField).CTYPE = {'RA---TPV','DEC--TPV'};
                    W(Ih).(WCSField).CUNIT = {'deg','deg'};
                else
                    W(Ih).(WCSField).Status = true;
                    
                    % read Keywords from the KeysN list
                    KeyNname = cell(1,Nn.*Naxes);
                    K = 0;
                    for In=1:1:Nn
                        for Iaxis=1:1:Naxes
                            K = K + 1;
                            KeyNname{K} = sprintf('%s%d',KeysN{In},Iaxis);
                        end
%=======
                    end
                %end
                
                    % read Keywords from the KeysN list
                    KeyNname = cell(1,Nn.*Naxes);
                    K = 0;
                    for In=1:1:Nn
                        for Iaxis=1:1:Naxes
                            K = K + 1;
                            KeyNname{K} = sprintf('%s%d',KeysN{In},Iaxis);
                        end
                    end
                    ValN = mgetkey(H(Ih),KeyNname);
                    ValN = reshape(ValN,2,Nn);
                    for In=1:1:Nn
                        % fixing a bug found by Na'ama
                        switch lower(KeysN{In})
                            case 'cunit'
                                if (any(isnan(ValN{1,In})))
                                    % CUNIT is not populated in header
                                    % set to default
                                    ValN{1,In} = Default.CUNIT;
                                    ValN{2,In} = Default.CUNIT;
                                end

                            case 'ctype'
                                if (any(isnan(ValN{1,In})))
                                    % CUNIT is not populated in header
                                    % set to default
                                    ValN{1,In} = Default.CTYPE1;
                                    ValN{2,In} = Default.CTYPE2;
                                end
                        end

                        if (iscellstr(ValN(:,In)))
                            W(Ih).(WCSField).(KeysN{In}) = ValN(:,In).';
                        else
                            W(Ih).(WCSField).(KeysN{In}) = cell2mat(ValN(:,In)).';
    %>>>>>>> d3d1fd3e53a5851582211798c8cdcd679ba36ecd
                        end
                        ValN = mgetkey(H(Ih),KeyNname);
                        ValN = reshape(ValN,2,Nn);
                        for In=1:1:Nn
                            % fixing a bug found by Na'ama
                            switch lower(KeysN{In})
                                case 'cunit'
                                    if (any(isnan(ValN{1,In})))
                                        % CUNIT is not popuklated in header
                                        % set to degault
                                        ValN{1,In} = Default.CUNIT;
                                        ValN{2,In} = Default.CUNIT;
                                    end
                            end

                            if (iscellstr(ValN(:,In)))
                                W(Ih).(WCSField).(KeysN{In}) = ValN(:,In).';
                            else
                                W(Ih).(WCSField).(KeysN{In}) = cell2mat(ValN(:,In)).';
                            end

                        end


                        % Read The CD/PC matrix
                        KeysCD = cell(1,Naxes.^2);
                        KeysPC = cell(1,Naxes.^2);
                        K = 0;
                        for Iaxes1=1:1:Naxes
                            for Iaxes2=1:1:Naxes
                                K = K + 1;
                                KeysCD{K} = sprintf('CD%d_%d',Iaxes1,Iaxes2);
                                KeysPC{K} = sprintf('PC%d_%d',Iaxes1,Iaxes2);
                            end
                        end

                        ValCD = mgetkey(H(Ih),KeysCD);
                        K = 0;
                        CD = nan(Naxes,Naxes);
                        for Iaxes1=1:1:Naxes
                            for Iaxes2=1:1:Naxes
                                K = K + 1;
                                CD(Iaxes1,Iaxes2) = ValCD{K};
                            end
                        end

                        % bug fix - treat cases in whic not all CD keywords are
                        % provided - assume no rotation.
                        if (any(isnan(CD(:))) && ~all(isnan(CD(:))))
                            CD(isnan(CD)) = 0;
                        end



                        if (any(isnan(CD(:))) || isempty(CD))
                            % CD is empty try to read PC
                            ValCD = mgetkey(H(Ih),KeysPC);
                            K = 0;
                            ScaleName = sprintf('CDELT');
                            for Iaxes1=1:1:Naxes
                                %ScaleName = sprintf('CDELT%d',Iaxes1);
                                for Iaxes2=1:1:Naxes
                                    K = K + 1;
                                    CD(Iaxes1,1) = ValCD{K}.*W(Ih).(WCSField).(ScaleName)(Iaxes1);
                                end
                            end

                        end
                        W(Ih).(WCSField).CD = CD;


                        % Read distortions

                        % look for PV coeficients
                        FlagMatchPV = ~Util.cell.isempty_cell(regexp(H(Ih).(HeaderField)(:,1),'PV\d+\_\d+','match'));


                        Names  =regexp(H(Ih).(HeaderField)(FlagMatchPV,1), 'PV(?<D1>\d+)\_(?<D2>\d+)','names');
                        Nnames = numel(Names);
                        PV_Ind = zeros(Nnames,2);
                        for Inames=1:1:Nnames
                            PV_Ind(Inames,:) = [str2double(Names{Inames}.D1), str2double(Names{Inames}.D2)];
                        end

                        W(Ih).(WCSField).PV.Ind     = PV_Ind;
                        W(Ih).(WCSField).PV.KeyVal  = H(Ih).(HeaderField)(FlagMatchPV,2);
                        W(Ih).(WCSField).PV.KeyName = H(Ih).(HeaderField)(FlagMatchPV,1);

                        % look for SIP coeficients
                        % TBD
                    end
                    
                end

            end
                
        end
    
        
        
        function H=sip2head(StructSIP,OutType)
            % Write a SIP structure into an HEAD object
            % Package: @ClassWCS
            % Description: Write SIP structure keyword values into an HEAD
            %              object.
            % Input  : - A SIP structure.
            %          - Output type {'cell','HEAD'}.
            %            Default is 'HEAD'.
            % Output : - An HEAD object with the SIP parameters
            % Example: H = ClassWCS.sip2head(StructSIP);
            
            if (nargin<2)
                OutType = 'HEAD';
            end
            
            Fields = {'A','B','AP','BP'};
            Nf     = numel(Fields);
            
            HeadCell = zeros(0,3);
            K = 0;
            if (StructSIP.Exist)
                for If=1:1:Nf
                
                    if (isfield(StructSIP.matrix,Fields{If}))
                        Mat.(Fields{If}) = StructSIP.matrix.(Fields{If});
                    else
                        Mat.(Fields{If}) = [];
                    end
                    
                    if (~isempty(Mat.(Fields{If})))
                        FieldOrderStr = sprintf('%s_ORDER',Fields{If});
                        K = K + 1;
                        [HeadCell{K,1:3}] = deal(FieldOrderStr,StructSIP.(FieldOrderStr),'SIP order');
                        
                        [MatSize] = size(Mat.(Fields{If}));
                        for Ix=1:1:MatSize(1)
                            for Iy=1:1:MatSize(2)
                                if ( Mat.(Fields{If})(Ix,Iy) ~= 0)
                                    FieldName = sprintf('%s_%d_%d',Fields{If},Ix-1,Iy-1);
                                    K = K + 1;
                                    [HeadCell{K,1:3}] = deal(FieldName, Mat.(Fields{If})(Ix,Iy), 'SIP distortion coef.');
                                end
                            end
                        end
                    end
                end
            end
            % format output
            switch lower(OutType)
                case 'head'
                    H = HEAD;
                    H = add_key(H,HeadCell);
                case 'cell'
                    H = HeadCell;
                otherwise
                    error('Unknown OutType option');
            end
        end

        
        % get SIP distortion parameters from header
        function StructSIP=get_sip(H)
            % Get SIP distortion keywords from header
            % Package: @ClassWCS
            % Description: Get the SIP distortion keywords from header.
            % Input  : - An header class object.
            % Output : - A structure array with the following fields
            %            .Exist - indicate if SIP exits
            %            .A     - A structure array of all SIP polynomials
            %                     of axis 1 for the detector to sky distortions.
            %            .B     - A structure array of all SIP polynomials
            %                     of axis 2 for the detector to sky distortions.
            %            .AP    - A structure array of all SIP polynomials
            %                     of axis 1 for the skey to detector distortions.
            %            .BP    - A structure array of all SIP polynomials
            %                     of axis 2 for the sky to detector distortions.
            %            .matrix - A structure array of A, B, AP, BP
            %                     polynomial coef in matrix form.
            %                    .matrix.A(I,J) corresponds to A_I_J
            %                    keyword.
            % Example: StructSIP = ClassWCS.get_sip(S1)

            Keys  = {'A','B','AP','BP'};
            Nkeys = numel(Keys);

%             Key_A_Order  = 'A_ORDER';   % sip polynomial order, axis 1, detector to sky
%             Key_B_Order  = 'B_ORDER';   % sip polynonial order, axis 2, detector to sky  
%             Key_AP_Order = 'AP_ORDER';   % sip polynomial order, axis 1, sky to detector
%             Key_BP_Order = 'BP_ORDER';   % sip polynonial order, axis 2, sky to detector

            KeyOrder = cell(1,Nkeys);
            for Ikeys=1:1:Nkeys
                KeyOrder{Ikeys} = sprintf('%s_ORDER',Keys{Ikeys});
            end

            Nh = numel(H);
            
            StructSIP = Util.struct.struct_def({'Exist', Keys{:}, 'matrix','A_ORDER','B_ORDER'},Nh,1);
            
            for Ih=1:1:Nh
                ValOrder = mgetkey(H,KeyOrder);
                for Ikeys=1:1:Nkeys
                    Order = ValOrder{Ikeys};

                    if (~isnan(Order))
                        % exist
                        StructSIP(Ih).Exist = true;
                        StructSIP(Ih).(KeyOrder{Ikeys}) = Order;
                        
                        Name = cell(Order+1, Order+1);
                        for Iorder1=0:1:Order
                            for Iorder2=0:1:Order
                                Name{Iorder1+1,Iorder2+1} = sprintf('%s_%d_%d',Keys{Ikeys},Iorder1,Iorder2);
                            end
                        end
                        Val = mgetkey(H,Name(:));
                        Val = cell2mat(Val);
                        Val(isnan(Val)) = 0;
                        CellVal = num2cell(Val);
                        StructSIP(Ih).(Keys{Ikeys}) = cell2struct(CellVal,Name(:),2);
                        StructSIP(Ih).matrix.(Keys{Ikeys}) = reshape(Val,Order+1,Order+1);
                    else
                        StructSIP(Ih).Exist = false;

                    end
                end
            end

        end
        
        
        
        
        % TranClass -> WCS
        function W=tranclass2wcs(TranC,varargin)
            % Convert TranClass object into a ClassWCS object
            % Package: @ClassWCS
            % Description: Convert a TranClass object into a ClassWCS
            %              object with the WCS transformation.
            % Input  : - A TranClass object
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'NormXY' - X/Y coordinates normalization factor.
            %                       Default is 1.
            %            'CooCenter' - Sky reference coordinates center [deg].
            %            'ImCenter'  - Image reference coordinates center [pix].
            %            'Scale'     - Pixel scale [arcsec/pix].
            % Output : - A ClassWCS object.
            
            
            WCSField = ClassWCS.WCSField;
            
            DefV.NormXY               = 1;
            DefV.CooCenter            = [0 0];
            DefV.ImCenter             = [0 0];
            DefV.Scale                = 1;
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            RAD = 180./pi;
            
            
            W = ClassWCS;
                        
            W.(WCSField).Status = true;  
            W.(WCSField).WCSAXES = 2;
            W.(WCSField).RADECSYS= 'ICRS';
            W.(WCSField).CRVAL   = InPar.CooCenter.*RAD;
            %W.(WCSField).CRPIX   = InPar.ImCenter + [Xshift, Yshift];
            W.(WCSField).CTYPE   = {'RA---TAN', 'DEC--TAN'};
            W.(WCSField).CUNIT   = {'deg', 'deg'};
            W.(WCSField).CRTYPE  = {'deg', 'deg'};
            %W.(WCSField).CD      = [Xrotx Xroty; Yrotx Yroty].*InPar.Scale./(InPar.NormXY.*3600);
            
            [PolyS1,OutOrder1,OutCoef1,ChPoly1]=tranclass2poly(TranC,1,true,InPar.NormXY);
            [PolyS2,OutOrder2,OutCoef2,ChPoly2]=tranclass2poly(TranC,2,true,InPar.NormXY);
            % Note the parameters normalization is done in fit_transform
            %[PolyS1,OutOrder1,OutCoef1,ChPoly1]=tranclass2poly(TranC,1,true,1);
            %[PolyS2,OutOrder2,OutCoef2,ChPoly2]=tranclass2poly(TranC,2,true,1);
            
            % find and store the shift parameters
            Fsx = all(OutOrder1==0);
            Fsy = all(OutOrder2==0);
            %ImCenter = InPar.ImCenter + [OutCoef1(Fsx), OutCoef2(Fsy)];
            W.(WCSField).CRPIX = InPar.ImCenter + [OutCoef1(Fsx), OutCoef2(Fsy)].*InPar.NormXY; %.*RAD.*3600./InPar.Scale;
            
            % remove the parameters from OutOrder/OutCoef not to confuse
            % next step
            OutOrder1 = OutOrder1(:,~Fsx);
            OutCoef1  = OutCoef1(~Fsx);
            OutOrder2 = OutOrder2(:,~Fsy);
            OutCoef2  = OutCoef2(~Fsy);
            
            
            
            % find and store the rotation parameters
            
            Fxx = OutOrder1(1,:)==1 & OutOrder1(2,:)==0;
            Fxy = OutOrder1(1,:)==0 & OutOrder1(2,:)==1;
            Fyx = OutOrder2(1,:)==1 & OutOrder2(2,:)==0;
            Fyy = OutOrder2(1,:)==0 & OutOrder2(2,:)==1;
            
            % Need to atke care of FunM*
            CD  = [OutCoef1(Fxx), OutCoef1(Fxy); OutCoef2(Fyx), OutCoef2(Fyy)];
            CD  = CD.*InPar.NormXY;
            
            % remove the parameters from OutOrder/OutCoef not to confuse
            % next step
            OutOrder1 = OutOrder1(:,~Fxx & ~Fxy);
            OutCoef1  = OutCoef1(~Fxx & ~Fxy);
            OutOrder2 = OutOrder2(:,~Fyy & ~Fyx);
            OutCoef2  = OutCoef2(~Fyy & ~Fyx);
            
            W.(WCSField).CD      = CD.*InPar.Scale./3600;
            
            if (~isempty(OutCoef1) || ~isempty(OutCoef2))
                % The rest of the coef. are in units of pixels
                % Note the minus sign
                OutCoef1 = -OutCoef1.*InPar.NormXY;  %.*InPar.Scale./3600;
                OutCoef2 = -OutCoef2.*InPar.NormXY;  %.*InPar.Scale./3600;
                % find and store the high order polynomials
                StructSIP = ClassWCS.poly2sip({OutOrder1, OutOrder2},{OutCoef1, OutCoef2});
                W.(WCSField).sip = StructSIP;

                W.(WCSField).CTYPE   = {'RA---TAN-SIP', 'DEC--TAN-SIP'};
            end
            
        end
        
        % TranClass -> WCS
        function W=tranclass2wcs_tpv(TranC,varargin)
            % Convert TranClass object into a ClassWCS object
            % Package: @ClassWCS
            % Description: Convert a TranClass object into a ClassWCS
            %              object with the WCS transformation.
            % Input  : - A TranClass object
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'NormXY' - X/Y coordinates normalization factor.
            %                       Default is 1.
            %            'CooCenter' - Sky reference coordinates center [deg].
            %            'ImCenter'  - Image reference coordinates center [pix].
            %            'Scale'     - Pixel scale [arcsec/pix].
            % Output : - A ClassWCS object.
            
            
            WCSField = ClassWCS.WCSField;
            
            DefV.NormXY               = 1;
            DefV.CooCenter            = [0 0];
            DefV.ImCenter             = [0 0];
            DefV.Scale                = 1;
            DefV.CD                   = [];
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            RAD = 180./pi;
            
            
            W = ClassWCS;
            W.(WCSField).Status = true;              
            W.(WCSField).WCSAXES = 2;
            W.(WCSField).RADECSYS= 'ICRS';
            W.(WCSField).CRVAL   = InPar.CooCenter.*RAD;
            W.(WCSField).CRPIX = InPar.ImCenter;
            %W.(WCSField).CRPIX   = InPar.ImCenter + [Xshift, Yshift];
            W.(WCSField).CTYPE   = {'RA---TPV', 'DEC--TPV'};
            W.(WCSField).CUNIT   = {'deg', 'deg'};
            W.(WCSField).CRTYPE  = {'deg', 'deg'};
            W.(WCSField).CD      = InPar.CD;
            %W.(WCSField).CD      = [Xrotx Xroty; Yrotx Yroty].*InPar.Scale./(InPar.NormXY.*3600);
            
            [PolyS1,OutOrder1,OutCoef1,ChPoly1]=tranclass2poly(TranC,1,true,1); %InPar.NormXY);
            [PolyS2,OutOrder2,OutCoef2,ChPoly2]=tranclass2poly(TranC,2,true,1); %InPar.NormXY);
            % Note the parameters normalization is done in fit_transform
            %[PolyS1,OutOrder1,OutCoef1,ChPoly1]=tranclass2poly(TranC,1,true,1);
            %[PolyS2,OutOrder2,OutCoef2,ChPoly2]=tranclass2poly(TranC,2,true,1);
            
            StructTPV = ClassWCS.poly2tpv({OutOrder1, OutOrder2},{OutCoef1, OutCoef2});
            
            W.(WCSField).tpv = StructTPV;
            
        end
        
        
        
        % Polynomials -> SIP structure
        function StructSIP=poly2sip(CellOrder,CellCoef)
            % Store polynomial orders and coeficents in SIP structure
            % Package: @ClassWCS
            % Description: Given a cell array of polynomial orders, and
            %              polynomial coef., return a structure with the
            %              SIP distortions parameters.
            % Input  : - Cell array of orders. Up to four elements cell
            %            array for the 'A', 'B', 'AP', BP' polynomials.
            %            Each cell contains a two rows matrix. First row
            %            for X powers, and 2nd rows for powers in Y.
            %          - Cell arrat of coef. (corresponding to orders).
            %            Each cell contains a vector of coef per columns in
            %            corresponding order matrix.
            % Exammple: CellOrder = {[1 2 3;2 1 1],[1 2 2 3; 1 2 3 3]};
            %           CellCoef={[1 1 1],[1 2 3 4]};
            %           StructSIP=ClassWCS.poly2sip(CellOrder,CellCoef)
            Keys = {'A','B','AP','BP'};
            Nkey = numel(CellOrder);
            
            StructSIP.Exist = true;
            for Ikey=1:1:Nkey
                Order    = CellOrder{Ikey};
                KeyOrder = sprintf('%s_ORDER',Keys{Ikey});
                % max order supplied
                StructSIP.(KeyOrder) = max(Order(:));
                
                Val = zeros(StructSIP.(KeyOrder)+1,StructSIP.(KeyOrder)+1);
                for Xp=0:1:StructSIP.(KeyOrder)
                    for Yp=0:1:StructSIP.(KeyOrder)
                        Ixy = find(Order(1,:)==Xp & Order(2,:)==Yp);
                        if (~isempty(Ixy))
                            Val(Xp+1,Yp+1) = sum(CellCoef{Ikey}(Ixy));
                        end
                        
                        KeyName = sprintf('%s_%d_%d',Keys{Ikey},Xp,Yp);
                        StructSIP.(Keys{Ikey}).(KeyName) = Val(Xp+1,Yp+1);
                    end
                end
                StructSIP.matrix.(Keys{Ikey}) = Val;
            end
            
            
        end
        
        function StructTPV=poly2tpv(CellOrder,CellCoef)
            % Store polynomial orders and coeficents in TPV structure
            % Package: @ClassWCS
            % Input  : - A cell array of matrices of polynomials orders.
            %            A cell array element per axis.
            %            The matrix contains two rows for X and Y orders.
            %          - A cell array of vectors of polynomials coeficents.
            %            A cell array element per axis.
            % Output : - A structure containing the TPV (PV) coeficents.
            
            
            [Tab{1},Tab{2}] = ClassWCS.tpv_polydef;
            Naxes = numel(CellOrder);
            
            K = 0;
            for Iaxis=1:1:Naxes
                
                % number of coef/orders
                N{Iaxis} = numel(CellCoef{Iaxis});
            
                for I=1:N{Iaxis}
                    Orders  = CellOrder{Iaxis}(:,I);
                    % requirement also that r-comp is 0.
                    IndexPV = find(all(Tab{Iaxis}(1:2,:)==Orders) & Tab{Iaxis}(3,:)==0) - 1;

                    if (isempty(IndexPV))
                        % skip
                    else
                        K = K + 1;
                        %StructTPV.KeyName{K} = sprintf('PV%d_%d',Iaxis,IndexPV);
                        %StructTPV.KeyVal{K}  = sum(CellCoef{Iaxis}(I));
                        
                        StructTPV.Ind(K,1:2) = [Iaxis,IndexPV];
                        StructTPV.KeyVal{K}  = sum(CellCoef{Iaxis}(I));
                        StructTPV.KeyName{K} = sprintf('PV%d_%d',Iaxis,IndexPV);
                    end
                end
            end
        end
        
        
        
        
    end
    
    % get fields/values
    methods
        % get CTYPE
        function CType=get_ctype(W,Index)
            % Get CTYPE header keyword from a ClassWCS object for one of the dimensions
            % Package: @ClassWCS
            % Input  : - A ClassWCS object.
            %          - The keyword dimension index. Default is 1.
            % Output : - A strcture array with one element per ClassWCS
            %            element, containing the following fields:
            %            .KeyVal
            %            .CooName
            %            .ProjName
            %            .IsSIP
            % Example: CType=get_ctype(W);
            
            WCSField = ClassWCS.WCSField;
            
            if (nargin<2)
                Index = 1;
            end
            Nw = numel(W);
            
            for Iw=1:1:Nw
                % Full CTYPE keyword value
                CType(Iw).KeyVal  = W(Iw).(WCSField).CTYPE{Index};
                % CTYPE coordinate
                CType(Iw).CooName = W(Iw).(WCSField).CTYPE{Index}(1:4);
                CType(Iw).CooName = strrep(CType(Iw).CooName,'-','');
                % CTYPE projection
                CType(Iw).ProjName = W(Iw).(WCSField).CTYPE{Index}(6:8);
                
                % is SIP
                CType(Iw).IsSIP = false;
                if (numel(W(Iw).(WCSField).CTYPE{Index})==12)
                    if (strcmp(W(Iw).(WCSField).CTYPE{Index}(10:12),'SIP'))
                        CType(Iw).IsSIP = true;
                    end
                end
                
                
            end
            
        end
        
        % get PV coef.
        function [Ind,KeyVal,KeyName]=get_pv(W)
            % Get PV coeficients from ClassWCS object
            % Package: @ClassWCS
            % Input  : - A ClassWCS object.
            % Output : - A two column matrix of PV indices.
            %          - A cell array of PV keyword values.
            %          - A cell array of PV keyword names.
            % Exampe: [Ind,KeyVal,KeyName]=get_pv(W);
            
            WCSField = ClassWCS.WCSField;
            
            if (numel(W)>1)
                error('ClassWCS object should contain a single element');
            end
            
            Iw = 1;
            if (isfield(W(Iw).(WCSField),'PV'))
                if (~isempty(W(Iw).(WCSField).PV))
                    Ind     = W(Iw).(WCSField).PV.Ind;
                    KeyVal  = W(Iw).(WCSField).PV.KeyVal;
                    KeyName = W(Iw).(WCSField).PV.KeyName;
                else
                    Ind     = [];
                    KeyVal  = {};
                    KeyName = {};
                end
            else
                Ind     = [];
                KeyVal  = {};
                KeyName = {};
            end
        end
        
        % get TPV polynomials
        function [SX,SY,Xi,Yi]=get_tpv(W,X,Y)
            % Get TPV polynomial coeficients and apply to coordinates
            % Package: @ClassWCS
            % Input  : - A single element ClassWCS object.
            %          - Intermidate world coordinates in X axis in the
            %            native units of the transformation (i.e., the 
            %            CUNIT keyword; e.g. 'deg'). Default is [].
            %          - Intermidate world coordinates in Y axis in the
            %            native units of the transformation (i.e., the 
            %            CUNIT keyword; e.g. 'deg'). Default is [].
            % Output : - A structure for the X axis orders and coef.
            %            Orders include three columns [X order, Y order, R
            %            order].
            %          - A structure for the X axis orders and coef.
            %          - Pixel coordinate X
            %          - Pixel coordinate Y
            % Example: [SX,SY]=get_tpv(W);
            %          [SX,SY,Xi,Yi]=get_tpv(W,rand(4,1),rand(4,1));
            
            WCSField = ClassWCS.WCSField;
            
            if (nargin<3)
                X = [];
                Y = [];
            end
            
            if (numel(W)>1)
                error('ClassWCS object should contain a single element');
            end
            Iw = 1;
            [Ind,KeyVal,~]=get_pv(W(Iw));
            
            if (isempty(Ind) || isempty(KeyVal))
                % Ind is empty
                SX.Orders = [];
                SX.Coef   = [];
                SY.Orders = [];
                SY.Coef   = [];
                Xi        = X;
                Yi        = Y;
            else
                % PV coef are defined
                FlagX = Ind(:,1)==1;
                FlagY = Ind(:,1)==2;

                [TabX,TabY] = ClassWCS.tpv_polydef;

                PV_IndexX = Ind(FlagX,2);
                PV_IndexY = Ind(FlagY,2);

                SX.Orders = TabX(:,PV_IndexX+1).';
                SX.Coef   = cell2mat(KeyVal(FlagX));

                SY.Orders = TabY(:,PV_IndexY+1).';
                SY.Coef   = cell2mat(KeyVal(FlagY));

                if (nargout>2)
                    R  = sqrt(X.^2 + Y.^2);
                    Xi = sum(SX.Coef(:) .* X(:).'.^SX.Orders(:,1) .* Y(:).'.^SX.Orders(:,2) .* R(:).'.^SX.Orders(:,3)).';
                    Yi = sum(SY.Coef(:) .* X(:).'.^SY.Orders(:,1) .* Y(:).'.^SY.Orders(:,2) .* R(:).'.^SY.Orders(:,3)).';
                end
            end
        end
        
    end
    
    % coordinate conversions
    methods
        function [RA,Dec]=xy2coo(W,XY,varargin)
            % Given ClassWCS object, convert X/Y pixel coordinates to RA/Dec
            % Package: @ClassWCS
            % Description: Convert X/Y pixel coordinates to RA/Dec using a
            %              ClassWCS object.
            %              If a single output argument is requested then
            %              the output is an AstCat object.
            %              If two output arguments are requested then the
            %              each output is a vector.
            % Input  : - A ClassWCS object.
            %            If this is a single element object then the
            %            input X/Y coordinates must be an AstCat object
            %            with the same number of elements.
            %          - A vector of X coordinates or an AstCat object.
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'OutUnits' - Output units. Default is 'rad'.
            %            'ColNameX' - Cell array of X coordinate column
            %                         names. If the second input is an AstCat object 
            %                         then these are X coordinate column
            %                         names that program will search in the
            %                         AstCat object. The column names that
            %                         appear first will be selected.
            %                         Default is
            %                         {'XWIN_IMAGE','X','X_IMAGE'}.
            %            'ColNameY' - Like ColNameX, but for the Y
            %                         coordinates. Default is
            %                         {'YWIN_IMAGE','Y','Y_IMAGE'}.
            %            'ColNameRA'- RA column name in the output AstCat
            %                         object. Default is 'RA'.
            %            'ColNameDec'-Dec column name in the output AstCat
            %                         object. Default is 'Dec'.
            % Output : - Vector of RA, or an AstCat object with [RA, Dec].
            %          - Vector of Dec.
            % Example: [RA,Dec]=xy2coo(W,[300,2700],'OutUnits','deg');
            %          [AC]=xy2coo(W,[300,2700]);
            
            DefV.OutUnits             = 'rad';
            DefV.ColNameX             = {'XWIN_IMAGE','X','X_IMAGE'};
            DefV.ColNameY             = {'YWIN_IMAGE','Y','Y_IMAGE'};
            DefV.ColNameRA            = 'RA';
            DefV.ColNameDec           = 'Dec';
            
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            %WCSField     = ClassWCS.WCSField;
            CatField     = AstCat.CatField;
            ColCellField = AstCat.ColCellField;
            
            if (AstCat.isastcat(XY))
                AC = XY;
            else
                AC = AstCat.array2astcat(XY,{InPar.ColNameX{1}, InPar.ColNameY{1}});
            end

                
            Nac = numel(AC);
            Nw  = numel(W);
            if ~((Nw==1 && Nac>0) || (Nac==1 && Nw>0) || (Nw==Nac))
                error('Size of ClassWCS and AstCat objects are not compatible');
            end
            
            N = max(Nw,Nac);
            if (nargout==1)
                OutAC = AstCat(N,1);
            end
            for I=1:1:N
                Iw  = max(Nw,I);
                Iac = max(Nac,I);
                % for each ClassWCS element
                
                
                
                [~,ColX] = select_exist_colnames(AC(Iac),InPar.ColNameX.');
                [~,ColY] = select_exist_colnames(AC(Iac),InPar.ColNameY.');
                
                PixCoo = AC(Iac).(CatField)(:,[ColX, ColY]);
                Nsrc   = size(PixCoo,1);
                
                
                if (~W(Iw).WCS.Status)
                    % ClassWCS is not containing relevant info
                    RA  = nan(Nsrc,1);
                    Dec = nan(Nsrc,1);
                
                else
                    
                    % Applay shift (via CRPIX) and rotation (via CD matrix)
                    % Include TPV distortions
                    Xi = pixel2intermediate(W(Iw),PixCoo); 

                    % Applay distortions (e.g., PV)
                    %Xd = pixel_distortion

                    % Applay sky projection
                    Angle = intermediate2native(W(Iw),Xi,'rad');

                    % Applay sky rotation
                    % native2celestial
                    % TBD

                    RA  = Angle(:,1);
                    Dec = Angle(:,2);

                    ConvFactor = convert.angular('rad',InPar.OutUnits);
                    RA         = RA.*ConvFactor;
                    Dec        = Dec.*ConvFactor;
                end
                
                if (nargout==1)
                    OutAC(I).(CatField)     = [RA, Dec];
                    OutAC(I).(ColCellField) = {InPar.ColNameRA, InPar.ColNameDec};
                    OutAC(I) = colcell2col(OutAC(I));
                end
                
            end
            if (nargout==1)
                RA = OutAC;
            end
            
        end
        
        
           function [X,Y]=coo2xy(W,RADec,varargin)
            % Given ClassWCS object, convert RA/Dec to pixel coordinates
            % Package: @ClassWCS
            % Description: Convert RA/Dec to X/Y pixel coordinates using a
            %              ClassWCS object.
            %              If a single output argument is requested then
            %              the output is an AstCat object.
            %              If two output arguments are requested then the
            %              each output is a vector.
            % Input  : - A ClassWCS object.
            %            If this is a single element object then the
            %            input X/Y coordinates must be an AstCat object
            %            with the same number of elements.
            %          - A two column vector of RA/Dec coordinates or
            %            an AstCat object.
            %          * Arbitrary number of pairs of arguments: ...,keyword,value,...
            %            where keyword are one of the followings:
            %            'InUnits'  - Input units. Default is 'rad'.
            %            'ColNameRA'- Cell array of X coordinate column
            %                         names. If the second input is an AstCat object 
            %                         then these are RA coordinate column
            %                         names that program will search in the
            %                         AstCat object. The column names that
            %                         appear first will be selected.
            %                         Default is
            %                         {'ALPHAWIN_J2000','RA'}.
            %            'ColNameDec'- Like ColNameX, but for the Dec
            %                         coordinates. Default is
            %                         {'DELTAWIN_J2000','Dec'}.
            %            'ColNameX' - X column name in the output AstCat
            %                         object. Default is 'X'.
            %            'ColNameY' - Y column name in the output AstCat
            %                         object. Default is 'Y'.
            % Output : - Vector of X, or an AstCat object with [X, Y].
            %          - Vector of Y.
            % Example: [X,Y]=coo2xy(W,[1,1],'InUnits','deg');
            %          [AC]=coo2xy(W,[1,1]);
            
            DefV.InUnits              = 'rad';
            DefV.ColNameRA            = {'ALPHAWIN_J2000','RA'};
            DefV.ColNameDec           = {'DELTAWIN_J2000','Dec'};
            DefV.ColNameX             = 'X';
            DefV.ColNameY             = 'Y';
            
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            %WCSField     = ClassWCS.WCSField;
            CatField     = AstCat.CatField;
            ColCellField = AstCat.ColCellField;
            
            if (AstCat.isastcat(RADec))
                AC = RADec;
            else
                AC = AstCat.array2astcat(RADec,{InPar.ColNameRA{1}, InPar.ColNameDec{1}});
            end

                
            Nac = numel(AC);
            Nw  = numel(W);
            if ~((Nw==1 && Nac>0) || (Nac==1 && Nw>0) || (Nw==Nac))
                error('Size of ClassWCS and AstCat objects are not compatible');
            end
            
            N = max(Nw,Nac);
            if (nargout==1)
                OutAC = AstCat(N,1);
            end
            for I=1:1:N
                Iw  = max(Nw,I);
                Iac = max(Nac,I);
                % for each ClassWCS element
                
                [~,ColRA]  = select_exist_colnames(AC(Iac),InPar.ColNameRA(:));
                [~,ColDec] = select_exist_colnames(AC(Iac),InPar.ColNameDec(:));
                
                Angle = AC(Iac).(CatField)(:,[ColRA, ColDec]);
                
                % Apply inverse sky projection
                Xi    = native2intermediate(W(Iw),Angle,InPar.InUnits);
                
                % Applay shift (via CRPIX) and rotation (via inv(CD) matrix)
                % Include inverse TPV distortions
                P = intermediate2pixel(W(Iw),Xi);
                X = P(:,1);
                Y = P(:,2);
                
                if (nargout==1)
                    OutAC(I).(CatField)     = [X, Y];
                    OutAC(I).(ColCellField) = {InPar.ColNameX, InPar.ColNameY};
                    OutAC(I) = colcell2col(OutAC(I));
                end
                
            end
            if (nargout==1)
                X = OutAC;
            end
            
        end
        
        
    end
    
    
    % transformations
    methods 
        %
        function ProjName=get_projection(W)
            % Get projection name from ClassWCS object
            % Package: @ClassWCS
            % Description: Get projection name from ClassWCS object
            % Input  : - A ClassWCS object.
            % Output : - A cell array of projection names. One per ClassWCS
            %            object element.
            % Example: 
            WCSField = ClassWCS.WCSField;
            
            Nw = numel(W);
            ProjName = cell(Nw,1);
            for Iw=1:1:Nw
                ProjName{Iw} = W(Iw).(WCSField).CTYPE{1}(6:8);
            end
            
        end
        
        function [Phi0,Theta0]=get_native_fiducial(W)
            %
            
            WCSField = ClassWCS.WCSField;
            
            Nw = numel(W);
            for Iw=1:1:Nw
                if (isfield(W(Iw).(WCSField),'PV1_1'))
                    Phi0   = W(Iw).(WCSField).PV1_1;
                    Theta0 = W(Iw).(WCSField).PV1_2;
                else
                    Phi0 = 0
                end
            end
            
        end
        
        function [LonPole,LatPole]=get_lonlatpole(W,OutUnits)
            % Get LONPOLE/LATPOLE from ClassWCS object
            % Package: @ClassWCS
            % Description: Get LONPOLE/LATPOLE from ClassWCS object.
            %              Default is defined for the 'TAN' projection.
            % Input  : - A ClassWCS object
            %          - Output units. Default is 'rad'.
            % Output : - A vector of LONPOLE.
            %          - A vector of LATPOLE.
            % Exaple: [LonPole,LatPole]=get_lonlatpole(W);
            
            if (nargin<2)
                OutUnits = 'rad';
            end
            
            WCSField = ClassWCS.WCSField;
            LonPoleField = 'LONPOLE';
            LatPoleField = 'LATPOLE';
            
            Nw = numel(W);
            LonPole = zeros(Nw,1);
            LatPole = zeros(Nw,1);
            for Iw=1:1:Nw
                %CRVal = W(Iw).(WCSField).CRVAL;
                
                if (isfield(W(Iw).(WCSField),LonPoleField))
                    LonPole(Iw) = W(Iw).(WCSField).(LonPoleField);
                else
                    LonPole(Iw) = convert.angular('deg',W(Iw).(WCSField).CUNIT{1},180);
                end
                if (isfield(W(Iw).(WCSField),LatPoleField))
                    LatPole(Iw) = W(Iw).(WCSField).(LatPoleField);
                else
                    LatPole(Iw) = convert.angular('deg',W(Iw).(WCSField).CUNIT{2},0);
                end
                
                % convert to output units
                LonPole(Iw) = convert.angular(W(Iw).(WCSField).CUNIT{1},OutUnits,LonPole(Iw));
                LatPole(Iw) = convert.angular(W(Iw).(WCSField).CUNIT{2},OutUnits,LatPole(Iw));
            end
            
            
        end
        
        function CUnit=get_cunit(W)
            % Get CUNIT values from ClassWCS object
            % Package: @ClassWCS
            % Input  : - A ClassWCS object
            % Output : - A cell array of CUNIT values
            
            WCSField = ClassWCS.WCSField;
            
            Nw = numel(W);
            if (Nw>1)
                error('get_cunit is defined only for a single element ClassWCS object');
            end
            for Iw=1:1:Nw
                CUnit = W(Iw).(WCSField).CUNIT;
            end
                
        end
        
        %
        function XY=pixel2intermediate(W,PixCoo)
            % Convert pixel coordinates to intermediate WCS coordinates.
            % Description: Convert pixels coordinates (p_j) to intermediate
            %              wcs (projection plane coordinates).
            %              This is achived by applying the ...
            % Package: @ClassWCS
            % Input  : - An ClassWCS object containing the WCS parameters.
            %          - A two column matrix containing the [X, Y] pixel
            %            coordinates in the image.
            % Output : - A two column matrix containing the projection
            %            plane coordinates (intermediate WCS coordinates).
            %            The units of these coordinates are as specified in
            %            the CUNIT WCS keyword parameter (usually 'deg').
            % Example: X=pixel2intermediate(W,[0 0;1 1;2 2])
            
            if (numel(W)>1)
                error('ClassWCS object should contain a single element');
            end
            
            Iw = 1;
            
            WCSField = ClassWCS.WCSField;
            XY = (W(Iw).(WCSField).CD * (PixCoo - W(Iw).(WCSField).CRPIX).').';
            
            % Applay TAN/PV distortions if relevant
            CType = get_ctype(W(Iw));
            
            switch lower(CType.ProjName)
                case 'tpv'
                    % applay TAN/PV distortions
                    [~,~,XY(:,1),XY(:,2)] = get_tpv(W(Iw),XY(:,1),XY(:,2));
                    
                otherwise
                    % do nothing
            end
            
        end
        
        function [X,Y]=inverse_tpv(W,Xt,Yt)
            % Applay the inverse tpv polynomial transformation
            % Package: @ClassWCS
            % Description: Applay the inverse tpv polynomial transformation
            % Input  : - A single element ClassWCS object.
            %          - Target X coordinate.
            %          - Target Y coordinate.
            % Output : - X coordinate of the inverse transformation.
            %          - Y coordinate of the inverse transformation.
            % Example: [~,~,Xt,Yt] = get_tpv(W,1,1);
            %          [X,Y] = inverse_tpv(W,Xt,Yt)
            ThreshConv = 1e-12;
            
            if (numel(W)>1)
                error('ClassWCS object should contain a single element');
            end
            
            Iw = 1;
            
            Xt_O = Xt;
            Yt_O = Yt;
            
            [~,~,Xtt_1,Ytt_1] = get_tpv(W(Iw),Xt,Yt);
            
            
            DX_1 = Xt - Xtt_1;
            DY_1 = Yt - Ytt_1;
            Conv = false;
            while (~Conv)
                
                X_1  = Xt + DX_1;
                Y_1  = Yt + DY_1;

                [~,~,Xt_1,Yt_1] = get_tpv(W(Iw),X_1,Y_1);
                DX_1 = Xt_O - Xt_1;
                DY_1 = Yt_O - Yt_1;
                %[DX_1, DY_1]
                Xt = X_1;
                Yt = Y_1;
                if (max(abs(Xt_O-Xt_1))<ThreshConv && max(abs(Yt_O-Yt_1))<ThreshConv)
                    Conv = true;
                end
                
            end
            X = X_1;
            Y = Y_1;
            
        end
        
        function PixCoo=intermediate2pixel(W,XY)
            % Convert intermediate WCS coordinates to pixel coordinates.
            % Description: Convert intermediate wcs (projection plane
            %              coordinates) to pixels coordinates (p_j).
            % Package: @ClassWCS
            % Input  : - An ClassWCS object containing the WCS parameters.
            %          - A two column matrix containing the projection
            %            plane coordinates (intermediate WCS coordinates).
            %            The units of these coordinates are as specified in
            %            the CUNIT WCS keyword parameter (usually 'deg').
            % Output : - A two column matrix containing the [X, Y] pixel
            %            coordinates in the image.
            % Example: P=intermediate2pixel(W,[0 0;1 1;2 2])
            
            if (numel(W)>1)
                error('ClassWCS object should contain a single element');
            end
            
            %error('in progress')
            
            Iw = 1;
            
            WCSField = ClassWCS.WCSField;
            
            % Applay TAN/PV distortions if relevant
            CType = get_ctype(W(Iw));
            
            switch lower(CType.ProjName)
                case 'tpv'
                    % applay TAN/PV distortions
                    % invers tpv function
                    % regular: [~,~,XY(:,1),XY(:,2)] = get_tpv(W(Iw),XY(:,1),XY(:,2));
                    
                    % inverse TPV...
                    [X,Y] = inverse_tpv(W,XY(:,1),XY(:,2));
                    XY = [X(:), Y(:)];
                    
                otherwise
                    % do nothing
            end
            
            
            PixCoo = [inv(W(Iw).(WCSField).CD) * XY.' + W(Iw).(WCSField).CRPIX(:)].';
            
            %XY = (W(Iw).(WCSField).CD * (PixCoo - W(Iw).(WCSField).CRPIX).').';
            
            
            
        end
        
        function Angle=intermediate2native(W,Xi,Output)
            % Apply spherical projection to intermediate coordinates
            % Description: Convert intermediate pixel coordinates to native
            %              spherical coordinates by applying the projection
            %              transformation.
            
            % Example: Angle=intermediate2native(W,X)
            
            if (nargin<3)
                Output = 'rad';
            end
            
            RAD = 180./pi;
            
            WCSField = ClassWCS.WCSField;
            
            % number of columns in the inpput coordinates - i.e.,
            % dimensionality
            Ncol = size(Xi,2);
            
            ProjName=get_projection(W);
            Nw = numel(W);
            if (Nw>1)
                error('native2celestial is defined only for a single element ClassWCS object');
            end
            for Iw=1:1:Nw
                % units
                Xr = zeros(size(Xi));
                CRVal = zeros(1,Ncol);
                for Icol=1:1:Ncol
                    Xr(:,Icol)  = convert.angular(W(Iw).(WCSField).CUNIT{Icol},'rad',Xi(:,Icol));
                    CRVal(Icol) = convert.angular(W(Iw).(WCSField).CUNIT{Icol},'rad',W(Iw).(WCSField).CRVAL(Icol));
                end
                
                
                
                switch lower(ProjName{Iw})
                    case {'tan','tpv'}
                        % Tangential (Gnomonic) projection
                        % is zenithal perspective with gamma=0, mu=0
                        X = Xr(:,1);
                        Y = Xr(:,2);
                        
                        [Long,Lat]=celestial.proj.pr_ignomonic(X,Y,CRVal);
                        %[Long,Lat]=celestial.proj.pr_ignomonic(X,Y,[0 0]);
                        Angle = [Long(:), Lat(:)];
                        
%                         Gamma = 0;
%                         Mu    = 0;
%                         
%                         R = sqrt(X.^2 + (Y.*cos(Gamma)).^2);
%                         Rho = R./(RAD.*(Mu + 1) + Y.*sin(Gamma));
%                         Psi = atan2(Rho,1);
%                         Omega = asin( Rho.*Mu./sqrt(1+Rho.^2) );
%                         Theta = Psi - Omega;
%                         Phi   = atan2(-Y.*cos(Gamma), X);
%                         Theta = [Theta, Phi];
                    otherwise
                        error('Unknown Projection name');
                end
            end
                        
            % Convert to output units
            Angle = convert.angular('rad',Output,Angle);
            
        end
        
        function Xi=native2intermediate(W,Angle,Input)
            % Apply inverse spherical projection to native coordinates
            % Description: Convert native spherical coordinates to 
            %              intermediate pixel coordinates by applying the
            %              inverse projection transformation.
            
            % Example: Xi=native2intermediate(W,Angle,'rad')
                        
            if (nargin<3)
                % Input coordinates
                Input = 'rad';
            end
            
            RAD = 180./pi;
            
            WCSField = ClassWCS.WCSField;
            
            % number of columns in the inpput coordinates - i.e.,
            % dimensionality
            Ncol = size(Angle,2);
            
            ProjName = get_projection(W);
            Nw = numel(W);
            if (Nw>1)
                error('native2intermidiate is defined only for a single element ClassWCS object');
            end
            
            for Iw=1:1:Nw
                AngleRad = convert.angular(Input,'rad',Angle); % [rad]
                % units
                CRVal = zeros(1,Ncol);
                for Icol=1:1:Ncol
                    CRVal(Icol) = convert.angular(W(Iw).(WCSField).CUNIT{Icol},'rad',W(Iw).(WCSField).CRVAL(Icol));
                end
                
                switch lower(ProjName{Iw})
                    case {'tan','tpv'}
                        % Tangential (Gnomonic) projection
                        % is zenithal perspective with gamma=0, mu=0
                       
                        [X,Y] = celestial.proj.pr_gnomonic(AngleRad(:,1),AngleRad(:,2),RAD,CRVal);
                        
                    otherwise
                        error('Unknown Projection name');
                end
                
                Xi = [X(:), Y(:)];
                
                % units
%                 Xr = zeros(size(Xi));
%                 CRVal = zeros(1,Ncol);
%                 for Icol=1:1:Ncol
%                     Xr(:,Icol)  = convert.angular(W(Iw).(WCSField).CUNIT{Icol},'rad',Xi(:,Icol));
%                     CRVal(Icol) = convert.angular(W(Iw).(WCSField).CUNIT{Icol},'rad',W(Iw).(WCSField).CRVAL(Icol));
%                 end
                
            end
              
        end
        
        
        
        function Alpha=native2celestial(W,Angle)
            % DOES NOT WORK
            % Description: Convert native spherical coordinates to
            %              celestial spherical coordinates.
            % Example: Alpha=native2celestial(W,Theta)
            
            WCSField = ClassWCS.WCSField;
            
            
            % Get LONPOLE and LATPOLE key values
            [Phi_P,Theta_P] = get_lonlatpole(W,'rad');
            
            Nw = numel(W);
            if (Nw>1)
                error('native2celestial is defined only for a single element ClassWCS object');
            end
            for Iw=1:1:Nw
                Phi     = Angle(:,2);
                Theta   = Angle(:,1);
                
                CUnit   = get_cunit(W);
                CRVal   = W(Iw).(WCSField).CRVAL;
                Alpha_P = 0; %convert.angular(CUnit{1},'rad',CRVal(1));
                Delta_P = 0; %convert.angular(CUnit{2},'rad',CRVal(2));
                
                Alpha = Alpha_P + atan2(sin(Theta).*cos(Delta_P) - cos(Theta).*sin(Delta_P).*cos(Phi - Phi_P), -cos(Theta).*sin(Phi-Phi_P));
                Delta = asin( sin(Theta).*sin(Delta_P) + cos(Theta).*cos(Delta_P).*cos(Phi-Phi_P) );

                Alpha = [Alpha(:), Delta(:)];
            end
            
        end
    end
    
    methods
        % WCS to Header
        function H=wcs2head(W,H)
            % Convert ClassWCS object to header WCS key/par
            % Package: @ClassWCS
            % Description:
            % Input  : - A ClassWCS object.
            %          - Optional header in which to concat the WCS
            %            keywords.
            % Output : - An HEAD object with the WCS keywords.
            
            
            
            WCSField = ClassWCS.WCSField;
            % debug: W.WCS.CD = W.WCS.CD .* [1 -1;1 1];
            if (nargin<2)
                H = [];
            end
            if (iscell(H))
                Tmp = H;
                H = HEAD;
                H = add_key(H,Tmp);
            end
            if (isempty(H))
                H = HEAD(size(W));
            end
            
            Nw = numel(W);
            for Iw=1:1:Nw
                CellHead = cell(0,3);
                Iline    = 0;
                FN = fieldnames(W(Iw).(WCSField));
                Nfn = numel(FN);
                for Ifn=1:1:Nfn
                    if (isstruct(W(Iw).(WCSField).(FN{Ifn})))
                        % struct contains distortions or meta data
                        switch lower(FN{Ifn})
                            case 'sip'
                                AddCell = ClassWCS.sip2head(W(Iw).(WCSField).sip,'cell');
                                CellHead= [CellHead; AddCell];
                                Iline    = Iline + size(AddCell,1);
                                
                                if (size(AddCell,1)>1)
                                    IsSIP    = true;
                                else
                                    IsSIP    = false;
                                end
                                
                            case 'tpv'
                                Nc = numel(W.WCS.tpv.KeyVal);
                                AddCell = cell(Nc,3);
                                [AddCell{:,3}]=deal('');
                                [AddCell{:,1}] = deal(W(Iw).(WCSField).tpv.KeyName{:});
                                [AddCell{:,2}] = deal(W(Iw).(WCSField).tpv.KeyVal{:});
                                
                                CellHead= [CellHead; AddCell];
                                Iline    = Iline + Nc;
                            case 'pv'
                                Nc = numel(W.WCS.PV.KeyVal);
                                AddCell = cell(Nc,3);
                                [AddCell{:,3}]=deal('');
                                [AddCell{:,1}] = deal(W(Iw).(WCSField).PV.KeyName{:});
                                [AddCell{:,2}] = deal(W(Iw).(WCSField).PV.KeyVal{:});
                                
                                CellHead= [CellHead; AddCell];
                                Iline    = Iline + Nc;
                                
                            otherwise
                                warning('Unknown structure type in a ClassWCS object');
                        end
                        
                        
                    else
                        % N1 N2 are the dimensions of the keyword value element
                        % e.g., CD is usually a 2x2 matrix...
                        if (ischar(W(Iw).(WCSField).(FN{Ifn})))
                            N1 = 1;
                            N2 = 1;
                        else
                            [N1,N2] = size(W(Iw).(WCSField).(FN{Ifn}));
                        end
                        
                        if (N1>1 && N2>1)
                            % field contains a matrix
                            for I1=1:1:N1
                                for I2=1:1:N2
                                    KeyName = sprintf('%s%d_%d',FN{Ifn},I1,I2);
                                    KeyVal  = W(Iw).(WCSField).(FN{Ifn})(I1,I2);
                                    Iline   = Iline + 1;
                                    CellHead(Iline,:) = {KeyName, KeyVal, ''};
                                end
                            end
                        else
                            % field contains a vector/scalar
                            if (N1>1 || N2>1)
                                % field contains a vector
                                for I1=1:1:max(N1,N2)
                                    KeyName = sprintf('%s%d',FN{Ifn},I1);
                                    KeyVal  = W(Iw).(WCSField).(FN{Ifn})(I1);
                                    Iline   = Iline + 1;
                                    if (iscell(KeyVal))
                                        KeyVal = KeyVal{1};
                                    end
                                    CellHead(Iline,:) = {KeyName, KeyVal, ''};
                                end
                            else
                                % field contain a scalar
                                Iline   = Iline + 1;
                                CellHead(Iline,:) = {FN{Ifn}, W(Iw).(WCSField).(FN{Ifn}), ''};
                            end
                            
                        end
                        
                    end
                end

                % delete old WCS keywords
                H(Iw) = delete_key(H(Iw),CellHead(:,1));
                % add new keywords
                H(Iw) = add_key(H(Iw),CellHead);
                
            end
            
            % end of main loop (Iw)
        end
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
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
        
        
        function [XiTab,EtaTab]=tpv_polydef
            % Return the TPV transformation polynomial orders definitions
            % Package: @ClassWCS
            % Description: Return the TPV transformation polynomial orders
            %              definitions.
            % Input  : null
            % Output : - A three column matrix of the xi' polynomial order
            %            in xi, eta, r,
            %            where PV index 0 refers to the first table line.
            %          - A three column matrix of the eta' polynomial order
            %            in xi, eta, r,
            %            where PV index 0 refers to the first table line.
            % Example: [XiTab,EtaTab]=ClassWCS.tpv_polydef
            
            XiTab  = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 2 0 0; 1 1 0; 0 2 0; 3 0 0; 2 1 0; 1 2 0; 0 3 0; 0 0 3; 4 0 0; 3 1 0; 2 2 0; 1 3 0; 0 4 0; 5 0 0; 4 1 0; 3 2 0; 2 3 0; 1 4 0; 0 5 0; 0 0 5; 6 0 0; 5 1 0; 4 2 0; 3 3 0; 2 4 0; 1 5 0; 0 6 0; 7 0 0; 6 1 0; 5 2 0; 4 3 0; 3 4 0; 2 5 0; 1 6 0; 0 7 0; 0 0 7].';
            EtaTab = [0 0 0; 0 1 0; 1 0 0; 0 0 1; 0 2 0; 1 1 0; 2 0 0; 0 3 0; 1 2 0; 2 1 0; 3 0 0; 0 0 3; 0 4 0; 1 3 0; 2 2 0; 3 1 0; 4 0 0; 0 5 0; 1 4 0; 2 3 0; 3 2 0; 4 1 0; 5 0 0; 0 0 5; 0 6 0; 1 5 0; 2 4 0; 3 3 0; 4 2 0; 5 1 0; 6 0 0; 0 7 0; 1 6 0; 2 5 0; 3 4 0; 4 3 0; 5 2 0; 6 1 0; 7 0 0; 0 0 7].';
        
        end
    end
    
    % transformations
    methods (Static)
        
        function CD=read_cd_mtarix(WCS)
            
           %if ('WCSAXES'
            
        end
        
        
        
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
            
            
        function pix2native(P,Proj,PV) 
            %
            
            switch lower(Proj)
                case 'tan'
                    
                otherwise
                    error('unknown transformation type');
            end
                    
        end
        
        function distortion_sip_poly2keys(PolyOrder,PolyCoef)
            
        end
        
        
        
    end
    
    
end

            
