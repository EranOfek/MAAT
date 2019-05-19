%--------------------------------------------------------------------------
% MASK class                                                         class
% Description: A class of for an image mask
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef MASK
    properties (SetAccess = public)
        Mask
        MaskDic
    end
    
    methods

        %-------------------------
        %--- Class constructor ---
        %-------------------------
        
        function Mask=MASK(N,M)
            % MASK object constructor
            % Package: @MASK
            % Description: MASK constructor method
            % Input  : - Number of rows or [rows columns].
            %          - Number of columns.
            % Output : - A MASK object
            % Reliable: 2
            
            MaskField = 'Mask';
            
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
                    Mask(I,J).(MaskField) = [];
                end
            end
        end
        
        
    end
    
    % Static class - Fields
    methods (Static)
    
        function Ans=ismask(Obj)
            % Return true if object is MASK
            % Description: Check if object is of MASK class.
            % Input  : - Object
            % Output : - {true|false}.
            %     By : Eran O. Ofek                    Feb 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: S=1; MASK.ismask(S);
            % Reliable: 2
            Ans = isa(Obj,'MASK');
        end
            
        function Name=MaskField
            % Return the MASK class Mask field name
            Name = 'Mask';
        end
      
    end
      
    % Static class - bit mask dictionaries
    methods (Static)
        function [Bit,Map]=def_bitmask_pipeline(Field)
            % Bit mask dictionary
            % Description: The pipeline bit mask definition.
            %              Given the Bit mask name return the bit mask index.
            % Input  : - Bit mask name.
            % Output : - Bit index. If Bit mask is empty then return empty.
            %          - A structure array of the bit mask dictionary.
            % Tested : Matlab R2011b
            %     By : Eran O. Ofek                    Feb 2014
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % See also: maskflag_check.m, maskflag_set.m
            % Example: Bit=MASK.def_bitmask_pipeline('Bit_ImSaturated');
            %          [~,Map]=MASK.def_bitmask_pipeline
            % Reliable: 2
            %--------------------------------------------------------------------------

            if (nargin==0)
                Field = [];
            end

            Ind = 0;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_ImSaturated';   % 1     SIM/flag_saturated
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_BiasNoisy';     % 2     SIM/bias
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_BiasNoise0';    % 3     SIM/bias
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_BiasAnom';      % 4     SIM/bias
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_BiasFlare';     % 4     SIM/bias
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_BiasVariable';     % 4   SIM/flag_bias_variable
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_BadPixel';      % 5
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_FlatNaN';       % 6
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_FlatLowNim';    % 7
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_FlatHighStd';   % 8
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_FlatLow';       % 9 - flat low value (7-64) [in sim_flat.m]
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_StarFlatHighStd';    % 10
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_StarFlatAnom';       % 11
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_Hole';          % Hole (below background) in corrected image
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_CR';                % 12
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_CR_MEXTRACTOR';     % 13
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_CR_LACOS';           % 14
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_Edge';              % 15
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_Spike';            % 16
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_Ghost';             % 17
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_NonLin';        % 18
            Map(Ind).Ind  = Ind;
            
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_BackGrad';          % Large gradients in background
            Map(Ind).Ind  = Ind;
            
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_BackDiff';          % Background differ from ref image background
            Map(Ind).Ind  = Ind;
            
            %-----------------------------
            %--- sources related flags ---
            %-----------------------------

            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_SN30';     %  18
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_SN10';     % 19
            Map(Ind).Ind  = Ind;
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_SN5';       % 20
            Map(Ind).Ind  = Ind;
            % segmented sources (S/N>3)
            Ind = Ind + 1;
            Map(Ind).Name = 'Bit_ConnSrc';       % 21 
            Map(Ind).Ind  = Ind;


            if (isempty(Field))
                Bit = [];
            else

                Ind = find(strcmp({Map.Name},Field));
                Bit = Map(Ind).Ind;
            end
        end
        
        
    end
    
    % 
    methods (Static)
        function [Flag,Count]=isbit_flag(Mask,BitMask,Comparison,BitType,Dic)
            % Search for specific open bits in bitmask matrix
            % Package: @MASK
            % Description: Search for specific open bits in bitmask matrix
            % Input  : - A matrix or vector of bitmask flags.
            %          - A cell array of bit names as appear in the
            %            dictionary, or an exact bit mask value.
            %          - Comparison method:
            %            'bitand'   - Bit and operation. Default.
            %            'exactval' - By exact value.
            %          - Bit type. If the BitMask is numeric then this parameter is
            %            used to specify if the bit mask is a vector of indices
            %            ('index') or a bit mask value ('value').
            %            Default is 'value'.
            %          - Bit mask dictionary.
            %            Default is @MASK.def_bitmask_pipeline
            % Output : - A logical array indicating if bit mask was
            %            identified in each array element. 
            %          - Count of identified bit masks.
            % Example: [Flag,Count]=MASK.isbit_flag(Sm.Cat(:,40),{'Bit_CR'})
            
            
            Def.Comparison = 'bitand';  % 'bitand'|exactval' 
            Def.BitType    = 'Value';   % 'value'|'index'
            Def.Dic        = @MASK.def_bitmask_pipeline;
            if (nargin==2)
                Comparison = Def.Comparison;
                BitType    = Def.BitType;
                Dic        = Def.Dic;
            elseif (nargin==3)
                BitType    = Def.BitType;
                Dic        = Def.Dic;
            elseif (nargin==4)
                Dic        = Def.Dic;
            elseif (nargin==5)
                % do nothing
            else
                error('Input arguments: [MaskIm,BitMask,[Comparison,BitType]]');
            end


            
            
            if ischar(BitMask)
                BitMask = {BitMask};
            end
            
            if (iscell(BitMask))
                % get bit map value from names
                M = MASK;
                M.MaskDic  = Dic;
                BitMaskVal = bitname2val(MASK,BitMask);

            else
                switch lower(BitType)
                    case 'value'
                        % assume bit mask is provided by its decimal value
                        BitMaskVal = BitMask;
                    case 'index'
                        BitMaskVal = sum(bitset(0,BitMask));
                    otherwise
                        error('Unown BitType option');
                end
            end

            switch lower(Comparison)
                case 'exactval'
                    Flag = Mask == BitMaskVal;
                case 'bitand'
                    Flag = bitand(uint32(Mask),uint32(BitMaskVal),'uint32')>0;
                otherwise
                    error('Unknown Comparison option');
            end

            if (nargout>1)
                Count = sum(Flag(:));
            end
    
        end
    end
    
end
            
