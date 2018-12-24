% A class for stack objects
% Description: A class implementation of FIFO stack
%              .St contains the stack - first element is the first in,
%                  while the last element is the last in.
%              .Data - the data.
%              .StackType - stack type {'FIFO'}
% Tested : Matlab R2018a
%     By : Eran O. Ofek                    Nov 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 2
%--------------------------------------------------------------------------

classdef stack < handle
    properties (SetAccess = public)
        St
        StackType
        RecentEntry = 'last';     
        UserData
    end
    
    properties (Hidden)
        Data
        DataType
        Ind
        Size
        LastInd
        Full
        SI
    end

    % constructor
    methods
        
        function S=stack(Data,StackType)
            % stack constructor method
            % Package: @stack
            % Input  : - A cell, struct array,
            %            or array vector of some size.
            %            Default is nan(100,1).
            %          - Type {'fifo','filo'}. Default is 'fifo'
            
            
            Def.Data      = nan(100,1);
            Def.StackType = 'fifo';
            
            if (nargin<2)
                StackType = Def.StackType;
                if (nargin<1)
                    Data = Def.Data;
                end
            end
              
            S.StackType = StackType;
            if (iscell(Data))
                S.DataType = 'cell';
            else
                if (isstruct(Data))
                    S.DataType = 'struct';
                else
                    if (isnumeric(Data))
                        S.DataType = 'array';
                    else
                        error('Unknown Data Type');
                    end
                end
            end
            
            S.Data    = Data;
            S.Size    = numel(S.Data);
            S.LastInd = 0;
            S.Full    = false;
           
        end

    end
    
    % get/set
    methods
        function Array=get.St(S)
            % get St property from stack (i.e., sorted stack)
            
            sorted_ind(S);
            Array = S.Data(S.SI,:);
        end
    end
    
  
    % add
    methods
        function S=add(S,NewEl)
            % add new element into stack object
            
            if (S.LastInd<S.Size)
                S.LastInd = S.LastInd + 1;
                switch lower(S.DataType)
                    case {'array'}
                        S.Data(S.LastInd,:) = NewEl;
                    case 'cell'
                        S.Data{S.LastInd} = NewEl;
                    otherwise
                        error('Unknown DataType option');
                end
            else
                % FIFO
                switch lower(S.StackType)
                    case 'fifo'
                        S.Full    = true;
                        S.LastInd = 1;
                        switch lower(S.DataType)
                            case {'array','struct'}
                                S.Data(S.LastInd,:) = NewEl;
                            case 'cell'
                                S.Data{S.LastInd} = NewEl;
                            otherwise
                                error('Unknown DataType option');
                        end
                    otherwise
                        error('Unknown StackType option');
                end
            end
            % update SI
            S = sorted_ind(S);
            
        end
        
        
        function S=sorted_ind(S)
            % get vector of sorted indices (first in is 1)
            
            if (~S.Full)
                SIv = (1:1:S.LastInd)';
            else
                SIv = [(S.LastInd+1:S.Size), (1:S.LastInd)]';
            end
            S.SI = SIv;
        end
        
        
    end
    
 
end

            
