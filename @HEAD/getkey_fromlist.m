function BestVal=getkey_fromlist(Head,KeyList,Method,varargin)
% Search for the first existing haeder keyword in a list of keywords.
% Package: @HEAD
% Description: Get a specific keyword value from an HEAD object, searching
%              it in multiple keywords and returning the best result.
% Input  : - An HEAD object (e.g., a SIM).
%          - A string or a cell array of strings. Each string is a keyword
%            name. The function looks for these keywords in the header
%            and retrive their value. Then, the best value is selected.
%            Alternatively this can be a numeric value. In this case,
%            the numeric value will be return as is (one per HEAD element)
%            in a cell array.
%          - One of the following methods by which to select the keywords
%            value:
%            'first' - get the first not NaN value. Default.
%            'last'  - get the last not NaN value.
%            'max'   - get the max val.
%            'min'   - get the min val.
%            'mean'  - get the mean non NaN val.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'SpaceDel' - Delete spaces from key value. Default is true.
%            'Conv2Num' - A logical flag indicating if to attempt to convert the value
%                         to a number (if possible). Default is false.
% Output : - A cell vector of the best selected key values. One per header.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: BestVal = getkey_fromlist(Head,{'IMGTYPE','IMTYPE','IMGTYP'});
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin<3)
    Method = 'first';
end

DefV.SpaceDel          = true;
DefV.Conv2Num          = false;
if (numel(varargin)>0)
    InPar = InArg.populate_keyval(DefV,varargin,mfilename);
    %InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
else
    InPar = DefV;
end

Nh      = numel(Head);
if (isnumeric(KeyList))
    % Key is already given as a numeric value - return as is
    BestVal = KeyList;
    if (numel(BestVal)==1)
        % make the same length as Head
        BestVal = ones(Nh,1).*BestVal;
    else
        % convert to a column vector
        BestVal = BestVal(:);
    end
    % convert the vector into a cell array so the output will be always
    % a cell array.
    BestVal = num2cell(BestVal);
else
    Val     = mgetkey(Head,KeyList,InPar.SpaceDel,InPar.Conv2Num); 
    IsNaN   = Util.cell.isnan_cell(Val);
    
    BestVal = cell(Nh,1);
    for Ih=1:1:Nh
        % for each Header
        switch lower(Method)
            case 'first'
                % first not NaN value
                Ind = find(~IsNaN(Ih,:),1,'first');
                if (isempty(Ind))
                    BestVal{Ih} = NaN;
                else
                    BestVal{Ih} = Val{Ih,Ind};
                end
            case 'last'
                % first not NaN value
                Ind = find(~IsNaN(Ih,:),1,'last');
                if (isempty(Ind))
                    BestVal{Ih} = NaN;
                else
                    BestVal{Ih} = Val{Ih,Ind};
                end
            case 'max'
                % max value
                BestVal{Ih} = max([Val{Ih,:}],[],2);
            case 'min'
                % min value
                BestVal{Ih} = min([Val{Ih,:}],[],2);
            case 'mean'
                % mean value
                BestVal{Ih} = nanmean([Val{Ih,:}],2);
            otherwise
                error('Unknown Method option');
        end
    end

end

    
