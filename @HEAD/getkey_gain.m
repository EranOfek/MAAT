function Gain=getkey_gain(Head,varargin)
%--------------------------------------------------------------------------
% getkey_gain function                                         class/@HEAD
% Description: Search and get GAIN header keyword value from header.
% Input  : - An HEAD object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'GainKeys' - Cell array containin a list of header keywords in
%                         which to look for the GAIN value.
%                         Default is {'GAIN'}.
%            'SelectionMethod' - Method by which to select the value out
%                         of all possible keywords. Default is 'first'.
%                         This uses the getkey_fromlist.m function.
%            'OutType' - Gain output type:
%                        'mat' - Return a vector of gains, one per HEAD
%                                element. Default.
%                        'cell'- Retuen a cell array of gains.
% Output : - The gain values.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Gain=getkey_gain(Sim);
% Reliable: 2
%--------------------------------------------------------------------------


DefV.GainKeys           = {'GAIN'};
DefV.SelectionMethod    = 'first';
DefV.OutType            = 'mat';
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% get gain
% output is in a cell array
Gain = getkey_fromlist(Head,InPar.GainKeys,InPar.SelectionMethod,'Conv2Num',true);

% output type
switch lower(InPar.OutType)
    case 'mat'
        % convert the cell array into vector
        Gain = cell2mat(Gain);
    case 'cell'
        % do nothing
    otherwise
        error('Unknown OutType option');
end

