function KVstr=construct_keyval_string(varargin)
%--------------------------------------------------------------------------
% construct_keyval_string function                                 General
% Description: Construct a ...,keyword, value,... string from pairs of
%              parameters.
% Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            The values will be transformed to strings and will be
%            appended after their keyword (e.g., 'c','4',...)
%            Alternatively, this can be a cell array of pairs, or a
%            single string.
%            If a single string then return the string as is.
% Output : - String of key,val pairs.
%            A "-" will be added before eack keyword string.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jan 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: KVstr=construct_keyval_string('c 4 m first');
%          KVstr=construct_keyval_string({'c',4,'m','first'});
%          KVstr=construct_keyval_string('c',4,'m','first');
% Reliable: 2
%--------------------------------------------------------------------------

Narg = numel(varargin);
if (Narg==1),
    if (ischar(varargin{1})),
        % do nothing
        KVstr = varargin{1};
    elseif (iscell(varargin{1})),
        KVstr = build_kv_string(varargin{1}{:});
    else
        error('Unknown input argument');
    end
else
    if (Narg.*0.5~=floor(Narg.*0.5)),
        error('Input must be even number of arguments');
    end
    KVstr = build_kv_string(varargin{:});
end


function KVstr=build_kv_string(varargin)

Narg = numel(varargin);
KVstr = '';
for Iarg=1:2:Narg-1,
    if (isnumeric(varargin{Iarg+1})),
        varargin{Iarg+1} = num2str(varargin{Iarg+1});
    end
    KVstr = sprintf('%s -%s %s',KVstr,varargin{Iarg},varargin{Iarg+1});
end


    