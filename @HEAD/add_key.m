function Head=add_key(Head,varargin)
% Add new keyword, value and comment lines to an Header object.
% Package: @HEAD
% Description: Add new keyword, value and comment lines to an Header
%              object without checking if keyword exist.
%              See replace_key.m for keyword replacment/addition.
% Input  : - An Header object.
%          * Either an arbitrary number of triplets:
%            ...,key,val,comment,...;
%            or an arbitrary number of cell arrays, each cell array
%            containing 2 or 3 columns of {key,val,comment}, where the
%            comment field is optional;
%            or a single element HEAD object.
% Output : - The Header object with the added lines.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: H.add_key('A','2','aa','B',3,'com');
%          H.add_key({'A','2';'aa','B',3});
%          H = add_key(S,S(1));
% Reliable: 2
%--------------------------------------------------------------------------

HeaderField = 'Header';


if (~iscell(varargin{1}) && ~HEAD.ishead(varargin{1}))
    % assume varargin contains triplets of key,val,comments
    Ntriplets   = numel(varargin)./3;
    Tmp = reshape(varargin,3,Ntriplets).';
    clear varargin;
    varargin{1} = Tmp;
end

Narg = numel(varargin);

Nh = numel(Head);
for Ih=1:1:Nh
    for Iarg=1:1:Narg
        if (HEAD.ishead(varargin{Iarg}))
            varargin{Iarg} = varargin{Iarg}.(HeaderField);
        end
        if (size(varargin{Iarg},2)==2)
            % add comment column
            Nrow  = size(varargin{Iarg},1);
            Tmp = cell(Nrow,1);
            [Tmp{1:end}]=deal('');
            varargin{Iarg} = [varargin{Iarg}, Tmp];
        end
        
        Head(Ih).(HeaderField) = [Head(Ih).(HeaderField); varargin{Iarg}];
    end
end
    
    
