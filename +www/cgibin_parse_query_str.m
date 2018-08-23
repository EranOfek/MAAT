function [ST,CT]=cgibin_parse_query_str(QS,Decode,ConvertToNum)
% Break a URL parameters query string to parameter names and values.
% Package: www
% Description: Break a URL parameters query string to parameter names
%              and values.
% Input  : - URL query string. If empty, then get the query string from
%            the QUERY_STRING environment variable. Default is empty.
%          - decode strings {false|true}. Default is true.
%          - A boolean scalar or vector indicating if the value strings
%            in the structure
%            should be transformed to numerical values. Default is false;
% Output : - Structure in which the field names represent parameter
%            names, and the field value (content) represent the parameter
%            value.
%          - Cell array in which each cell contains a cell array of the
%            {Par,Value} pair.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Dec 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% [ST,CT]=www.cgibin_parse_query_str('q=parse+post+get+cgi-bin+matlav&oq=parse+post+get+cgi-bin+matlav');
% Reliable: 2
%--------------------------------------------------------------------------

Def.QS     = [];
Def.Decode = true;
Def.ConvertToNum = false;
if (nargin==0),
    QS       = Def.QS;
    Decode   = Def.Decode;
    ConvertToNum = Def.ConvertToNum;
elseif (nargin==1),
    Decode   = Def.Decode;
    ConvertToNum = Def.ConvertToNum;
elseif (nargin==2),
    ConvertToNum = Def.ConvertToNum;
elseif (nargin==3),
    % do nothing
else
    error('Illegal number of input arguments');
end

if (isempty(QS)),
    QS = getenv('QUERY_STRING');
end

%'q=parse+post+get+cgi-bin+matlav&oq=parse+post+get+cgi-bin+matlav';
%S=regexp(['&',QS,'&'],'\&(?<Par>.+)\=(?<Value>.+)\&','names')

CT=regexp(regexp(QS,'&','split'),'=','split');
Ns = length(CT);

if (length(ConvertToNum)==1),
    if (ConvertToNum),
        ConvertToNum = true(Ns,1);
    else
        ConvertToNum = false(Ns,1);
    end
end

if (~isempty(CT{1}{1})),
    
    for Is=1:1:Ns,
        if (Decode),
            ST.(CT{Is}{1}) = urldecode(CT{Is}{2});
        else
            ST.(CT{Is}{1}) = CT{Is}{2};
        end
        if (ConvertToNum(Is)),
            ST.(CT{Is}{1}) = str2double(ST.(CT{Is}{1}));
        end
    end
else
    ST = [];
end