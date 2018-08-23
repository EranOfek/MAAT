function Answer=isfield_check(Struct,Field,Fun)
% If field exist run a function of fiel.
% Package: Util.struct
% Description: Check if a field name exist in structure, run a function
%              on the field content and return the function output.
%              
% Input  : - Structure.
%          - String containing field name.
%          - Function handle. Default is @all, so it will check if all
%            the field values are true.
% Output : - Return the function output. Retuern false, if the field
%            doesn't exist.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Feb 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Answer=Util.struct.isfield_check(Cat,'isHTM');
% Reliable: 
%--------------------------------------------------------------------------

if (nargin==2)
    Fun = @all;
end

Answer = false;
if (isfield(Struct,Field))
    Answer = Fun(Struct.(Field));
end
