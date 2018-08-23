function VArg=struct2varargin(Struct)
% Structure field name and content to cell array of key,val pairs
% Package: Util.struct
% Description: Given a structure, prepare a cell array of all the
%              field_names, field_values,...
%              This is useful for converting InPar to varargin input.
% Input  : - A structure.
% Output : - A cell array of all the field_names, field_values,...
% Tested : Matlab R2013a
%     By : Eran O. Ofek                    Aug 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: VArg=Util.struct.struct2varargin(InPar);
% Reliable: 2
%--------------------------------------------------------------------------

FN = fieldnames(Struct);
FV = struct2cell(Struct);

VArg = cell(numel(FN).*2,1);
VArg(1:2:end-1) = FN;
VArg(2:2:end)   = FV;


