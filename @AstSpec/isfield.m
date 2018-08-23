function obj=isfield(Sim,Field)
% Check if a string is a field in an AstSpec object.
% Package: @AstSpec
% Description: Check if a string is a field in an AstSpec object.
% Input  : - An AstSpec object.
%          - A field name.
% Output : - True if a field name, false otherwise.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: isfield(AstS,'Wave');
% Reliable: 2
%--------------------------------------------------------------------------

obj = any(strcmp(fieldnames(Sim),Field));