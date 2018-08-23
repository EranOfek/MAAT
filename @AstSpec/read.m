function AS=read(AS,File,varargin)
% Read spectra using spec_read.m and add it to an AstSpec object.
% Package: @AstSpec
% Description: Read spectra using spec_read.m and add it to an AstSpec
%              object.
% Input  : - AstSpec class object.
%          - Files from which to read spectra and concat to
%            input AstSpec.
%            See spec_read.m for options.
%          * Additional arguments to pass to spec_read.m
% Output : - AstSpec class object.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AS=read(AS,'file.fits');
% Reliable: 2
%--------------------------------------------------------------------------

AS = [AS; spec_read(File,varargin{:})];

