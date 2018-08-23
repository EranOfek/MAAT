function Path=cdsclient_path(varargin)
% Return the path of the local cdsclient directory
% Package: VO.VizieR
% Description: Return the path of the local cdsclient directory.
%              Edit in order to change the parameters.
% Input  : * Optional keyword 'Path' and its value.
%            Default is '~/matlab/bin/cdsclient-3.80/'.
% Output : - Local cdsclient directory path.
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Path=VO.VizieR.cdsclient_path('Path','~/matlab/bin/cdsclient-3.80/');
% Reliable: 2

DefV.Path                 = '~/matlab/bin/cdsclient-3.80/';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);
Path = InPar.Path;
