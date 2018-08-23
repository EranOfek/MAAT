function DirCon=cdsclient_prog_names
% Return the list of programs in the cdsclient directory
% Package: VO.VizieR
% Input  : null
% Output : - Strcture array of cdsclient directory content.
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: DirCon=VO.VizieR.cdsclient_prog_names; {DirCon.name}
% Reliable: 2

Path = VO.VizieR.cdsclient_path;

DirCon = dir(Path);