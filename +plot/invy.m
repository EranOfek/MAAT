function invy
% Invert the y-axis of the current axis.
% Package: plot
% Description: Invert the y-axis of the current axis.
%              This is a shortcut command to axis('ij') and
%              set(gca,'YDir','Inverse').
% Input  : null
% Output : null
% Tested : Matlab 4.2
%     By : Eran O. Ofek                    Feb 1994
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 1
%------------------------------------------------------------------------------
axis('ij');
