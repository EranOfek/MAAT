function L=contour_levels(C)
% Break a contour map into seperate contour levels.
% Package: plot
% Description: Given a contour map (the first output argument of contour)
%              return a structure array in which each element contains one
%              seperated contour.
% Input  : - Countour map.
% Output : - Structure array of leveles.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: L=plot.contour_levels(C)
% Reliable: 
%--------------------------------------------------------------------------



C = C.';
% collect all levels
N = size(C,1);
IndL = 0;
I    = 1;
L    = Util.struct.struct_def({'Level','Npt','Cont'},0,1);
if N>0
    Finish = false;
    while ~Finish
        IndL = IndL + 1;
        L(IndL).Level = C(I,1);
        L(IndL).Npt   = C(I,2);
        L(IndL).Cont  = C(I+1:I+L(IndL).Npt,:);
        
        I = I + L(IndL).Npt + 1;
        if (I>N)
            Finish = true;
        end
    end
end