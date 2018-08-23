function plot_invchildren(H)
%--------------------------------------------------------------------------
% plot_invchildren function                                       plotting
% Description: Invert the order of the childrens under a given handle.
% Input  : - Handle. Default is gca.
% Output : null
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Sep 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: plot_invchildren;
% Reliable: 2
%--------------------------------------------------------------------------

Def.H = [];
if (nargin==0),
    H = Def.H;
end
if (isempty(H)),
    H = gca;
end

Hc = get(H,'Children');
set(H,'Children',flipud(Hc));


    