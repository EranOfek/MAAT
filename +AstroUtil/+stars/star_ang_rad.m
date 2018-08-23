function [Theta,FlagInRange]=star_ang_rad(Mag,Color,MagBand,ColorBand)
% Calculate empirical angular radii of stars from magnitude and colors.
% Package: AstroUtil.stars
% Description: Empirical angular radii of stars based on their magnitude
%              and colors.
% Input  : - Magnitude.
%          - Color.
%          - Band in which magnitude is provided. Default is 'g'.
%          - Color name in which color is provided. Default is 'g-r'.
% Output : - Angular radii of stars [mas].
%          - Flag indicating if angular radii is in range (1) or out of
%            range (0) and therefore unreliable.
% Reference: http://arxiv.org/abs/1311.4901
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Jan 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Theta,FlagInRange]=AstroUtil.stars.star_ang_rad(10,0.4,'g','g-i');
% Reliable: 2
%--------------------------------------------------------------------------


% http://arxiv.org/abs/1311.4901
Table = {'V','B-V',[-0.02, 1.73], [0.49612 0.00111],[1.11136 0.00939],[-1.18694 0.02541],[ 0.91974 0.02412],[-0.19526 0.00738], 7.8;...
         'g','g-r',[-0.23, 1.40], [0.66728 0.00203],[0.58135 0.01180],[ 0.88293 0.03470],[-1.41005 0.04331],[0.67248 0.01736],  9.7;...
         'g','g-i',[-0.43, 2.78], [0.69174 0.00125],[0.54346 0.00266],[-0.02149 0.00097],[0 0],            [0 0],              9.2};
     


Imag   = find(strcmp(MagBand,Table(:,1)) & strcmp(ColorBand,Table(:,2)));

if (~isempty(Imag)),
    Data = [Table{Imag,3:end}];
    % check if color in range
    FlagInRange = Color>Data(:,1) & Color<Data(:,2);
    Ai   = Data(3:2:end-1);
    Veci = (0:1:length(Ai)-1);
    %Theta0 = 10.^(sum(bsxfun(@times, Ai, bsxfun(@power,Color,Veci)),2))
    
    Theta0 = 10.^(sum(bsxfun(@times, Ai, bsxfun(@power,Color, Veci)),2));   % ang diamater at mag=0
    Theta  = 10.^(-0.2.*Mag).*Theta0 .*0.5;   % [mas] angular radius
else
    Theta = NaN;
end
    
    
    