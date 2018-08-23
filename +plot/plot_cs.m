function H=plot_cs(X,Y,varargin)
%--------------------------------------------------------------------------
% plot_cs function                                                plotting
% Description: Plot colored symbols.
% Input  : - Vector of X.
%          - Vector of Y.
%          * Arbitrary number of pairs of input parameters: key,val,...
%            Available keywords are:
%            'Color' - Symbols edge color. This may be a 3-element
%                      raw vector of a single color to apply to all
%                      the data points; A 3-column matrix in which each
%                      raw is the color of one data point; or a column
%                      vector of numbers in the range 0 to 1 that will
%                      be map to colors using the provided ColorMap.
%                      Default is [0 0 1].
%            'MarkerFaceColor' - Symbols face colors. Like Color, but for
%                      the symbol face. Default is [0 0 1].
%            'ColorMap' - 3 column matrix of color map. Default is the
%                      default color map.
%            'Marker' - Marker type, or a cell array of marker types
%                      per data point. Default is 'o'. Alternatively, this
%                      can be an index for a marker type indicated by
%                      'MarkerMap'.
%            'MarkerMap' - If 'Marker' is numeric than it will be
%                      taken from MarkerMap. Default is
%                      {'.'; 'o'; 'x'; '+'; '*'; 's'; 'd'; 'v'; '^'; '<';
%                      '>'; 'p'; 'h'}.
%                      For Example Marker=2 will set the Marker to 'o'.
%            'MarkerSize' - Marker size, or a vector of marker sizes per 
%                      element. Default is 8.
%            'InterpMethod' - Color map interpolation method.
%                      Default is 'linear'.
%            'Hold' - hold on current figure. Default is false. This will
%                     not override hold on made prior to the function call.
%            'HoldEnd' - hold on the final figure. Default is false.
% Output : - Vecor of handels for each symbol.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    May 2007
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: H=plot_cs(rand(100,1),rand(100,1),'Color',rand(100,3),'MarkerFaceColor',rand(100,3),'MarkerSize',rand(100,1).*10)
% Reliable: 2
%--------------------------------------------------------------------------

DefV.Color           = [0 0 1];
DefV.MarkerFaceColor = [0 0 1];
DefV.Marker          = 'o';
DefV.MarkerMap       = {'.'; 'o'; 'x'; '+'; '*'; 's'; 'd'; 'v'; '^'; '<'; '>'; 'p'; 'h'};
DefV.MarkerSize      = 8;
DefV.ColorMap        = colormap;
DefV.InterpMethod    = 'linear';
DefV.Hold            = false;     % hold on current figure
DefV.HoldEnd         = false;     % hold on at the end
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

if (size(InPar.Color,2)==1),
    Ncm = size(ColorMap,1);
    InPar.Color           = interp1((0:1./(Ncm-1):1).',InPar.ColorMap,InPar.Color,InPar.InterpMethod);
    InPar.MarkerFaceColor = interp1((0:1./(Ncm-1):1).',InPar.ColorMap,InPar.MarkerFaceColor,InPar.InterpMethod);
end

if (isnumeric(InPar.Marker)),
    InPar.Marker = InPar.MarkerMap(InPar.Marker);
end
    
if (~iscell(InPar.Marker)),
   InPar.Marker = {InPar.Marker};
end

if (InPar.Hold)
    hold on;
end


N = numel(X);
H = zeros(N,1);
for I=1:1:N,
    
   H(I) = plot(X(I),Y(I),'.');
   hold on;
   set(H(I),'Marker',InPar.Marker{min(I,numel(InPar.Marker))},...
            'MarkerFaceColor',InPar.MarkerFaceColor(min(I,size(InPar.MarkerFaceColor,1)),:),...
            'MarkerSize',InPar.MarkerSize(min(I,numel(InPar.MarkerSize))));
end

if (InPar.HoldEnd)
    hold on;
else
    hold off;
end