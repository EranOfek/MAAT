function spec_select_allaper(Image,varargin)
%--------------------------------------------------------------------------
% spec_select_allaper function                                      ImSpec
% Description: 
% Input  : -
% Output : - 
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Feb 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


ImageField  = 'Im';
HeaderField = 'Header';
FileField   = 'ImageFileName';
MaskField   = 'Mask';
BackImField = 'BackIm';
ErrImField  = 'ErrIm';


DefV.DispDir       = 'x';         % Dispersion direction
DefV.Range         = [];          % range to collapse
DefV.Ncollapse     = 3;           % number of collpases in range 
DefV.AdjPlotBack   = 'div';       % adjust plot background {'div'|'sub'|'no'}
% find_contrast_peaks.m 
DefV.PeakAlgo      = 'RMSmin';
DefV.ContrastRMS   = 8;
% spec_collapse_dispaxis.m parameters
DefV.CollapseAlgo  = 'optimal'; % {'optimal'|'mean'|'median','quantile'}
DefV.RN            = 7;   % [e-]
DefV.Prct          = 0.9; 
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

% deal with Image input
Sim = image2sim(Image);

% rotate image so dispersion direction is x-axis
switch lower(InPar.DispDir)
    case 'y'
        Sim.(ImageField) = [Sim.(ImageField)].';
    otherwise
        % do nothing
end

% Define collapse range
Size = size(Sim.(ImageField));
if (isempty(InPar.Range)),
    InPar.Range = [1 Size(2)];
end
    
% Collapse of entire image range
Collapse = spec_collapse_dispaxis(Sim.(ImageField),varargin{:},'DispDir','x','Range',InPar.Range);
VecSpat  = (1:1:Size(1)).';
VecDisp  = (1:1:Size(2)).';

% Initial search for sources
[SelectedMaxima,Contrast_UnitsRMS]=find_contrast_peaks(VecSpat,Collpase{1},InPar.ContrastRMS,varargin{:});

got here


    ColSize = range(InPar.Range)./InPar.Ncollapse;
    InPar.Range = round([[InPar.Range(1):ColSize:(InPar.Range(2)-ColSize)].',[InPar.Range(1)+ColSize:ColSize:InPar.Range(2)].']);
    
    Colors = generate_colors(InPar.Ncollapse);
    String = cell(1,InPar.Ncollapse);
    for Ic=1:1:InPar.Ncollapse,
       Collapse = spec_collapse_dispaxis(Sim.(ImageField),varargin{:},'DispDir','x','Range',InPar.Range(Ic,:));
       switch lower(InPar.AdjPlotBack)
           case 'div'
               PlotVector = Collapse{1}./nanmedian(Collapse{1});
           case 'sub'
               PlotVector = Collapse{1} - nanmedian(Collapse{1});
           case {'no','none'}
               PlotVector = Collapse{1};
           otherwise
               error('Unknown AdjPlotBack option');
       end
       H=plot(VecSpatI,PlotVector); %, '-','Color',Colors(Ic,:))
       set(H,'Color',Colors(Ic,:));
       hold on;
       String{Ic} = sprintf('Disp %d..%d',InPar.Range(Ic,:));
    end
    legend(String{:});
    H = xlabel('Spatial position [pix]');
    set(H,'FontSize',18);
    H = ylabel('Normalized counts');
    set(H,'FontSize',18);
