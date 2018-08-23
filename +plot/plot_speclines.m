function plot_speclines(SeriesName,RangeY,varargin)
%--------------------------------------------------------------------------
% plot_speclines function                                         plotting
% Description: Overplot series of spectral lines on top of a spectrum.
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jun 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

TextShift = 0.010;

DefV.Connect          = true;
DefV.Z                = 0;
DefV.LinesPar         = {'Color',[0.8 0.8 0.8],'LineStyle','-'};
DefV.TextPar          = {'FontSize',12};
DefV.Ser              = '';
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});


if (isnumeric(SeriesName)),
    Lines = SeriesName;
    Ser   = InPar.Ser;
else
    switch lower(SeriesName)
        case 'balmer5'
            Lines = [6564.6; 4862.7; 4341.7; 4102.9; 3971.2];
            Ser   = 'Balmer';
        case 'balmer4'
            Lines = [6564.6; 4862.7; 4341.7; 4102.9];
            Ser   = 'Balmer';
        case 'balmer3'
            Lines = [6564.6; 4862.7; 4341.7];
            Ser   = 'Balmer';
        case 'balmer2'
            Lines = [6564.6; 4862.7];
            Ser   = 'Balmer';
        case 'balmer1'
            Lines = [6564.6];
            Ser   = 'Balmer';
        case 'hei'
            Lines = [3888.6 4026.2 4471.5 5015.7 5875.6 6678.2 7065.2 7281.4]';
            Ser   = 'HeI';
        case 'heii'
            Lines = [3203.1 4685.7 5411.5 6560.1]';
            Ser   = 'HeII';
        
        otherwise
            error('Unknown SeriesName option');
    end
end

if (~isempty(InPar.Ser)),
    % Override Ser name
    Ser = InPar.Ser;
end

Lines = Lines.*(1+InPar.Z);

for Ilines=1:1:numel(Lines),
    Hp = plot([Lines(Ilines);Lines(Ilines)],[RangeY(1); RangeY(2)],InPar.LinesPar{:});
end

XLim = get(gca,'XLim');
if (~isnan(InPar.Connect)),
    MinL = min(Lines);
    MaxL = max(Lines);
    if (InPar.Connect>MaxL),
        plot([MinL;InPar.Connect],[RangeY(1); RangeY(1)],InPar.LinesPar{:});
        H = text(InPar.Connect.*(1+TextShift),RangeY(1),Ser);
        set(H,InPar.TextPar{:});
    elseif (InPar.Connect<MinL),
        plot([MaxL;InPar.Connect],[RangeY(1); RangeY(1)],InPar.LinesPar{:});
        H = text(InPar.Connect.*(1-TextShift),RangeY(1),Ser);
        set(H,'HorizontalAlignment','right');
    else
        error('Connect point must be smaller or larger than all lines wavelength');
    end
end

