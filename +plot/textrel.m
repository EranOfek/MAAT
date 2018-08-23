function H=textrel(Xrel,Yrel,Text,varargin);
%--------------------------------------------------------------------------
% textrel function                                                plotting
% Description: Write a text in current axis, where the position of the
%              text is specified as a relative position in respect to
%              the axis limits.
% Input  : - Relative X position.
%          - Relative Y position.
%          - Text
%          * Arbitrary number of input arguments to pass to the text.m
%            command.
% Output : - Text handle.
% Tested : Matlab 2011b
%     By : Eran O. Ofek                     May 2013
%    URL : http://weizamann.ac.il/home/eofek/matlab/
% Example: H=textrel(0.1,0.9,'Hello');
% Reliable: 2
%--------------------------------------------------------------------------

XLim = get(gca,'XLim');
YLim = get(gca,'YLim');
XScale = get(gca,'XScale');
YScale = get(gca,'YScale');

switch lower(XScale)
    case 'linear'
        X = XLim(1) + (XLim(2)-XLim(1)).*Xrel;
    case 'log'
        X = XLim(1).*10.^(log10(XLim(2)./XLim(1)).*Xrel);
    otherwise
        error('Unknown XScale option');
end
        
switch lower(XScale)
    case 'linear'
        Y = YLim(1) + (YLim(2)-YLim(1)).*Yrel;
    case 'log'
        Y = YLim(1).*10.^(log10(YLim(2)./YLim(1)).*Yrel);
    otherwise
        error('Unknown XScale option');
end

H = text(X,Y,Text,varargin{:});

