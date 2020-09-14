function [H,Hp]=loglogneg(X,Y,varargin)
% A loglog plot in which negative values are ploted at -log(abs(val))
% Package: plot
% Description: A log log plot that allows plotting negative numbers.
%              For example 1,-1 will be plotted at log10(1), -log10(abs(-1))
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [H]=plot.loglogneg(AccR,AccD,'.','MarkerSize',1);
% Reliable: 
%--------------------------------------------------------------------------




Fpp = X>0 & Y>0;
Fpn = X>0 & Y<0;
Fnp = X<0 & Y>0;
Fnn = X<0 & Y<0;


LXpp = log10(X(Fpp));
LYpp = log10(Y(Fpp));
LXpn = log10(X(Fpn));
LYpn = log10(abs(Y(Fpn)));
LXnp = log10(abs(X(Fnp)));
LYnp = log10(abs(Y(Fnp)));
LXnn = log10(abs(X(Fnn)));
LYnn = log10(abs(Y(Fnn)));


IsH = ishold;

H=plot.subplot1(2,2,'Gap',[0.0 0.0]);

% top-left
plot.subplot1([1 1])
Hp(1) = plot(LXnp,LYnp,varargin{:});
H(1).XDir = 'reverse';
H(1).YDir = 'normal';
axis([-2 2 -2 2])

% top -right
plot.subplot1([1 2])
Hp(2) = plot(LXpp,LYpp,varargin{:});
H(3).XDir = 'normal';
H(3).YDir = 'normal';
axis([-2 2 -2 2])

% bottom-left
plot.subplot1([2 1])
Hp(3) = plot(LXnn,LYnn,varargin{:});
H(2).XDir = 'reverse';
H(2).YDir = 'reverse';
axis([-2 2 -2 2])

% bottom-right
plot.subplot1([2 2])
Hp(4) = plot(LXpn,LYpn,varargin{:});
H(4).XDir = 'normal';
H(4).YDir = 'reverse';
axis([-2 2 -2 2])


if (~IsH)
    hold off;
end

