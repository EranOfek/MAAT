function H=hist_stairs(X,N,Marker,varargin)
% Plot an histogram using stairs plot.
% Package: plot
% Description: Plot an histogram using stairs plot.
% Input  : - X.
%          - Y.
%          - Marker, default is 'k-'.
%          * Arbitrary number of pairs of ...,keyword,value,...
%            Avaliable keywords are:
%            'Type'  - {'v'|'h'} - for horizontal or vertical histogram.
%                      Default is 'v'.
% Output : - Handle to stairs plot.
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    Jul 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
%-----------------------------------------------------------------------------
Def.Marker = 'k-';
if (nargin==2)
   Marker = Def.Marker;
end

DefV.Type  = 'v';
%InPar = set_varargin_keyval(DefV,'y','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


SD = sortrows([X,N],1);
X  = SD(:,1);
N  = SD(:,2:end);

Nv   = length(X);
MidX = 0.5.*(X(2:end) + X(1:end-1));
XX   = zeros(Nv.*2,1);
Xm   = zeros((Nv-1).*2,1);
Xm(1:2:(Nv-1).*2-1) = MidX;
Xm(2:2:(Nv-1).*2) = MidX;
XX   = [X(1); Xm; X(end)];
NN   = zeros(size(XX));
NN(1:2:end-1) = N;
NN(2:2:end)   = N;

switch lower(InPar.Type)
 case 'v'
    H = plot(XX,NN,Marker);
 case 'h'
    H = plot(NN,XX,Marker);
 otherwise
    error('Unknown Type option');
end

