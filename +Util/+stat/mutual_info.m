function MI=mutual_info(X,Y,varargin)
% Calculate the mutual information of two vectors (degree of independency)
% Package: Util.stat
% Description: Calculate the mutual information of two vectors.
% Input  : - A vector (X).
%          - A vector (Y).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'XLim' - 
%            'YLim' - 
%            'NbinX'-
%            'NbinY'-
% Output : - Mutual information
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: MI=Util.stat.mutual_info(randn(1000,1),rand(1000,1))
% Reliable: 
%--------------------------------------------------------------------------



DefV.XLim                = [];
DefV.YLim                = [];
DefV.NbinX               = 10;
DefV.NbinY               = 10; 
DefV.Kernel              = 'poisson';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if isempty(InPar.XLim)
    InPar.XLim = [min(X)-10.*eps, max(X)+10.*eps];
end

if isempty(InPar.YLim)
    InPar.YLim = [min(Y)-10.*eps, max(Y)+10.*eps];
end
    
Xrange = InPar.XLim(2) - InPar.XLim(1);
Yrange = InPar.YLim(2) - InPar.YLim(1);

Xedges = (InPar.XLim(1):Xrange./InPar.NbinX:InPar.XLim(2));
Yedges = (InPar.YLim(1):Yrange./InPar.NbinY:InPar.YLim(2));

Pxy = histcounts2(X,Y,Xedges,Yedges);
Pxy = Pxy./sum(Pxy(:));

% Kernel
switch lower(InPar.Kernel)
    case 'none'
        % do nothing
    case 'poisson'
        % Poisson matched filter (Ofek & Zackay 2018)
        [Ny,Nx] = size(Pxy);
        GK = Kernel2.gauss(2,2,0,Nx,Ny); % Gaussian Kernel
        PK = log(1+GK); % Poisson noise Gaussian Kernel
        Pxy = ifft2(fft2(Pxy).*(fft2(PK)));
        
    otherwise
        error('Unknown Kernel option');
end



Px  = sum(Pxy,1);
Py  = sum(Pxy,2);

LgP = log2(Pxy./(Px.*Py));
FlagN0 = Pxy>0;

MI = Pxy(FlagN0).*LgP(FlagN0);
MI = sum(MI(:));

