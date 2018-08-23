function [Flat,Lambda]=flat_optimal(Cube,varargin)
% SHORT DESCRIPTION HERE
% Package: ImUtil.Im
% Description: 
% Input  : - Cube of Nim X Ny X Nx
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Nov 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------


if (nargin==0)
    % simulation mode
    Level  = 20000;
    LambdaSim = rand(10,1).*0.3 + 1;
    Nim    = 10;
    Ny     = 100;
    Nx     = 100;
    Cube   = zeros(Nim,Ny,Nx);
    for Iim=1:1:Nim
        Cube(Iim,:,:) = poissrnd(round(ones(Ny,Nx).*Level.*LambdaSim(Iim)));
    end
end

%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

for Iim=1:1:Nim
    Tmp = Cube(Iim,:,:);
    Lambda(Iim) = mean(Tmp(:));
end
Lambda = Lambda./mean(Lambda);


Size = size(Cube);
Nim  = Size(1);

%Lambda = ones(Nim,1);
%Flat   = ones(Size(2:3));

Conv = false;
while (~Conv)

    Flat   = sum(Lambda.*Cube,1)./Nim;

    Lambda = squeeze(sum(sum(Flat.*Cube,2),3))./squeeze(sum(sum(Cube.^2,2),3));
    
    Lambda./LambdaSim
    
end

