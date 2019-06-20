function [AstC,Col]=aper_phot(Sim,Coo,varargin)
% Aperture photometry for defined positions in SIM images.
% Package: @SIM
% Description: Perform aperture photometry on a SIM object.
% Input  : - A SIM object.
%          - Coordinates around which to perform aperture photometry.
%            Either a two column matrix [X,Y], or an AstCat object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            
% Output : - 
% See also: aperphot.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 2
%--------------------------------------------------------------------------

CatField    = AstCat.CatField;

DefV.AperRad            = [2 4 8 12 16];
DefV.Annulus            = [16 22];
DefV.ColXY              = {'XWIN_IMAGE','YWIN_IMAGE'};
DefV.BackFun            = @median; %@rmean;           %@median;
DefV.BackFunPar         = {}; %{1,[0.1 0.1]};    % {};
DefV.BackErrFun         = @std;
DefV.BackErrFunPar      = {};
DefV.Back               = [];
DefV.Field              = SIM.ImageField;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Naper    = numel(InPar.AperRad);
if (AstCat.isastcat(Coo))
    XY = col_get(Coo,InPar.ColXY);
else
    XY = Coo;
end

Col.X       = 1;
Col.Y       = 2;
Col.Back    = 3;
Col.BackErr = 4;
Col.Aper    = (5:1:4+Naper);
Col.AperErr = (max(Col.Aper)+1:1:(max(Col.Aper)+Naper));


Nsim = numel(Sim);
AstC = AstCat(size(Sim));
for Isim=1:1:Nsim
    % for each SIM element
    %Image = Sim(Isim).(InPar.Field);
    Image = double(Sim(Isim).(InPar.Field)); % Na'ama, 20180828
    
    SizeIm = size(Image);
    %[MatX,MatY] = meshgrid((1:1:SizeIm(2)),(1:1:SizeIm(1)));

    MaxR = max(max(InPar.AperRad),max(InPar.Annulus));

    AperRad2 = InPar.AperRad.^2;
    Annulus2 = InPar.Annulus.^2;

    [CI,CR2] = ImUtil.Im.find_within_radius_cell(SizeIm,XY(:,1),XY(:,2),MaxR,false);

    Nback = numel(InPar.Back);

    Nsrc = size(XY,1);
    AstC(Isim).(CatField)  = zeros(Nsrc,max(Col.Aper));
    AstC(Isim).(CatField)(:,[Col.X, Col.Y]) = XY;
    for Isrc=1:1:Nsrc
        %sum(Image(CI{Isrc}));

        % background per pixels
        IndBack  = CI{Isrc}(CR2{Isrc}>Annulus2(1) & CR2{Isrc}<Annulus2(2));
        BackArea = numel(IndBack);
        AstC(Isim).(CatField)(Isrc,Col.Back) = InPar.BackFun(Image(IndBack),InPar.BackFunPar{:}); %./BackArea;  % already per pixel

        % background noise
        % Bug fix - sometimes background is empty/NaN
        Tmp = InPar.BackErrFun(Image(IndBack),InPar.BackErrFunPar{:});
        if isempty(Tmp)
            AstC(Isim).(CatField)(Isrc,Col.BackErr) = NaN;
        else
            AstC(Isim).(CatField)(Isrc,Col.BackErr) = Tmp;
        end

        if (isempty(InPar.Back))
            % User didnt supply background - use background as estimated from
            % annulus
            SrcBack = AstC(Isim).(CatField)(Isrc,Col.Back);
        else
            % user supplied background - user user background
            % however, The output catalog will contain the background measured
            % in the annulus
            SrcBack = InPar.Back(min(Isrc,Nback));
        end


        % aperture photometry
        for Iaper=1:1:Naper
            IndAper = CI{Isrc}(CR2{Isrc}<AperRad2(Iaper));
            AperArea = numel(IndAper);
            SubImage = Image(IndAper);
            %AperArea = sum(Image(R2(:)<InPar.AperRad(Iaper).^2));
            %Cat(Isrc,Col.Aper(Iaper)) = sum(Image(R2(:)<InPar.AperRad(Iaper).^2)) - AperArea.*Cat(Isrc,Col.Back);
            AstC(Isim).(CatField)(Isrc,Col.Aper(Iaper))   = sum(SubImage) - AperArea.*SrcBack;
            AstC(Isim).(CatField)(Isrc,Col.AperErr(Iaper)) = sqrt(sum(SubImage));
        end

    end
end