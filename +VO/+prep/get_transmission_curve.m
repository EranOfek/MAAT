function AstF=get_transmission_curve(Telescope)
% Read astronomical filters from WWW into an AstFilter object
% Package: VO.prep
% Description: Read astronomical filters from the filter profile service
%              into an AstFilter object.
% Input  : - Telescope name (see http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=Swift
%            for options).
% Output : - An AstFilter object.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Nov 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstF=VO.prep.get_transmission_curve
% Reliable: 2
%--------------------------------------------------------------------------


%DefV. = 
%InPar = InArg.populate_keyval(DefV,varargin,mfilename);

%Telescope = 'HST'


URL1 = sprintf('http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?mode=browse&gname=%s',Telescope);

L = urlread(URL1);

I1=strfind(L,'filters');
I2=min(strfind(L,'Filter ID'));
L1 = L(I1:I2);
Tmp = regexp(L1,'gname2=(?<aa>\w+)"','names');
InstList = {Tmp.aa};
Ninst = numel(InstList);
AstF = AstFilter;
K = 0;
% all instruments
if (Ninst==0)
    Ninst=-1;
end
for Iinst=1:1:abs(Ninst)
    %Iinst, InstList{Iinst}
    if (Ninst==-1)
        URL11 = sprintf('http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=&&mode=browse&gname=%s',Telescope);
        Li    = urlread(URL11);
        Str = sprintf('%s/(?<ii>\\w+)\\.(?<aa>\\w+)</a></td><',Telescope);
        Tmp = regexp(Li,Str,'names');
        InstList{Iinst}=Tmp(1).ii;
        FiltList = {Tmp.aa};
        Nf = numel(FiltList);
    else
        URL11 = sprintf('http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=&&mode=browse&gname=%s&gname2=%s#filter',Telescope,InstList{Iinst});
        Li    = urlread(URL11);
        Str = sprintf('%s/%s\\.(?<aa>\\w+)</a></td><',Telescope,InstList{Iinst});
        Tmp = regexp(Li,Str,'names');
        FiltList = {Tmp.aa};
        Nf = numel(FiltList);
    end
    
    for If=1:1:Nf
        %If, FiltList{If}
        K = K+1;
        URL2 = sprintf('http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=%s/%s.%s&&mode=browse&gname=%s&gname2=%s#filter',...
                        Telescope, InstList{Iinst},FiltList{If},Telescope,FiltList{If});
        Lf = urlread(URL2);
        Tmp1 = regexp(Lf,'fps.php\?ID=(?<aa>.+)">VOTable</a>','names');

        URL3 = sprintf('http://svo2.cab.inta-csic.es/svo/theory/fps3/getdata.php?format=ascii&id=%s',Tmp1.aa);
        TC = urlread(URL3);
        Tmp4 = textscan(TC,'%f %f\n','CommentStyle','#');
        TC   = [Tmp4{1},Tmp4{2}];
        TC   = [TC(1,1)-0.01 0; TC; TC(end,1)+0.01 0];

        S=regexp(Tmp1.aa,'\.','split');
        Family = S{1};
        Band   = S{2};
        if (Ninst==-1)
           CC=regexp(Family,'/','split');
           if (strcmp(CC{1},CC{2}))
               Family = CC{1};
           end
        end
        AstF(K).family = Family;
        AstF(K).band   = Band;
        AstF(K).nT     = [TC(:,1),TC(:,2)./max(TC(:,2))];
        In0 = TC(:,2)>0;
        AstF(K).min_wl = min(TC(In0,1));
        AstF(K).max_wl = max(TC(In0,1));
        AstF(K).source = URL3;

     
        % calculate eff_wl
        AstF(K).eff_wl = nansum(AstF(K).nT(:,1).*AstF(K).nT(:,2))./nansum(AstF(K).nT(:,2));

        % calculate filter half_width
        % the width tranmit half the flux
       
        CumSum = cumsum(TC(:,2));
        CumSum = CumSum./CumSum(end);
        CumSum = CumSum + 10.*eps.*(1:1:numel(TC(:,1)))';
        AstF(K).half_width = interp1(CumSum,TC(:,1),0.75)-interp1(CumSum,TC(:,1),0.25);
        
        % calculate fwhm
        % the width at which the transmission drops by 1/2 relative
        % to max
        MaxnT = max(AstF(K).nT(:,2));
        I1 = find(AstF(K).nT(:,2)>(0.5.*MaxnT),1,'first');
        I2 = find(AstF(K).nT(:,2)>(0.5.*MaxnT),1,'last');
        AstF(K).fwhm = AstF(K).nT(I2,1) - AstF(K).nT(I1,1);

        [~,II] = max(TC(:,2));
        %AstF(K).peak = TC(II,1);
    end
end



%
if (1==0)
F = AstFilter.get;
unique({F.family})
F1=VO.prep.get_transmission_curve('HST');

F2=VO.prep.get_transmission_curve('Swift');
F3=VO.prep.get_transmission_curve('GAIA');
F4=VO.prep.get_transmission_curve('Kepler');
F5=VO.prep.get_transmission_curve('PAN-STARRS');
F6=VO.prep.get_transmission_curve('SkyMapper');
F7=VO.prep.get_transmission_curve('SOFIA');
F8=VO.prep.get_transmission_curve('Spitzer');
F9=VO.prep.get_transmission_curve('WISE');
F10=VO.prep.get_transmission_curve('LSST');
F11=VO.prep.get_transmission_curve('LasCumbres');

FF = [AstFilter.get('Johnson');...
     AstFilter.get('Cousins');...
     AstFilter.get('SDSS');...
     AstFilter.get('Johnson1975');...
     AstFilter.get('Bessel');...
     AstFilter.get('Stromgren');...
     AstFilter.get('2MASS');...
     AstFilter.get('DENIS');...
     AstFilter.get('DES');...
     AstFilter.get('GALEX');...
     AstFilter.get('IRAS');...
     AstFilter.get('KAIT');...
     AstFilter.get('LIM');...
     AstFilter.get('POSS/I');...
     AstFilter.get('POSS/II');...
     AstFilter.get('PTF');...
     AstFilter.get('ROSAT/WFC');...
     AstFilter.get('Ratir');...
     F1(:); F2(:); F3(:); F4(:); F5(:); F6(:); F7(:); F8(:); F9(:); F10(:); F11(:)];
 
end 
     