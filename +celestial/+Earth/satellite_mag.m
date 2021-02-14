function [Res]=satellite_mag(varargin)
% Satellite apparent magnitude
% Package: celestial
% Description: Satellite apparent magnitude
% Input  : * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Dist' - Satellite distance [km]. Default is 500.
%            'Area' - Satellite area [cm^2]. Default is 1.
%            'Albedo' - Albedo. Default is 0.1.
%            'SunAbsMag' - Sun abs. mag. Default is -26.74 (V Vega).
% Output : - Structure containing the following fields:
%            'Mag' - Satellite apparent magnitude.
%            'AngV' - Satellite angular velicity assuming Keplerian orbit
%                     with no projection [arcsec/s].
%            'MagResEl' - Magnitude per resolution element (FWHM-length).
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Dec 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=celestial.Earth.satellite_mag
% Reliable: 2
%--------------------------------------------------------------------------

RAD = 180./pi;

DefV.Dist                 = 500;     % km
DefV.Area                 = 15;       % [cm^2]
DefV.Albedo               = 0.1;
DefV.FWHM                 = 2.5;
DefV.SunAbsMag            = -26.74;  % V
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Res.Mag = InPar.SunAbsMag - 2.5.*log10(InPar.Area ./ (4.*pi.*(InPar.Dist.* 1e5)^2)) - 2.5.*log10(InPar.Albedo);

K = celestial.Kepler.kepler3law(constant.EarthM,'a',constant.EarthR+InPar.Dist.*1e5);
Res.AngV = K.v./(InPar.Dist.*1e5).*RAD.*3600;  % ''/s

Res.MagResEl = Res.Mag + 2.5.*log10(Res.AngV./InPar.FWHM);


%%

if 1==0
    ExpTime = 15;
    PixScale = 1.2;
    FWHM     = 2.5;

    Ntelescope = 48;
    % plots
    DistV=logspace(log10(300),log10(40000),10)';
    %DistV=logspace(log10(300),log10(1e7),100)';
    Nd   = numel(DistV);
    for Id=1:1:Nd
        Res=celestial.Earth.satellite_mag('Area',1,'Dist',DistV(Id));
        %Res.Mag
        LengthPix  = min(6000,Res.AngV.*ExpTime./PixScale);
        LengthFWHM = min(6000./FWHM,Res.AngV.*ExpTime./FWHM);
        
        
%         [SN,Op]=telescope.sn.sn_calc('SN',7,'Aper',27,'FocalLength',62,'FWHM',sqrt(FWHM.*LengthFWHM),'ExpTime',ExpTime,'PixSize',3.67,'RN',1.5);
%         AreaStreak(Id) = 10.^(0.4.*(Res.MagResEl - SN.Mag));
%         
        ExpTimeStationary = FWHM./Res.AngV;
        [SN,Op]=telescope.sn.sn_calc('SN',5,'Aper',27,'FocalLength',62,'FWHM',FWHM,'ExpTime',ExpTimeStationary,'PixSize',3.67,'RN',1.5);
        AreaPoint(Id) = 10.^(0.4.*(Res.Mag - SN.Mag));
        
        
        MaxResEl = min(Res.AngV.*ExpTime,10000)./FWHM;
        Npt = min(ExpTime./ExpTimeStationary, MaxResEl);
        AreaStreak(Id) = AreaPoint(Id)./sqrt(Npt);
        
        
    end
    
    loglog(DistV,sqrt(AreaStreak),'k-','LineWidth',2)    
    hold on
    loglog(DistV,sqrt(AreaStreak./Ntelescope),'k--','LineWidth',2)    
    
    loglog(DistV,sqrt(AreaPoint),'-','Color',[0.8 0.8 0.8],'LineWidth',2) 
    loglog(DistV,sqrt(AreaPoint./Ntelescope),'--','Color',[0.8 0.8 0.8],'LineWidth',2) 
    legend('Streak detection','Sreak x 48 tel','Point detection','Point x 48 tel','Location','SouthEast');
    axis([300 40000 0.3 2e2]);
    H=xlabel('Distance [km]');
    H.FontSize=18;
    H.Interpreter = 'latex';
    H=ylabel('Size [cm]');
    H.FontSize=18;
    H.Interpreter = 'latex';
    %print SatSizeDist.eps -depsc2
    
end
    

if 1==0
    ExpTime = 15;
    PixScale = 1.2;
    FWHM     = 2.5;

    Ntelescope = 48;
    
    VecDist = logspace(3,7,40);
    VecVel  = logspace(log10(1),log10(70),10);
    Ndist   = numel(VecDist);
    Nvel    = numel(VecVel);
    
    AngVel  = VecVel./VecDist.' .*RAD.*3600;  % "/s
    AreaStreak = nan(Ndist,Nvel);
    
    for Idist=1:1:Ndist
        for Ivel=1:1:Nvel
            
            Res=celestial.Earth.satellite_mag('Area',1,'Dist',VecDist(Idist));
            %Res.Mag
            LengthPix  = min(1024,AngVel(Idist,Ivel).*ExpTime./PixScale);
            LengthFWHM = min(1024./FWHM,AngVel(Idist,Ivel).*ExpTime./FWHM);
            FracOfLight = 1024./(AngVel(Idist,Ivel).*ExpTime./PixScale);
            if FracOfLight>1
                FracOfLight = 1;
            end
            
            if LengthFWHM<1
                LengthFWHM = 1;
            end
            
            [SN,Op]=telescope.sn.sn_calc('SN',5,'Aper',27,'FocalLength',62,'FWHM',FWHM,'ExpTime',ExpTime,'PixSize',3.67,'RN',1.5.*sqrt(2.*LengthPix));
            MagLimit = SN.Mag-2.5.*log10(sqrt(LengthFWHM)) + 2.5.*log10(FracOfLight);
            Radius(Idist,Ivel) = sqrt((10.^(0.4.*(Res.Mag - MagLimit)))./pi);
        end
    end
    
    [C,H]=contour(VecDist,VecVel,(Radius'./100),[0.03 0.1 0.3 1 3 10 30 100])
    clabel(C,H,'LabelSpacing',100,'Color','k','FontWeight','bold')
    set(gca,'XS','log','YS','linear')     
    %colorbar
    axis([1e3 1e7 1 70]);
    
    % plot GeoStat ang vel
    GeoSatAngVel = 15; % arcsec/sec
    LowSatAngVel = 7.7./500.*RAD.*3600;
    VelGeo = VecDist.*GeoSatAngVel./(RAD.*3600);
    VelLow = VecDist.*LowSatAngVel./(RAD.*3600);
    hold on;
    plot(VecDist,VelGeo,'k-','LineWidth',2)
    plot(VecDist,VelLow,'k-','LineWidth',2)
    H = text(4.4e5,40,'Geo-sat ang. vel.')
    H.Rotation = 80;
    H = text(3.2e3,40,'Low-sat ang. vel.')
    H.Rotation = 78;
    % plot Keplerian orbit around the Earth
    Res=celestial.Kepler.kepler3law(5.9e27,'a',(6371+VecDist).*1e5);
    plot(VecDist,Res.v./1e5,'r--','LineWidth',2)
    
    
    H=xlabel('Distance [km]');
    H.FontSize = 18;
    H.Interpreter = 'latex';
    H=ylabel('Projected Velocity [km/s]');
    H.FontSize = 18;
    H.Interpreter = 'latex';
    
    % print DetectabilitySmallNEO.eps -depsc2

end

            
    
    
    

