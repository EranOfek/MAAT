function [Data]=get_adata(RA,Dec,Radius,varargin);
%--------------------------------------------------------------------------
% get_adata function                                             Catalogue
% Description: Get 'everything' for a given celestial poistion.
% Input  : - J2000.0 RA in [rad], [H M S], or sexagesimal string.
%          - J2000.0 Dec in [rad], [Sign D M S], or sexagesimal string.
%          - Search radius [arcsec].
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            For most keywords the value can be 0 or >0, where 0 means
%            do not get data, and >0 means get data. For the catalogs,
%            >0 values used as the search radius in arcse.
%            The keywords can be on of the followings:
%            'GetIm'   - Copy FITS images from various sources {1|0},
%                        default is 0.
%            'GetCat'  - Search various catalogue {1|0}, default is 1.
%                        If 0, then return only links.
%            'Coo'     - Coordinates structure {1|0}, default is 1.
%            'SDSS'    - {1|0}, default is Radius.
%            'GALEX'   - {1|0}, default is Radius.
%            'POSS'    - {1|0}, default is 1.
%            'MASS'    - (2MASS) {1|0}, default is Radius.
%            'NED'     - {1|0}, default is Radius.
%            'FIRST'   - {1|0}, default is Radius.
%            'NVSS'    - {1|0}, default is Radius.
%            'MaxBCG'  - Look for MaxBCG clusters {1|0}, default is Radius.
%            'Abell'   - {1|0}, default is Radius.
%            'PGC'     - {1|0}, default is Radius.
%            'ROSAT'   - {1|0}, default is Radius.
%            'XMM'     - {1|0}, default is Radius.
%            'SNR'     - {1|0}, default is 1. - UNSEPORTED
%            'Pulsar'  - {1|0}, default is 1. - UNSEPORTED
%            'HST'     - {1|0}, default is 1. - UNSEPORTED
% Output : - Structure containing data output.
%            The fields are the keywords the user asked for.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                     April 2007
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: Data=get_adata([10 0 0],[1 41 0 0],60.*5);
% Reliable: 2
%--------------------------------------------------------------------------
RAD    = 180./pi;

DSS_FOV  = [6 6];

% convert to radians
RA     = convertdms(RA,'gH','r');
Dec    = convertdms(Dec,'gD','R');
Radius = Radius;

Get.GetIm    = 0;
Get.GetCat   = 1;
Get.SDSS     = Radius;
Get.GALEX    = Radius;
Get.POSS     = 1;
Get.MASS     = Radius;
Get.NED      = Radius;
Get.FIRST    = Radius;
Get.NVSS     = Radius;
Get.MaxBCG   = Radius;
Get.Abell    = Radius;
Get.PGC      = Radius;
Get.ROSAT    = Radius;
Get.XMM      = Radius;
Get.SNR      = 1;
Get.Pulsar   = 1;
Get.HST      = 1;


Narg = length(varargin);
for Iarg=1:2:Narg-1,
   if (isfield(Get,varargin{Iarg})==1),
      % field exist
      Get = setfield(Get,varargin{Iarg},varargin{Iarg+1});
   else
      error(sprintf('Unknown keyword: %s option',varargin{Iarg}));
   end
end

%-------------------
%--- Coordinates ---
%-------------------
Data.Coo.Eq       = [RA.*RAD, Dec.*RAD];
Data.Coo.EqSex{1} = convertdms(RA,'r','SH');
Data.Coo.EqSex{2} = convertdms(Dec,'r','SD');
GalCoo            = coco([RA, Dec],'j2000.0','g');
Data.Coo.Gal      = GalCoo.*RAD;
EclCoo            = coco([RA, Dec],'j2000.0','e');
Data.Coo.Ecl      = EclCoo.*RAD;



%------------
%--- SDSS ---
%------------
if (Get.SDSS>0),
   Radius  = Get.SDSS./(RAD.*3600);
   [Run,MJD,Dist]=coo2run(RA,Dec);
   
   if (isempty(Run{1})==1),
      Data.SDSS.Exist = 0;
   else
      Data.SDSS.Exist = 1;
   end
   Data.SDSS.Run = Run;
   Data.SDSS.MJD = MJD;
   if (Data.SDSS.Exist==1),
      [FC,Im,NaviFC] = get_sdss_finding_chart(RA,Dec);
      Data.SDSS.FC = NaviFC;
      switch Get.GetIm
       case 0
          I = 1;    % only one image
          Data.SDSS.Images = get_sdss_corrim(Run{I},'y','ugriz');
       case 1
          I = 1;    % only one image
          [Data.SDSS.Images,FN] = get_sdss_corrim(Run{I},'y','ugriz');
       otherwise
          error('Unknown Get.GetIm option');
      end
   
      if (Get.GetCat==1),   
         Q{1} = {'ra','dec','modelMag_u','modelMagErr_u','modelMag_g','modelMagErr_g','modelMag_r','modelMagErr_r','modelMag_i','modelMagErr_i','modelMag_z','modelMagErr_z','primTarget','Photoz.z','Photoz.zErr'};
         Q{2}={'PhotoPrimary','Photoz'};
         MinRA  = RA - Radius./cos(Dec);
         MaxRA  = RA + Radius./cos(Dec);
         MinDec = Dec - Radius;
         MaxDec = Dec + Radius;
         Q{3}=sprintf('(PhotoPrimary.ObjID=Photoz.ObjID) and (ra between %f and %f) and (dec between %f and %f)',MinRA.*RAD,MaxRA.*RAD,MinDec.*RAD,MaxDec.*RAD);
         
         [Cat,Msg] = run_sdss_sql(Q,[],'cell');
   
         Data.SDSS.CatHeader = Q{1};
         Data.SDSS.Cat       = Cat;

         Q{1} = {'ra','dec','modelMag_u','modelMagErr_u','modelMag_g','modelMagErr_g','modelMag_r','modelMagErr_r','modelMag_i','modelMagErr_i','modelMag_z','modelMagErr_z'};
         Q{2}={'PhotoPrimary'};
         MinRA  = RA - Radius./cos(Dec);
         MaxRA  = RA + Radius./cos(Dec);
         MinDec = Dec - Radius;
         MaxDec = Dec + Radius;
         Q{3}=sprintf('(ra between %f and %f) and (dec between %f and %f)',MinRA.*RAD,MaxRA.*RAD,MinDec.*RAD,MaxDec.*RAD);
         
         [Cat,Msg] = run_sdss_sql(Q,[],'cell');
   
         Data.SDSS.CatPhotHeader = Q{1};
         Data.SDSS.CatPhot       = Cat;
      
      
         % sdss spectroscopy
         load SDSS_DR5_AllSpec.mat
         [L,Dist,PA] = cat_search(SDSS_DR5_AllSpec(:,1:2)./RAD,[1 2],[RA Dec],Radius','circle','Dec','sphere');
         Data.SDSS.Spec = [SDSS_DR5_AllSpec(L,[1:2]),Dist,PA, SDSS_DR5_AllSpec(L,[3 6 7:13])];
         Data.SDSS.SpecHeader = {'RA','Dec','Dist','PA','MJD','Z','Zerr','Zconf','Zstatus','Zwarning','SpecClass','VelDisp','VelDispErr'};
      
         clear SDSS_DR5_AllSpec
      end
   end
end


%-------------
%--- GALEX ---
%-------------
if (Get.GALEX>0),
   Radius  = Get.GALEX./(RAD.*3600);
   load GALEX_Images_DR2p3.mat
   L = cat_search(GALEX_Images_DR2p3(:,1:2)./RAD,[1 2],[RA Dec],0.625./RAD,'circle');
   if (isempty(L)==0),
      Data.GALEX.Exist = 1;
   else   
      Data.GALEX.Exist = 0;
   end
end


%------------
%--- POSS ---
%------------
if (Get.POSS>0),
   switch Get.GetIm
    case 0
       Data.POSS.Im2R = get_dss(RA,Dec,DSS_FOV,'y','2R');
       Data.POSS.Im2B = get_dss(RA,Dec,DSS_FOV,'y','2B');
       Data.POSS.Im2I = get_dss(RA,Dec,DSS_FOV,'y','2I');
       Data.POSS.Im1R = get_dss(RA,Dec,DSS_FOV,'y','1R');
       Data.POSS.Im1B = get_dss(RA,Dec,DSS_FOV,'y','1B');
       Data.POSS.ImQV = get_dss(RA,Dec,DSS_FOV,'y','QV');
    case 1
       [Data.POSS.Im2R,FN] = get_dss(RA,Dec,DSS_FOV,'y','2R');
       [Data.POSS.Im2B,FN] = get_dss(RA,Dec,DSS_FOV,'y','2B');
       [Data.POSS.Im2I,FN] = get_dss(RA,Dec,DSS_FOV,'y','2I');
       [Data.POSS.Im1R,FN] = get_dss(RA,Dec,DSS_FOV,'y','1R');
       [Data.POSS.Im1B,FN] = get_dss(RA,Dec,DSS_FOV,'y','1B');
       [Data.POSS.ImQV,FN] = get_dss(RA,Dec,DSS_FOV,'y','QV');
    otherwise
       error('Unknown Get.GetIm option');
   end

end


%-------------
%--- 2MASS ---
%-------------
if (Get.MASS>0),
   Radius  = Get.MASS./(RAD.*3600);
   if (Get.GetCat==1),
      try
         [Data.MASS.Cat,Data.MASS.CatHeader] = search2mass(RA,Dec,Radius.*RAD.*60);
      catch
         Data.MASS.Cat = NaN;
         Data.MASS.CatHeader = NaN;
      end
   end
end

%-----------
%--- NED ---
%-----------
if (Get.NED>0),
   Radius  = Get.NED./60;
   SignHTML = {'-','%2B'};
   DecSign  = SignHTML{round(sign(Dec).*0.5+1.5)};

   Data.NED.Link = sprintf('http://nedwww.ipac.caltech.edu/cgi-bin/nph-objsearch?in_csys=Equatorial&in_equinox=J2000.0&lon=%fd&lat=%s%fd&radius=%f&search_type=Near+Position+Search&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=Distance+to+search+center&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES&z_constraint=Unconstrained&z_value1=&z_value2=&z_unit=z&ot_include=ANY&nmp_op=ANY',RA.*RAD,DecSign,abs(Dec).*RAD,Radius);

end


%-------------
%--- FIRST ---
%-------------
if (Get.FIRST>0),
   Radius  = Get.FIRST./(RAD.*3600);
   Data.FIRST.Link      = sprintf('http://third.ucllnl.org/cgi-bin/firstcutout?RA=%10.7f&Dec=%10.6f&ImageSize=1&MaxInt=1&ImageType=GIF',RA.*RAD./15,Dec.*RAD);

   if (Get.GetCat==1),   
      load first_030411.mat;
      [L,Dist,PA]          = cat_search(first_030411,[1 2],[RA Dec],Radius,'circle','Dec','sphere');
      Data.FIRST.Cat       = [first_030411(L,[1:2]), Dist,PA, first_030411(L,[3:12])];
      Data.FIRST.CatHeader = {'RA','Dec','Dist','PA','Sidelobes','PeakFlux','IntFlux','RMS','DeconMajor','DeconMinor','DeconPA','Major','Minor','PA'};
      clear first_030411;
   end

end



%------------
%--- NVSS ---
%------------
if (Get.NVSS>0),
   Radius  = Get.NVSS./(RAD.*3600);
   Data.NVSS.Link = sprintf('http://www.cv.nrao.edu/cgi-bin/postage.pl?RA=%10.7f&Dec=%10.6f&Size=0.083333+0.083333&Equinox=2000&Type=application/postscript',RA.*RAD./15,Dec.*RAD);

   if (Get.GetCat==1),
      load NVSS.mat
			       [L,Dist,PA] = cat_search(NVSS,[1 2],[RA Dec],Radius,'circle','Dec','sphere');

      Data.NVSS.Cat = NVSS(L,:);
      Data.NVSS.CatHeader = {'RA','Dec','ErrRA','ErrDec','Flux','FluxErr','Major','Minor','PA','ErrMajor','ErrMinor','ErrPA'};

      clear NVSS
   end
end


%--------------
%--- MaxBCG ---
%--------------
if (Get.MaxBCG>0),
   Radius  = Get.MaxBCG./(RAD.*3600);
   if (Get.GetCat==1),
      load MaxBCG.mat;
      [L,Dist,PA] = cat_search(MaxBCG,[1 2],[RA Dec],Radius,'circle','Dec','sphere');
      Data.MaxBCG.Cat = MaxBCG(L,[1:4, 9]);
      Data.MaxBCG.CatHeader = {'RA','Dec','photoz','z','Ngal'};
      clear MaxBCG;
   end
end


%-------------
%--- Abell ---
%-------------
if (Get.Abell>0),
   Radius  = Get.Abell./(RAD.*3600);
   if (Get.GetCat==1),
      load abell.mat
      [L,Dist,PA] = cat_search(abell,[3 4],[RA Dec],Radius,'circle','Dec','sphere');

      Data.Abell.Cat = abell(L,[1 3:7 15 22]);
      Data.Abell.CatHeader = {'N','RA','Dec','Mag10','DistGr','RichGr','Radius','z'};
      clear abell
   end
end


%-----------
%--- PGC ---
%-----------
if (Get.PGC>0),
   Radius  = Get.PGC./(RAD.*3600);
   if (Get.GetCat==1),
      load pgc.mat
      [L,Dist,PA] = cat_search(pgc,[1 2],[RA Dec],Radius,'circle','Dec','sphere');
      Data.PGC.Cat = pgc(L,:);
      Data.PGC.CatHeader = {'RA','Dec','Type','logD25','ErrlogD25','LogAxisRatio','ErrLogAxisRatio','PA','ErrPA','NumNames'};

      clear pgc
   end
end


%-------------
%--- ROSAT ---
%-------------
if (Get.ROSAT>0),
   Radius  = Get.ROSAT./(RAD.*3600);
   if (Get.GetCat==1),
      load rosat_faint_xs.mat
      [L,Dist,PA]=cat_search(rosat_faint_xs(:,1:2)./RAD,[1 2],[RA Dec],Radius,'circle','Dec','sphere');
      Data.ROSAT.FaintCat = rosat_faint_xs(L,[1 5]);
      Data.ROSAT.FaintCatHeader = {'RA','Dec','PosErr','CTS','ErrCTS'};
      clear rosat_faint_xs
   end
end


%-----------
%--- XMM ---
%-----------
if (Get.XMM>0),
   Radius  = Get.XMM./(RAD.*3600);
   if (Get.GetCat==1),
      load XMM_SSC2.mat
      [L,Dist,PA] = cat_search(XMM_SSC2,[1 2],[RA Dec],Radius,'circle','Dec','sphere');
      Data.XMM.Cat = XMM_SSC2(L,:); 
      Data.XMM.CatHeader = {'RA','Dec','PosErr','MJD','ExpTime','Flux','FluxErr','HR1','ErrHR1','HR2','ErrHR2','HR3','ErrHR3','HR4','ErrHR4'};
      clear XMM_SSC2
   end
end



%-----------
%--- SNR ---
%-----------
if (Get.SNR>0),
   Radius  = Get.SNR./(RAD.*3600);
   if (Get.GetCat==1),

   end
end



%--------------
%--- Pulsar ---
%--------------
if (Get.Pulsar>0),
   Radius  = Get.Pulsar./(RAD.*3600);
   if (Get.GetCat==1),

   end
end


%-----------
%--- HST ---
%-----------
if (Get.HST>0),
   Radius  = Get.HST./(RAD.*3600);
   if (Get.GetCat==1),

   end
end
