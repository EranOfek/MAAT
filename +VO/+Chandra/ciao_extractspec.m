function [Cnt]=ciao_extractspec(EvtFile,RA,Dec,varargin)
% Prepare the ARF and MRF Chandra files required for X-ray spectroscopy.
% Package: VO.Chandra
% Description: Use Chandra/CIAO tasks to prepare the ARF and MRF files
%              required for X-ray spectroscopy of a source.
%              Instellation: ciao is required.
% Input  : - The Chandra level 2 event file.
%          - J2000.0 RA of the source for which the spectrum should be
%            extracted.
%          - J2000.0 Dec of the source for which the spectrum should be
%            extracted.
%          * Arbitrary pairs of ...,key,val,... arguments.
%            The following keywords are available:
%            'Aper'    - Object extraction radius [pixels]. Default is 2.5.
%            'Annulus' - sky annulus in for which to calculate the
%                        ARF/MRF files. Default is [40 100] pixels.
%            'CooSys'  - Coordinate system of source {'eq','xy'}.
%                        Default is 'eq'.
%            'OutRoot' - Output directory. Default is 'spec'.
% Output : - Number of counts within aperture.
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Mar 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: VO.Chandra.ciao_extractspec('acisf11122N002_evt2.fits',X,Y,'CooSys','XY')
% Reliable: 2 (Modification were not tested)
%--------------------------------------------------------------------------



% start FTOOLS/HEASOFT/CIAO
CIAO_CSH = 'source /home/eran/bin/ciao/ciao-4.5/bin/ciao.csh';

DefV.Aper          = 2.5;
DefV.Annulus       = [40 100];
DefV.CooSys        = 'eq';
DefV.ObjGroupNum   = 15;
DefV.BackGroupNum  = 100;
DefV.Prefix        = '';
DefV.OutRoot       = 'spec';
DefV.Weight        = 'no';
DefV.Correct       = 'yes';
DefV.GroupType     = 'NUM_CTS';
DefV.SecondaryDir  = '../secondary/';
DefV.UseSpecExtract= 'y'; %
DefV.EnergyGrid    = '0.3:11.0:0.01';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);



% Instructions: http://asc.harvard.edu/ciao/threads/pointlike/index.html

Verbose = 0;
Clobber = 'yes';

%EvtFile = 'acisf11122N002_evt2.fits'
%RA      = '09:42:53.3'
%Dec     = '+09:29:42.0'
Annulus = [40 80];


% converty RA/Dec to pixel coordinates
switch lower(InPar.CooSys)
    case 'eq'
        RA  = celestial.coo.convertdms(RA,'gH','r');
        Dec = celestial.coo.convertdms(Dec,'gD','R');

        [X,Y] = swxrt_coo2xy(EvtFile,RA,Dec,'chandra');
    case 'xy'
        X = RA;
        Y = Dec;
    otherwise
        error('unknown CooSys option');
end


ObjFITS       = sprintf('%s%s',InPar.Prefix,'Object.pi');
BackFITS      = sprintf('%s%s',InPar.Prefix,'Object_bg.pi');
ObjectRMF     = sprintf('%s%s',InPar.Prefix,'Object_mkacisrmf.rmf');
BackRMF       = sprintf('%s%s',InPar.Prefix,'Object_bg_mkacisrmf.rmf');
ObjectARF     = sprintf('%s%s',InPar.Prefix,'Object.arf');
ObjectAspHist = sprintf('%s%s',InPar.Prefix,'Object.asphist');
BackARF       = sprintf('%s%s',InPar.Prefix,'Object_bg.arf');
ObjectGRP     = sprintf('%s%s',InPar.Prefix,'Object_grp.pi');
BackGRP       = sprintf('%s%s',InPar.Prefix,'Object_bg_grp.pi');


switch lower(InPar.UseSpecExtract)
    case 'y'
       %--- Use the specextract tool ---
      Command = 'punlearn specextract';
      system(sprintf('%s; %s',CIAO_CSH,Command));


      FilesE = dir('pcad*asol1.fits');
      ASolFile = FilesE.name;

      FilesE = dir(sprintf('%s*msk1.fit*',InPar.SecondaryDir));
      MskFile  = FilesE.name;

      FilesE = dir('*bpix1.fit*');
      BpixFile = FilesE.name;

      FilesE = dir(sprintf('%s*pbk0.fit*',InPar.SecondaryDir));
      PbkFile  = FilesE.name;

      Command = sprintf('specextract verbose=%d infile="%s[sky=circle(%08.2f,%08.2f,%f)]" outroot=%s bkgfile="%s[sky=annulus(%08.2f,%08.2f,%d,%d)]" weight=%s correct=%s asp=%s mskfile=%s%s badpixfile=%s pbkfile=%s%s grouptype=%s binspec=%d clobber=%s combine=%s bkgresp=%s',...
                   Verbose,EvtFile,X,Y,InPar.Aper,InPar.Prefix,...
                   EvtFile,X,Y,InPar.Annulus(1),InPar.Annulus(2),...
                   InPar.Weight,InPar.Correct,ASolFile,...
                   InPar.SecondaryDir,MskFile,...
                   BpixFile,...
                   InPar.SecondaryDir,PbkFile,...
                   InPar.GroupType,InPar.ObjGroupNum,...
                   Clobber,'no','yes');

      system(sprintf('%s; %s',CIAO_CSH,Command));

      % verify that the energy grids in ARF and RMF files are identical so that
      % we can use them in XSPEC

      [Table,Col] = FITS.read_table(EvtFile);
      %[~,~,~,Col,Table]=get_fitstable_col(EvtFile);
      D=Util.Geom.plane_dist(Table{Col.x},Table{Col.y},X,Y);
      I=find(D<InPar.Aper);
      Cnt = length(I);
      if (1==0)

      UnCCDID = unique(Table{Col.ccd_id}(I));
      if (length(UnCCDID)~=1)
         error('Near CCD edge - abort');
      end
    
      ChipX = mean(Table{Col.chipx}(I));
      ChipY = mean(Table{Col.chipy}(I));

      CCDmapping = {'0','ACIS-I0';'1','ACIS-I1';'2','ACIS-I2';'3','ACIS-I4';'4','ACIS-S0';'5','ACIS-S1';'6','ACIS-S2';'7','ACIS-S3';'8','ACIS-S4';'9','ACIS-S5'};
      I=find(Util.cell.isempty_cell(strfind(CCDmapping(:,1),num2str(UnCCDID)))==0);
      DetSubSystem = CCDmapping{I,2};


      Command = 'punlearn dmhistory';
      system(sprintf('%s; %s',CIAO_CSH,Command));
               
      Command = 'dmhistory spec.arf tool=all';
      system(sprintf('%s; %s',CIAO_CSH,Command));
                   
      Command = 'punlearn mkarf';
      system(sprintf('%s; %s',CIAO_CSH,Command));
                
      Command = sprintf('mkarf asphistfile="%s" outfile="%s" sourcepixelx="%f" sourcepixely="%f" engrid="%s" obsfile="%s" pbkfile="%s%s" dafile="%s" mirror="%s" detsubsys="%s" grating="NONE" maskfile="%s%s" ardlibparfile="%s" geompar="geom" verbose=%d clobber=%s',...
          ObjectAspHist,ObjectARF,X,Y,...   
          InPar.EnergyGrid,EvtFile,...
          InPar.SecondaryDir,PbkFile,...
          'CALDB',...
          'HRMA',...
          DetSubSystem,...
          InPar.SecondaryDir,MskFile,...
          'ardlib.par',...
          Verbose,Clobber);
      
      system(sprintf('%s; %s',CIAO_CSH,Command));

          
      Command = 'dmhistory spec.rmf tool=all';
      system(sprintf('%s; %s',CIAO_CSH,Command));
                   
      Command = 'punlearn mkacisrmf';
      system(sprintf('%s; %s',CIAO_CSH,Command));
               
      Command = sprintf('mkacisrmf infile="%s" outfile="%s" wmap="%s" energy="%s" channel"%s" chantype="%s" ccd_id="%d" chipx="%f" chipy="%f" gain="%s" asolfile="" obsfile="%s" logfile="" contlvl="%d" geompar="%s" thresh="%f" clobber=%s verbose=%d',...
                'CALDB',ObjectRMF,'none',InPar.EnergyGrid,...
                '1:1024:1','PI',UnCCDID,ChipX,ChipY,...
                'CALDB',EvtFile,100,'geom',1e-6,...
                Clobber,Verbose);
      system(sprintf('%s; %s',CIAO_CSH,Command));
      end
    case 'n'
        
%--- dmextract object ---
Command = 'punlearn dmextract';
system(sprintf('%s; %s',CIAO_CSH,Command));

Command = sprintf('dmextract verbose=%d infile="%s[sky=circle(%08.2f,%08.2f,%f)][bin pi]" outfile=%s clobber=%s',...
                   Verbose,EvtFile,X,Y,InPar.Aper,ObjFITS,Clobber);
system(sprintf('%s; %s',CIAO_CSH,Command));



%--- dmextract background ---
Command = 'punlearn dmextract';
system(sprintf('%s; %s',CIAO_CSH,Command));

Command = sprintf('dmextract verbose=%d infile="%s[sky=annulus(%08.2f,%08.2f,%d,%d)][bin pi]" outfile=%s clobber=%s',...
                   Verbose,EvtFile,X,Y,InPar.Annulus(1),InPar.Annulus(2),BackFITS,Clobber);
system(sprintf('%s; %s',CIAO_CSH,Command));


% get the chipx/chipy coordinates
[Table,Col] = FITS.read_table(EvtFile);
%[~,~,~,Col,Table]=get_fitstable_col(EvtFile);
D=Util.Geom.plane_dist(Table{Col.x},Table{Col.y},X,Y);
I=find(D<InPar.Aper);
length(I)

UnCCDID = unique(Table{Col.ccd_id}(I));
if (length(UnCCDID)~=1)
    error('Near CCD edge - abort');
end
    
ChipX = mean(Table{Col.chipx}(I));
ChipY = mean(Table{Col.chipy}(I));

DirCALDB = '/home/eran/bin/ciao/ciao-4.5/CALDB';

KeyVal = FITS.get_keys(EvtFile,{'gainfile'},1);
%KeyVal = get_fits_keyword(EvtFile,{'gainfile'},1,'BinaryTable');
GainFile = KeyVal{1};

GainFile = sprintf('%s%s%s%s',DirCALDB,filesep,'data/chandra/acis/det_gain/',GainFile);


% choose latest response file
DirResponse  = sprintf('%s%s%s',DirCALDB,filesep,'data/chandra/acis/p2_resp');
File = dir(sprintf('%s%s*.fits',DirResponse,filesep));
[~,Im]=max([File.datenum]);
ResponseFile = sprintf('%s%s%s',DirResponse,filesep,File(Im).name);

%--- run mkacisrmf for object ---

Command = 'punlearn mkacisrmf';
system(sprintf('%s; %s',CIAO_CSH,Command));


Energy = [0.25 10.0 0.01];  % min max binning_step

Command = sprintf('mkacisrmf infile=%s outfile=%s energy=%f:%f:%f channel=1:1024:1 chantype=PI wmap=none ccd_id=%d chipx=%f chipy=%f gain=%s clobber=%s',...
    ResponseFile,ObjectRMF,Energy(1),Energy(2),Energy(3),UnCCDID,ChipX,ChipY,GainFile,Clobber);
system(sprintf('%s; %s',CIAO_CSH,Command));


%--- run mkacisrmf for background ---

I=find(D>InPar.Annulus(1) & D<InPar.Annulus(2));
length(I)
UnCCDIDb = unique(Table{Col.ccd_id}(I));
if (length(UnCCDIDb)~=1)
    error('Near CCD edge - abort');
end
    
ChipXb = mean(Table{Col.chipx}(I));
ChipYb = mean(Table{Col.chipy}(I));

Command = sprintf('mkacisrmf infile=%s outfile=%s energy=%f:%f:%f channel=1:1024:1 chantype=PI wmap=none ccd_id=%d chipx=%f chipy=%f gain=%s clobber=%s',...
    ResponseFile,BackRMF,Energy(1),Energy(2),Energy(3),UnCCDIDb,ChipXb,ChipYb,GainFile,Clobber);
system(sprintf('%s; %s',CIAO_CSH,Command));

%--- run asphist (aspect histogram) ---
Command = 'punlearn asphist';
system(sprintf('%s; %s',CIAO_CSH,Command));

FilesE = dir('pcad*asol1.fits');
ASolFile = FilesE.name;


Command = sprintf('asphist infile=%s outfile=%s evtfile="%s[ccd_id=%d]" clobber=%s',...
        ASolFile,ObjectAspHist,EvtFile,UnCCDID,Clobber);
system(sprintf('%s; %s',CIAO_CSH,Command));

%--- run mkarf for source --- 
Command = 'punlearn mkarf';
system(sprintf('%s; %s',CIAO_CSH,Command));

CCDmapping = {'0','ACIS-I0';'1','ACIS-I1';'2','ACIS-I2';'3','ACIS-I4';'4','ACIS-S0';'5','ACIS-S1';'6','ACIS-S2';'7','ACIS-S3';'8','ACIS-S4';'9','ACIS-S5'};
I=find(Util.cell.isempty_cell(strfind(CCDmapping(:,1),num2str(7)))==0);
DetSubSystem = CCDmapping{I,2};


% look for mask file
FileM = dir('../secondary/*_msk*.fits*');

% look for dead area file (pbk file)
FileP = dir('../secondary/*_pbk*.fits*');

Command = sprintf('mkarf verbose=%d grating=NONE detsubsys=%s outfile=%s asphistfile="%s[ASPHIST]" obsfile="%s[EVENTS]" dafile=CALDB engrid="grid(%s[cols ENERG_LO,ENERG_HI])" sourcepixelx=%f sourcepixely=%f maskfile="../secondary/%s" pbkfile="../secondary/%s" clobber=%s',...
                  Verbose,DetSubSystem,ObjectARF,ObjectAspHist,EvtFile,ObjectRMF,X,Y,FileM(1).name,FileP(1).name,Clobber);
system(sprintf('%s; %s',CIAO_CSH,Command));

%--- run mkarf for background --- 
Command = 'punlearn mkarf';
system(sprintf('%s; %s',CIAO_CSH,Command));


Command = sprintf('mkarf verbose=%d grating=NONE detsubsys=%s outfile=%s asphistfile="%s[ASPHIST]" obsfile="%s[EVENTS]" dafile=CALDB engrid="grid(%s[cols ENERG_LO,ENERG_HI])" sourcepixelx=%f sourcepixely=%f maskfile="../secondary/%s" pbkfile="../secondary/%s" clobber=%s',...
                  Verbose,DetSubSystem,BackARF,ObjectAspHist,EvtFile,ObjectRMF,X,Y,FileM(1).name,FileP(1).name,Clobber);
system(sprintf('%s; %s',CIAO_CSH,Command));

%--- run arfcorr ---
%Command = 'punlearn arfcorr';
%system(sprintf('%s; %s',CIAO_CSH,Command));
%
%RegionFile = 'Object_Region.reg';
%Command = sprintf('dmmakereg region="circle(%f,%f,%f)" outfile=%s clobber=%s',...
%        X,Y,InPar.Aper,RegionFile,Clobber);
%system(sprintf('%s; %s',CIAO_CSH,Command));

%ObjectARFcorr = 'Object_corr.arf';
%Command = sprintf('arfcorr infile="%s[sky=circle(%08.2f,%08.2f,%f)][bin pi]" arf=%s outfile=%s x=%f y=%f region="region(%s)" verbose=%d clobber=%s',...
%                  EvtFile,X,Y,InPar.Aper,ObjectARF,ObjectARFcorr,X,Y,RegionFile,Verbose,Clobber);
%              Command
%system(sprintf('%s; %s',CIAO_CSH,Command));



%--- run dmgroup for source ---
Command = 'punlearn dmgroup';
system(sprintf('%s; %s',CIAO_CSH,Command));

Command = sprintf('dmgroup binspec="" infile=%s outfile=%s grouptype=NUM_CTS grouptypeval=%d xcolumn=channel ycolumn=counts verbose=%d clobber=%s',...
                  ObjFITS,ObjectGRP,InPar.ObjGroupNum,Verbose,Clobber);
system(sprintf('%s; %s',CIAO_CSH,Command));


%--- run dmgroup for background ---
Command = 'punlearn dmgroup';
system(sprintf('%s; %s',CIAO_CSH,Command));

Command = sprintf('dmgroup binspec="" infile=%s outfile=%s grouptype=NUM_CTS grouptypeval=%d xcolumn=channel ycolumn=counts verbose=%d clobber=%s',...
                  ObjFITS,BackGRP,InPar.BackGroupNum,Verbose,Clobber);
system(sprintf('%s; %s',CIAO_CSH,Command));

%--- update headers ---
Command = sprintf('punlearn dmhedit');
system(sprintf('%s; %s',CIAO_CSH,Command));


InFile = ObjFITS;
Key    = 'BACKFILE';
Value  = BackFITS;
Command = sprintf('dmhedit infile=%s filelist="" operation=add key=%s value=%s',...
                  InFile,Key,Value);
system(sprintf('%s; %s',CIAO_CSH,Command));

InFile = ObjFITS;
Key    = 'RESPFILE';
Value  = ObjectRMF;
Command = sprintf('dmhedit infile=%s filelist="" operation=add key=%s value=%s',...
                  InFile,Key,Value);
system(sprintf('%s; %s',CIAO_CSH,Command));

InFile = ObjFITS;
Key    = 'ANCRFILE';
Value  = ObjectARF;
Command = sprintf('dmhedit infile=%s filelist="" operation=add key=%s value=%s',...
                  InFile,Key,Value);
system(sprintf('%s; %s',CIAO_CSH,Command));


InFile = ObjectGRP;
Key    = 'BACKFILE';
Value  = BackFITS;
Command = sprintf('dmhedit infile=%s filelist="" operation=add key=%s value=%s',...
                  InFile,Key,Value);
system(sprintf('%s; %s',CIAO_CSH,Command));

InFile = ObjectGRP;
Key    = 'RESPFILE';
Value  = ObjectRMF;
Command = sprintf('dmhedit infile=%s filelist="" operation=add key=%s value=%s',...
                  InFile,Key,Value);
system(sprintf('%s; %s',CIAO_CSH,Command));

InFile = ObjectGRP;
Key    = 'ANCRFILE';
Value  = ObjectARF;
Command = sprintf('dmhedit infile=%s filelist="" operation=add key=%s value=%s',...
                  InFile,Key,Value);
system(sprintf('%s; %s',CIAO_CSH,Command));


InFile = BackFITS;
Key    = 'RESPFILE';
Value  = BackRMF;
Command = sprintf('dmhedit infile=%s filelist="" operation=add key=%s value=%s',...
                  InFile,Key,Value);
system(sprintf('%s; %s',CIAO_CSH,Command));

InFile = BackFITS;
Key    = 'ANCRFILE';
Value  = BackARF;
Command = sprintf('dmhedit infile=%s filelist="" operation=add key=%s value=%s',...
                  InFile,Key,Value);
system(sprintf('%s; %s',CIAO_CSH,Command));

InFile = BackGRP;
Key    = 'RESPFILE';
Value  = BackRMF;
Command = sprintf('dmhedit infile=%s filelist="" operation=add key=%s value=%s',...
                  InFile,Key,Value);
system(sprintf('%s; %s',CIAO_CSH,Command));

InFile = BackGRP;
Key    = 'ANCRFILE';
Value  = BackARF;
Command = sprintf('dmhedit infile=%s filelist="" operation=add key=%s value=%s',...
                  InFile,Key,Value);
system(sprintf('%s; %s',CIAO_CSH,Command));



    otherwise
        errorr('Unknown UseSpecExtract option');
end
