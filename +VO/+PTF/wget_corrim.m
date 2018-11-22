function  [Link,TableSort,Col,FileName]=wget_corrim(RA,DEC,varargin)
% Get PTF corrected images from the IPAC archive
% Package: VO.PTF
% Description:  Given coordinates or PTF field+CCDID, get the link to the
%               images and associated catalogs and mask files. The program
%               can also retrieve the files.
% Input  : - J2000.0 RA in [rad] or [H M S] or sexagesimal string.
%            Alternatively, if the second argument is an empty matrix then

%            then the first input argument will be regarded as
%            [FieldID CCDID].
%          - J2000.0 Dec in [rad] or [Sign D M S] or sexagesimal string. 
%            Alternatively, if this is an empty matrix then the first
%            argument is [FieldID CCDID].
%            Spatial queries (search by RA,Dec): return images that are
%            related to a search region (specified via the POS and SIZE
%            parameters) according to a spatial predicate (specified via
%            the INTERSECT parameter).        
%          * Arbitrary number of pairs of ...,keyword,value,...
%            The following keywords are available:
%            'SIZE'  - [width[,height]]
%                      (e.g., 'SIZE', 0.1  or  'SIZE', [0.1, 0.2]).
%                      The SIZE parameter consists of one or two values in
%                      decimal degrees. The first SIZE value is taken to be
%                      the full-width of the search region along the east 
%                      axis at POS, and the second is taken to be the
%                      full-height along the north axis. 
%                      If only one size value is specified, it is used as
%                      both the full-width and full-height. Negative sizes
%                      are illegal, and a width and height of zero indicate
%                      that the search region is a point.
%                      Default is 0.
%                      If INTERSECT=CENTER (see below), SIZE is ignored.
%             'Inetesect'- Specify the spatial search method.
%                          Allowed values are:
%                          'Covers' - X must completely contain S. 
%                                     Equivalent to CENTER and OVERLAPS
%                                     if S is a point.
%                          'Enclosed'-S must completely contain X. If S is
%                                     a point, the query will always return
%                                     an empty image table.
%                          'Center'  -X must contain the center of S. If S
%                                     is a point, this is equivalent to
%                                     COVERS and OVERLAPS.
%                          'Overlaps'-The intersection of S and X is
%                                     non-empty. If S is a point, this is
%                                     equivalent to CENTER and COVERS.
%                          Default is 'Overlaps'.
%             'JD' - Return results only in the specified JD range
%                    [min_JD, max_JD]. If single value is given then
%                    max_JD will be set to Inf.       
%             'Filter' - Cell array of filters. Return only specified
%                        filters. Default is all filters.
%             'Seeing' - Return only images which seeing is in the
%                        specified range [minSeeing, maxSeeing].
%                        If single value is given then minSeeing is 0.
%                        Default is [0 Inf].
%             'Save'   - Combination of chracter that specifiy which files
%                        to save.
%                        'i' - corrected images.
%                        'm' - mask files.
%                        'c' - catalogs.
%                        'p' - psf phot catalog
%                        'none' - do not save.
%                        e.g., 'i' or 'im'. Default is 'ic'.
%             'MaxGet' - Number of files to download in parallel.
%                        Default is 5.
%             'Extra'  - A string of extra parameters to pass to the wget
%                        command. See www.pwget.m for options.
%                        Default is '-q -nc'.
%             'User'   - User name, or cell array containg path for
%                        user/pass information.
%             'Pass'   - User password. If empty use user/pass info.
%                        Default is empty.
% Output : - Struct array of all the files found.
%            The following fields are available:
%            .JD      - Julian day
%            .seeing  - seeing
%            .filter  - filter
%            .Image   - Link to image.
%            .Mask    - Link to mask file.
%            .Cat     - Link to catalog
%          - A cell array of the IPAC table returned by specified search.
%          - Structure in which each field contains the column index in
%            the IPAC table.
%          - A structure with the following fields:
%            .Files_Image - Cell array of downloaded images.
%            .Files_Mask  - Cell array of downloaded masks.
%            .Files_Cat   - Cell array of downloaded catalogs.
% Tested : Matlab R2013b
%     By : Yifat Dzigan                    Jan 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Examples: 
%  %To retrive the Table and the Link to the images run with Save','none':
%   [Link,TableSort,Col,Filename]=VO.PTF.wget_corrim([100019 0],[], 'Save','none');
%  %To search with Ra, Dec: 
%   [Link,TableSort,Col,Filename]=VO.PTF.wget_corrim(RA,DEC,'Filter',{'ha663','ha656'},'intersect','OVERLAPS','SIZE',[0.2,0.1],'save','im', 'MaxGet',15)
% Reliable: 2
%--------------------------------------------------------------------------




RAD             = 180./pi;
Narg            = length(varargin);
ImTable         = tempname;

DefV.Save       = 'ic'; %or only 'i', 'm','c' or 'none'
DefV.SIZE       = 0;
DefV.Intersect ='OVERLAPS';
DefV.Seeing     = [0,inf];
DefV.JD         = [-inf,inf];
DefV.Filter     = {'g', 'r', 'ha656', 'ha663'};
DefV.MaxGet     = 5;
DefV.Extra      = '-nc -v'; % '-q -nc'; % Option to supply the extra wget switch parameters
DefV.User       = {'~/matlab/passwords/ptf_archive_pass'};
DefV.Pass       = [];

%InPar = set_varargin_keyval(DefV,'y','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (iscell(InPar.User))
    [InPar.User,InPar.Pass]=Util.code.read_user_pass_file(InPar.User{1});
end

Extrac = sprintf('--save-cookies=ptf.txt -O /dev/null');  % -O =>  write to file, If file already exists, it will be overwritten. 
%str = sprintf('http://irsa.ipac.caltech.edu/account/signon/login.do?josso_cmd=login&josso_username=eran@astro.caltech.edu&josso_password=A@anohma');
str = sprintf('http://irsa.ipac.caltech.edu/account/signon/login.do?josso_cmd=login&josso_username=%s&josso_password=%s',...
    InPar.User,InPar.Pass);

www.pwget({str},Extrac,1);



% Search and Download 
Fid             = [];
CCDid           = [];
if (isempty(DEC))
    % Search by Field and CCD id  %In this case the table has only 39 Col, arranged differently!
    Fid   = RA(1);      % Enter a PTF field id to search, format example: 4419
    CCDid = RA(2);      % CCDid in the range 1-11 
    %InPar = set_varargin_keyval(DefV,'y','use',varargin{:});
end

Extra1 = sprintf('--load-cookies=ptf.txt -O %s',ImTable);       %obligatory, -O =>  write to file, If file already exists, it will be overwritten.
if (~isempty(InPar.Extra))
    Extra1 = sprintf('%s %s',InPar.Extra,Extra1);               %combine 
end

if (~isempty(Fid))
    % get the table of images
  %  system(sprintf('wget %s"http://kanaloa.ipac.caltech.edu/ibe/search/ptf/dev/process?where=ptffield%%3D%%27%g%%27+AND+ccdid%%3D%%27%g%%27" -O %s', url_pass,Fid,CCDid,ImTable));
%    system(sprintf('wget %s"http://kanaloa.ipac.caltech.edu/ibe/search/ptf/dev_ims/ptf_l1_img?where=ptffield%%3D%%27%g%%27+AND+ccdid%%3D%%27%g%%27" -O %s', url_pass,Fid,CCDid,ImTable));
    str = sprintf('http://irsa.ipac.caltech.edu/ibe/search/ptf/images/level1?where=ptffield%%3D%%27%g%%27+AND+ccdid%%3D%%27%g%%27',Fid,CCDid);
    
    www.pwget({str},Extra1);
    %%%%system(sprintf('wget --load-cookies=ptf.txt "http://irsa.ipac.caltech.edu/ibe/search/ptf/images/level1?where=ptffield%%3D%%27%g%%27+AND+ccdid%%3D%%27%g%%27" -O %s',Fid,CCDid,ImTable));    
else
    RA  = celestial.coo.convertdms(RA, 'gH','r');
    DEC = celestial.coo.convertdms(DEC,'gD','R');
    %InPar = set_varargin_keyval(DefV,'y','use',varargin{:});
   
    InPar.Intersect = upper(InPar.Intersect);
    
    if length(InPar.SIZE)==2
            S_width  = InPar.SIZE(1);
            S_height = InPar.SIZE(2);
    else
            S_width  = InPar.SIZE(1);
            S_height = S_width;
    end
    
    % get the table of images
%    system(sprintf('wget %s"http://kanaloa.ipac.caltech.edu/ibe/search/ptf/dev_ims/ptf_l1_img?POS=%g,%g&SIZE=%g,%g&INTERSECT=%s" -O %s', url_pass,RA*RAD,DEC*RAD,S_width,S_height,InPar.Intersect,ImTable));
    str = sprintf('http://irsa.ipac.caltech.edu/ibe/search/ptf/images/level1?POS=%g,%g&SIZE=%g,%g&INTERSECT=%s',RA*RAD,DEC*RAD,S_width,S_height,InPar.Intersect);
    www.pwget({str},Extra1,InPar.MaxGet);
    %%%system(sprintf('wget --load-cookies=ptf.txt "http://irsa.ipac.caltech.edu/ibe/search/ptf/images/level1?POS=%g,%g&SIZE=%g,%g&INTERSECT=%s" -O %s',RA*RAD,DEC*RAD,S_width,S_height,InPar.Intersect,ImTable));    
end

%if (~isempty(ImTable)),
FIDim=fopen(ImTable);
Stop=false;
I=0;
while (~Stop)
    I=I+1;
    Line1 = fgetl(FIDim);
    if ~strcmp(Line1(1),'\')
        Stop = true;
    end
end
%Line1=fgetl(FIDim);
Line2=fgetl(FIDim);                                             

RegExp1=regexp(Line1,'\|','split');
RegExp1=RegExp1(2:end-1);
RegExp1=strtrim(RegExp1);

RegExp2=regexp(Line2,'\|','split');
RegExp2=RegExp2(2:end-1);
RegExp2=strtrim(RegExp2);

[~,s1]=regexp(Line1,'\|','match');
fclose(FIDim);



Format = '';
for I1=1:1:length(s1)-1
    switch lower(RegExp2{I1})
       case{'int','double','long'} 
           Format='%f';
           ColReg(I1)=cellstr(sprintf('%s',lower(RegExp1{I1})));
       case{'char'}
           Format='%s';
           ColReg(I1)=cellstr(sprintf('%s',lower(RegExp1{I1})));
    
    end
    F      =  {s1(I1)+1, s1(I1+1)-1, Format};
    Table{:,I1} =  Util.IO.read_formatted(ImTable,F,'comment','\');
end
delete(ImTable);  % delete temp file

for I2=1:1:length(Table)
    TableCell{I2}=Table{I2}(5:end);
end
Col = cell2struct(num2cell([1:length(ColReg)]),ColReg,2);

FileName = [];
if (length(TableCell{Col.obsdate}))==0
    TableSort   = TableCell;
    Link        = struct([]);
    Files_Image =[];
    Files_Mask  =[];
    Files_Cat   =[];

else
    %--- convert to JD ---
    DateVec = datevec(TableCell{Col.obsdate},'yyyy-mm-ddHH:MM:SS');
    JulDay  = celestial.time.julday(DateVec(:,[3 2 1 4 5 6]));
    TableCell{Col.obsdate} = JulDay;
    
    %--- Define values to sort by, min, max (seeing, obsdate-by JD) ---
    if length(InPar.Seeing)==1
        InPar.Seeing=[0, InPar.Seeing];
    end

    if length(InPar.JD)==1
       InPar.JD=[InPar.JD, inf];
    end

    %--- Sort by Filter, seeing, obsdate(JD) ---
    Filter=cellstr(TableCell{Col.filter});
    SetI=[];
    FilterNum=length(InPar.Filter);
    for J=1:1:FilterNum
        FilterName = upper(InPar.Filter(J));
        tempI      = find(strcmp(Filter,FilterName) & cell2mat(TableCell{Col.seeing})>min(InPar.Seeing) & cell2mat(TableCell{Col.seeing})<max(InPar.Seeing) & TableCell{Col.obsdate}>min(InPar.JD) & TableCell{Col.obsdate}<max(InPar.JD));
        SetI       = [SetI; tempI];
    end
    SetI=sort(SetI);
    Files_Image =[];
    Files_Mask  =[];
    Files_Cat   =[];

    if isempty(SetI)==1
        disp('Sorting parameters do not apply to catalog search');
        TableSort =[];
        Link      = struct([]);
        FileName=struct([]);
     
    else
        for I4=1:1:length(TableCell)
          TableSort{I4}=TableCell{I4}(SetI);
        end

        %--- END: Sort by Filter, seeing, obsdate(JD) ---
    
        % Link = 1x4 struct array with fields:
        %    JD
        %    seeing
        %    filter
        %    Image (*.fits)
        %    Mask
        %    Cat
        %    CatPSF
        
        BaseURL='http://irsa.ipac.caltech.edu/ibe/data/ptf/images/level1/';
        
        for j1=1:1:length(TableSort{Col.obsdate})
           L_obsdate{j1}=TableSort{Col.obsdate}(j1);
        end
        Link_length=length(L_obsdate);
        [Link(1:Link_length).JD]=deal(L_obsdate{:});

        L_seeing=TableSort{Col.seeing};
        [Link(1:Link_length).seeing]=deal(L_seeing{:}); % arrange the structure!

        L_filter=TableSort{Col.filter};
        [Link(1:Link_length).filter]=deal(L_filter{:}); % Link.filter

        L_pfilename=TableSort{Col.pfilename};
        [Link(1:Link_length).Image]=deal(L_pfilename{:}); 

        L_Image=TableSort{Col.pfilename};
        L_Image=strcat(BaseURL,L_Image); 
        [Link(1:Link_length).Image]=deal(L_Image{:}); 

        L_Mask=TableSort{Col.afilename1};
        L_Mask=strcat(BaseURL,L_Mask);
        [Link(1:Link_length).Mask]=deal(L_Mask{:}); 

        L_Cat=TableSort{Col.afilename3};
        L_Cat=strcat(BaseURL,L_Cat);
        [Link(1:Link_length).Cat]=deal(L_Cat{:}); 

        L_CatPSF=TableSort{Col.afilename4};
        L_CatPSF=strcat(BaseURL,L_CatPSF);
        [Link(1:Link_length).CatPSF]=deal(L_CatPSF{:}); 
        
       % BaseURL=sprintf('%s http://kanaloa.ipac.caltech.edu/ibe/data/ptf/dev/process/',url_pass);
         %BaseURL='http://irsa.ipac.caltech.edu/ibe/data/ptf/images/level1/';
         Extra2 = sprintf('--load-cookies=ptf.txt');
         if (~isempty(InPar.Extra))
             Extra2 = sprintf('%s %s',InPar.Extra,Extra2);               %combine 
         end
         
       % BaseURL=sprintf('%s http://kanaloa.ipac.caltech.edu/ibe/data/ptf/dev_ims/ptf_l1_img/',url_pass);
                              
        switch isempty(strfind(lower(InPar.Save),'i'))
            case 0                
               FileName.Files_Image=www.pwget({Link.Image},Extra2,InPar.MaxGet); %,BaseURL);
               %%%%FileName.Files_Image=system(sprintf('wget --load-cookies=ptf.txt "http://irsa.ipac.caltech.edu/ibe/data/ptf/images/level1/Link(i).Image"'));
            otherwise
               % do nothing
        end
        switch isempty(strfind(lower(InPar.Save),'m'))
               case 0
                FileName.Files_Mask=www.pwget({Link.Mask},Extra2,InPar.MaxGet); %,BaseURL);
                %%%%FileName.Files_Mask=system(sprintf('wget --load-cookies=ptf.txt "http://irsa.ipac.caltech.edu/ibe/data/ptf/images/level1/Link(i).Mask"'));
             otherwise
                 % do nothing
         end
         switch isempty(strfind(lower(InPar.Save),'c'))
             case 0
                FileName.Files_Cat=www.pwget({Link.Cat},Extra2,InPar.MaxGet); %,BaseURL);
                %%%%FileName.Files_Cat=system(sprintf('wget --load-cookies=ptf.txt "http://irsa.ipac.caltech.edu/ibe/data/ptf/images/level1/Link(i).Cat"'));

             otherwise
                 % do nothing
         end
         
         switch isempty(strfind(lower(InPar.Save),'p'))
             case 0
                FileName.Files_CatPSF=www.pwget({Link.CatPSF},Extra2,InPar.MaxGet); %,BaseURL);
                %%%%FileName.Files_Cat=system(sprintf('wget --load-cookies=ptf.txt "http://irsa.ipac.caltech.edu/ibe/data/ptf/images/level1/Link(i).Cat"'));

             otherwise
                 % do nothing
         end
    end
    %A@anohma;
end
       
