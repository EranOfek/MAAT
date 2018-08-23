function wget_all_hsc(varargin)
% SHORT DESCRIPTION HERE
% Package: VO
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Apr 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: VO.prep.wget_all_hsc('DecRange',DecRange,'Problem',Problem,'UnDetector',UnDetector);
% Reliable: 
%--------------------------------------------------------------------------

RAD = 180./pi;
CatField = AstCat.CatField;

%DefV.SaveInd              = 0;
DefV.MaxSrcQ              = 5000;
DefV.DecRange             = zeros(0,3);
DefV.Problem              = [];
DefV.UnDetector           = {};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


%%

HSTimages = cats.sources.HSTimages;
Nim = size(HSTimages.(CatField),1);



%SaveInd  = InPar.SaveInd;
DecRange = InPar.DecRange;
Problem  = InPar.Problem;
UnDetector = InPar.UnDetector;
K = 0;
Problem = [];
for Iim=37501:1:Nim
    Iim
    
    Nsc = HSTimages.(CatField)(Iim,4);
    Nq  = ceil(Nsc./InPar.MaxSrcQ);
    if (Nq>1)
        [Nq, Nsc]
    end
    CONT = true;
    Iq = 0;
    while CONT
        Iq = Iq + 1;
        %for Iq=1:1:Nq
        StartLine = 1+(Iq-1).*InPar.MaxSrcQ;
        EndLine   = StartLine + InPar.MaxSrcQ - 1;
    
        Q1 = 'SELECT MatchRA, MatchDec, SourceRA, SourceDec, MatchID, MemID, SourceID, ImageID, Detector, WaveLength, ExposureTime, StartMJD, MagAper2, MagAuto, KronRadius';
        %Q2 = 'FROM DetailedCatalog';
        Q2 = sprintf('FROM (SELECT ROW_NUMBER() OVER(ORDER BY(SELECT NULL as noorder)) As RowNum, * FROM  DetailedCatalog WHERE ImageID=%d) as alias',HSTimages.(CatField)(Iim,3));
        
        Format = '%f %f %f %f %f %f %f %f %s %f %f %f %f %f %f\n';

        RA  = HSTimages.(CatField)(Iim,1);
        Dec = HSTimages.(CatField)(Iim,2);
        [RA1,RA2,Dec1,Dec2] = celestial.coo.coo2box(RA,Dec,3./60./RAD);
        %Q3 = sprintf('WHERE ImageID=%d and RowNum between 1 and 20',HSTimages.(CatField)(Iim,3));
        Q3 = sprintf('WHERE RowNum between %d and %d',StartLine,EndLine);
        QS = sprintf('%s %s %s',Q1,Q2,Q3);
        pause(0.1);
        [C,~,~,Out] = VO.MAST.query_casjobs(QS,'Table','HSCv2','SaveInTable',true,'FormatString',Format);
        if (isempty(C))
            Problem = [Problem; Iim];
            warning('Problem')
            save Problem.mat Problem
            %K = K + 1;
            CONT = false;
            pause(100);
        else
            
            Nsrc = size(C.(CatField),1);
            Nsrc
            if (Nsrc<InPar.MaxSrcQ)
                CONT = false;

            end
            %C.(CatField).Instrument = regexprep(C.(CatField).Instrument,'"','');
            C.(CatField).Detector = regexprep(C.(CatField).Detector,'"','');
            %C.(CatField).Filter = regexprep(C.(CatField).Filter,'"','');


            %if (Iim==1 && Iq==1)
            if (isempty(UnDetector))
                %UnInstrument = unique(C.(CatField).Instrument);
                UnDetector = unique(C.(CatField).Detector);
            else
                %UnInstrument = unique([UnInstrument; C.(CatField).Instrument],'stable');
                UnDetector   = unique([UnDetector;   C.(CatField).Detector],'stable');
            end

            Ndet = numel(UnDetector);
            Detector = nan(Nsrc,1);
            for Idet=1:1:Ndet
                Detector(strcmp(C.(CatField).Detector,UnDetector{Idet})) = Idet;
            end
            C.(CatField).Detector = Detector;

            C.(CatField) = table2array(C.(CatField));

            K = K + 1;
            AllC(K) = C;
        end
    end
   
    if (Iim./100==floor(Iim./100))
        %SaveInd = SaveInd + 1;
        
        MC = merge(AllC);
        MC = sortrows(MC,2);
        MinDec = MC.Cat(1,2);
        MaxDec = MC.Cat(end,2);
        
        FileName = sprintf('hsc_%d.mat',Iim);
        save('-v7.3',FileName,'MC');
        DecRange = [DecRange; [Iim, MinDec, MaxDec]];
        save DecRange.mat DecRange UnDetector
        
        clear AllC;
        K = 0;
    end
    
end


        MC = merge(AllC);
        MC = sortrows(MC,2);
        MinDec = MC.Cat(1,2);
        MaxDec = MC.Cat(end,2);
        
        FileName = sprintf('hsc_%d.mat',Iim);
        save('-v7.3',FileName,'MC');
        DecRange = [DecRange; [Iim, MinDec, MaxDec]];
        save DecRange.mat DecRange UnDetector
        
        clear AllC;
        K = 0;


save UnDetector.mat UnDetector
 
%%
DecRange = zeros(0,3);

Files = dir('hsc_*.mat');
Nf = numel(Files);
for If=1:1:Nf
    If
    A=regexp(Files(If).name,'hsc_(?<N>\d+).mat','names');
    
    Cat = Util.IO.load2(Files(If).name);
    DecRange = [DecRange; [str2double(A.N), min(Cat.Cat(:,2)), max(Cat.Cat(:,2))]];
    
end
save DecRange.mat DecRange

%%

ColMin = 2;
ColMax = 3;

DecVec = (-90:1:90)';
Nd = numel(DecVec);
for Id=1:1:Nd-1
    Id
    D1 = DecVec(Id);
    D2 = DecVec(Id+1);
    
    Flag = (DecRange(:,ColMax)>D1 & DecRange(:,ColMax)<D2) | ...
           (DecRange(:,ColMin)>D1 & DecRange(:,ColMin)<D2) | ...
           (DecRange(:,ColMin)<D1 & DecRange(:,ColMax)>D2);
       
    Iflag = find(Flag);
    Nf = numel(Iflag);
    Cat = AstCat(Nf,1);
    for If=1:1:Nf
        Iim = DecRange(Iflag(If),1);
        File = sprintf('hsc_%d.mat',Iim);
        Cat(If) = Util.IO.load2(File);
    end
    Cat = merge(Cat);
    
    Cat.ColUnits   = {'rad','rad','rad','rad','','','','','','Ang','s','MJD','mag','mag','arcsec'};
    
    if (~isempty(Cat.Cat))
        Cat = sortrows(Cat,2);
        FlagD = Cat.Cat(:,2)>=D1 & Cat.Cat(:,2)<D2;

        Cat.Cat = Cat.Cat(FlagD,:);
        if (~isempty(Cat.Cat))
            Cat.Cat(:,1:4) = Cat.Cat(:,1:4)./RAD;

            if (D1<0)
                D1s = 'm';
            else
                D1s = '';
            end
            if (D2<0)
                D2s = 'm';
            else
                D2s = '';
            end
            FileDec = sprintf('HSC_%s%d_%s%d.mat',D1s,abs(D1),D2s,abs(D2));
            save('-v7.3',FileDec,'Cat');
        end
    end
end



