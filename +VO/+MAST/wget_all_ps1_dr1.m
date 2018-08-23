function wget_all_ps1_dr1(CreateTile,LastDigit,MinNeededTiles)
% Prepare a local copy of the PS1-DR1 catalog
% Package: VO.MAST
% Example: VO.MAST.wget_all_ps1_dr1(false,0,300000)

DataDir = '/raid/eran/Catalogue/PS1/';

DataDir = '/raid/eran/Catalogue/PS1/missing';

cd(DataDir);

RAD = 180./pi;
StepSize = 0.2;

if (nargin<3)
    MinNeededTiles = 0;
end

if (CreateTile)
    [TileList,TileArea]=celestial.coo.tile_the_sky(360.*5,180./StepSize);
    % all tiles in TileList start at RA=0 (no RA=0 crossing...)
    TileList = TileList.*RAD;  % convert to deg
    % remove tiles below dec of -32
    TileList = TileList(TileList(:,2)>-32,:);
    save TileListPS1.mat TileList
else
    load('TileListPS1.mat');
end

%%
% Query{1} = {'raMean', 'decMean',...
%                     'raMeanErr', 'decMeanErr',...
%                     'raStack', 'decStack',...
%                     'epochMean', 'posMeanChisq', 'nDetections',... 
%                     'gPSFMag','gPSFMagErr','gpsfLikelihood',...
%                     'rPSFMag','rPSFMagErr','rpsfLikelihood',...
%                     'iPSFMag','iPSFMagErr','ipsfLikelihood',...
%                     'zPSFMag','zPSFMagErr','zpsfLikelihood',...
%                     'yPSFMag','yPSFMagErr','ypsfLikelihood'};
% 
% Query{2} = {'StackObjectView'};



% TableJoin = '';



Query{1} = {'StackObjectView.raStack', 'StackObjectView.decStack',...
                    'StackObjectView.raStackErr','StackObjectView.decStackErr',...
                    'StackObjectView.epochMean', 'StackObjectView.posMeanChisq',...
                    'gPSFMag','gPSFMagErr','gpsfLikelihood','gMeanPSFMagStd','gMeanPSFMagNpt','gMeanPSFMagMin','gMeanPSFMagMax',...
                    'rPSFMag','rPSFMagErr','rpsfLikelihood','rMeanPSFMagStd','rMeanPSFMagNpt','rMeanPSFMagMin','rMeanPSFMagMax',...
                    'iPSFMag','iPSFMagErr','ipsfLikelihood','iMeanPSFMagStd','iMeanPSFMagNpt','iMeanPSFMagMin','iMeanPSFMagMax',...
                    'zPSFMag','zPSFMagErr','zpsfLikelihood','zMeanPSFMagStd','zMeanPSFMagNpt','zMeanPSFMagMin','zMeanPSFMagMax',...
                    'yPSFMag','yPSFMagErr','ypsfLikelihood','yMeanPSFMagStd','yMeanPSFMagNpt','yMeanPSFMagMin','yMeanPSFMagMax'};

Query{2} = 'FROM StackObjectView left join MeanObjectView on StackObjectView.objID=MeanObjectView.objID';
        
TableJoin = 'StackObjectView.objID=MeanObjectView.objID and';
Query{3}  = 'StackObjectView.objID=MeanObjectView.objID';


%%
QueryBack{1} = {'raStack', 'decStack',...
                    'raStackErr','decStackErr',...
                    'epochMean', 'posMeanChisq',...
                    'gMeanPSFMag','gMeanPSFMagErr','null','gMeanPSFMagStd','gMeanPSFMagNpt','gMeanPSFMagMin','gMeanPSFMagMax',...
                    'rMeanPSFMag','rMeanPSFMagErr','null','rMeanPSFMagStd','rMeanPSFMagNpt','rMeanPSFMagMin','rMeanPSFMagMax',...
                    'iMeanPSFMag','iMeanPSFMagErr','null','iMeanPSFMagStd','iMeanPSFMagNpt','iMeanPSFMagMin','iMeanPSFMagMax',...
                    'zMeanPSFMag','zMeanPSFMagErr','null','zMeanPSFMagStd','zMeanPSFMagNpt','zMeanPSFMagMin','zMeanPSFMagMax',...
                    'yMeanPSFMag','yMeanPSFMagErr','null','yMeanPSFMagStd','yMeanPSFMagNpt','yMeanPSFMagMin','yMeanPSFMagMax'};

QueryBack{2} = 'FROM MeanObjectView';
        
TableJoin = '';
QueryBack{3}  = '';

Backup = true;  % use only if Stack failed
if (Backup)
    Query = QueryBack;
end


%Q{1}='SELECT o.raStack, o.decStack, o.raMean, o.decMean, o.raMeanErr, o.decMeanErr, o.epochMean, o.posMeanChisq, o.nDetections, m.gMeanPSFMag,m.gMeanPSFMagErr,m.gMeanPSFMagStd,m.gMeanPSFMagNpt,m.gMeanPSFMagMin,m.gMeanPSFMagMax,m.gFlags,s.gpsfLikelihood, m.rMeanPSFMag,m.rMeanPSFMagErr,m.rMeanPSFMagStd,m.rMeanPSFMagNpt,m.rMeanPSFMagMin,m.rMeanPSFMagMax,m.rFlags,s.rpsfLikelihood, m.iMeanPSFMag,m.iMeanPSFMagErr,m.iMeanPSFMagStd,m.iMeanPSFMagNpt,m.iMeanPSFMagMin,m.iMeanPSFMagMax,m.iFlags,s.ipsfLikelihood, m.zMeanPSFMag,m.zMeanPSFMagErr,m.zMeanPSFMagStd,m.zMeanPSFMagNpt,m.zMeanPSFMagMin,m.zMeanPSFMagMax,m.zFlags,s.zpsfLikelihood, m.yMeanPSFMag,m.yMeanPSFMagErr,m.yMeanPSFMagStd,m.yMeanPSFMagNpt,m.yMeanPSFMagMin,m.yMeanPSFMagMax,m.yFlags,s.ypsfLikelihood';
%Q{2}='FROM objectThin as o inner join MeanObject as m on o.objid=m.objid inner join StackObjectAttributes as s on o.objid=s.objid';
%Q{3}='o.decMean between 0 and 0.5 and o.raMean between 0 and 0.5';




Ntile = size(TileList,1);

%%

StopLoop = false;
while ~StopLoop
%for Itile=I1:1:I2,
    tic;
    % search for tile to work on:
    
    FailLoad = true;
    while FailLoad
        try
            SF = load('SummaryMaster.txt');
            FailLoad = false;
        catch
            fprintf('Failed loading SummaryMaster.txt - wait 3s\n');
            pause(3);
        end
    end
    fprintf('A: %f\n',toc);
    
    SF = sortrows(SF,1);
    % search for lowest tile index that wasn't query yet
    NeededTiles = setdiff((1:1:Ntile)',SF(:,1));
        
    II    = find((NeededTiles-LastDigit)./10==floor((NeededTiles-LastDigit)./10) & NeededTiles>MinNeededTiles,1);
    Itile = NeededTiles(II);
    if (isempty(Itile))
        error('--- End ---');
    end
    fprintf('B: %f\n',toc);
    Itile
    FileName = sprintf('PS_Tile_%06d.hdf5',Itile);
    
    if (~exist(FileName,'file'))
        % if file already exist then skip
        
        MinRA  = TileList(Itile,3);
        MaxRA  = TileList(Itile,4);
        MinDec = TileList(Itile,5);
        MaxDec = TileList(Itile,6);

        
        
        if (~Backup)
            tic;
            [Out,~,Status,ResultOrig,Query]=VO.MAST.query_casjobs_recur(Query,...
                            'Table','PanSTARRS_DR1',...
                            'boxcoo',[MinRA MaxRA MinDec MaxDec].*pi./180,...
                            'StrRA','StackObjectView.raMean',...
                            'StrDec','StackObjectView.decMean');
            TT=toc;
        else
            tic;
            [Out,~,Status,ResultOrig,Query]=VO.MAST.query_casjobs_recur(Query,...
                            'Table','PanSTARRS_DR1',...
                            'boxcoo',[MinRA MaxRA MinDec MaxDec].*pi./180,...
                            'StrRA','raMean',...
                            'StrDec','decMean');
            TT=toc;
        end
        
        fprintf('1: %f\n',toc);
        
        if (Status<0)
            error('Prolem unrecovered');
        end
        
        if (isempty(Out))
            Nsrc = 0;
        else
            Nsrc = size(Out.Cat,1);
        end
        
        fprintf('-----------------------------------------------\n');
        fprintf('   Completed tile number: %d\n',Itile);
        fprintf('   Run time: %f\n',TT);
        fprintf('   Number of sources: %d\n',Nsrc);
        
        if (~isempty(Out))
            Nout = size(Out.Cat,1);
            if (size(Out.Cat,1)>0)
                cd(DataDir);

                % save the files in a different directory to avoid large
                % number of files in the search path
                cd ./PS_Tile/
                Util.IO.saveh(FileName,Out.Cat);
                cd ..
            end
        else
            Nout = 0;
        end

        fprintf('2: %f\n',toc);
        
        % print summary
        FID = fopen('SummaryMaster.txt','a');
        fprintf(FID,'%7d  %6d\n',Itile,Nout);
        fclose(FID);
        
        fprintf('3: %f\n',toc);
        % force quit
        %pause(1);
        %if (exist('stop_wget_ps1_dr1','file')>0)
        if (java.io.File('stop_wget_ps1_dr1').exists)
            StopLoop = true;
            %delete('stop_wget_ps1_dr1');
        end
        fprintf('4: %f\n',toc);
    end   
end

                    
%         
%         QueryLine = Query;
%         %QueryLine{3} = sprintf('%s StackObjectView.decMean between %9.5f and %9.5f and StackObjectView.raMean between %9.5f and %9.5f',TableJoin,MinRA,MaxRA,MinDec,MaxDec);
%         % faster:
%         QueryLine{3} = sprintf('%s StackObjectView.decMean>%9.5f and StackObjectView.decMean<=%9.5f and StackObjectView.raMean>%9.5f and StackObjectView.raMean<=%9.5f',TableJoin,MinDec,MaxDec,MinRA,MaxRA);
%         %
% 
%         tic;
%         [Out,ColCell,Status,Result] = Catalog.MAST.query_casjobs(QueryLine,'Table','PanSTARRS_DR1');
%         T=toc;
%         
%         Out
%         [Itile, T, Status]
%         
%         
%         if (Status<0 && Status>-101)
%             fprintf('Problem 1 - wait 100s and try again\n');
%             Result
%             pause(100);
%         else
% 
%             if (Status==-101)
%                 % query is too big - break to 4x4
%                 Nsize = 7;
%                 K = 0;
%                 B = zeros(Nsize.*Nsize,4);
%                 for Isize=0:1:Nsize-1
%                     for Jsize=0:1:Nsize-1
%                         K = K + 1;
%                         B(K,:) = [Isize, Isize+1, Jsize, Jsize+1];
%                     end
%                 end
% 
%                 %B=[0 1 0 1;0 1 1 2; 0 1 2 3; 0 1 3 4;1 2 0 1; 1 2 1 2;1 2 2 3;1 2 3 4;2 3 0 1;2 3 1 2;2 3 2 3;2 3 3 4;3 4 0 1;3 4 1 2;3 4 2 3; 3 4 3 4];
%                 ListEdge = B.*StepSize./Nsize;
%                 %[ListEdge,ListCenter]=ImUtil.Im.image_blocks([201 201],[-4 -4]);
%                 %ListEdge = (ListEdge-1)./1000;
% 
%                 Nedge = size(ListEdge,1);
%                 clear OutE;
%                 for Iedge=1:1:Nedge
% 
%                     MinRA1   = MinRA  + ListEdge(Iedge,1);
%                     MaxRA1   = MinRA  + ListEdge(Iedge,2);
%                     MinDec1  = MinDec + ListEdge(Iedge,3);
%                     MaxDec1  = MinDec + ListEdge(Iedge,4);
% 
%                     %QueryLine{3} = sprintf('%s raMean>%9.5f and raMean<=%9.5f and decMean>%9.5f and decMean<=%9.5f',TableJoin,MinRA1,MaxRA1,MinDec1,MaxDec1);
%                     QueryLine{3} = sprintf('%s StackObjectView.decMean>%9.5f and StackObjectView.decMean<=%9.5f and StackObjectView.raMean>%9.5f and StackObjectView.raMean<=%9.5f',TableJoin,MinDec1,MaxDec1,MinRA1,MaxRA1);
% 
%                     tic;
%                     [OutE(Iedge),ColCell,Status,Result] = Catalog.MAST.query_casjobs(QueryLine,'Table','PanSTARRS_DR1');
%                     T=toc;
% 
%                     OutE(Iedge)
%                     [Itile, T, Status]
% 
%                     if (Status<0),
%                         [Itile, Iedge]
%                         Result
%                         error('Problem 2');
%                     end
% 
%                     Out = merge(OutE);
%                 end
%             end
%         
        
         
        %Report(Itile).Status = Status;
        %Report(Itile).time   = T;

