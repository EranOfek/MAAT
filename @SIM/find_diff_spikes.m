function [Sim,SpikesDB,SimFilt,SimFiltIndiv]=find_diff_spikes(Sim,varargin)
% Find diffraction spikes in image and set the bit mask image.
% Description: Find diffraction spikes in image and set the bit mask image.
%              Return the Spikes database (position and arm lengths), and
%              populated the SIM object mask.
% Input  : - A SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ColX'  - Cell array of column names that may contain the
%                      X-axis position of the sources in the image in the
%                      associated AstCat object. Will select the first
%                      existing column name.
%                      Default is {'XWIN_IMAGE','X_IMAGE','X'}.
%            'ColY'  - Like ColX but for the Y-axis.
%                      Default is {'YWIN_IMAGE','Y_IMAGE','Y'}.
%            'SpikesTheta' - The angles of the spikes measured relative to
%                      the X-axis [deg] in the range of 0 to 180 deg.
%                      Default is [0 90].
%            'SpikesFilterLength' - The spike arm length in the spike
%                      filter. Default is 31.
%            'SpikesFilterGap' - The length of the central gap in the
%                      spike. Default is 15.
%            'SpikesFilterWidth' - Spike width. Default is 1.
%            'MaxSpikeLength' - Maximum spike length. Default is 1000.
%            'SpikeLengthThresh' - Detection threshold (sigmas) for
%                      measuring spike length. Default is 2.
%            'SpikeFiltThresh' - Spike filter detection threshold (sigmas).
%                      Default is 8.
%            'SpikeFiltIndivThresh' - Spike single arm filter detection
%                      threshold (sigmas). Default is 5.
%            'StarMinFlux' - Minimum count at location of spike producing
%                      source. Default is 50000.
%            'Plot'  - Plot spikes locations in ds9. ds9 should be open.
%                      Default is false.
%            'BackPar' - Parameters to pass to SIM/background.
%                      Default is {}.
%            'SubBack' - Subtract background from image. Default is true.
%            'SrcExtractionProg' - Source extraction program.
%                      Default is @mextractor.
%            'SrcExtractionPar' - Parameters to pass to the source
%                      extraction program. Default is {}.
%            'BitName' - Bit name in bitmask dictionary for the diffraction
%                      spikes flagging. Default is 'Bit_Spike'.
%            'BitType' - Bit mask type (if not defined).
%                      Default is 'uint32'.
%            'CombFun' - Bit mask combining function. Default is @bitor.
%            'SpikeGap'- The distance from the spikes generated source that
%                      will not be flagged as spike region (as this is the
%                      spike creator). Default is 4 pix.
%            'MaxSpikeAssym'- The maximum allowed assymetry between spike
%                      length. If ratio of spikes length is larger than
%                      this factor than set the longest arm to this factor
%                      multiplies by the shortest spike length.
%                      Default is 2.
% Output : - A SIM object with the updated mask image.
%          - A spikes database.
%          - Spikes filtered image (single).
%          - Individual spikes filtered images.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Sep 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Sim,SpikesDB]=find_diff_spikes(Sim)
% Reliable: 2
%--------------------------------------------------------------------------

ImageField   = SIM.ImageField;
BackField    = SIM.BackField;
CatField     = AstCat.CatField;

DefV.ColX                  = {'XWIN_IMAGE','X_IMAGE','X'};
DefV.ColY                  = {'YWIN_IMAGE','Y_IMAGE','Y'};
DefV.SpikesTheta           = [0,90];    % 0 is X-axis
DefV.SpikesFilterLength    = 31;
DefV.SpikesFilterGap       = 15;
DefV.SpikesFilterWidth     = 1;
DefV.MaxSpikeLength        = 1000;
DefV.SpikeLengthThresh     = 2;
DefV.SpikeFiltThresh       = 8;
DefV.SpikeFiltIndivThresh  = 5;
DefV.StarMinFlux           = 50000;
DefV.Plot                  = false;  % plot in ds9
DefV.BackPar               = {};
DefV.SubBack               = true;    % subtract background
DefV.SrcExtractionProg     = @mextractor;
DefV.SrcExtractionPar      = {};
DefV.BitName               = 'Bit_Spike';
DefV.BitType               = 'uint32';
DefV.CombFun               = @bitor;
DefV.SpikeGap              = 4;  
DefV.MaxSpikeAssym         = 2;

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

SpikesArmTheta = [InPar.SpikesTheta(:); InPar.SpikesTheta(:)+180];
Narm           = numel(SpikesArmTheta);

ColX = select_exist_colnames(Sim,InPar.ColX(:));
ColY = select_exist_colnames(Sim,InPar.ColY(:));
ColX = ColX{1};
ColY = ColY{1};

% Check if catalog is populated
if (~isfield_populated(Sim,CatField))
    Sim = InPar.SrcExtractionProg(Sim,InPar.SrcExtractionPar{:});
end


Ntheta = numel(InPar.SpikesTheta);

SpikesFilter = zeros(InPar.SpikesFilterLength,InPar.SpikesFilterLength);
for Itheta=1:1:Ntheta
    Ker(Itheta).Filt=Kernel2.line(InPar.SpikesFilterLength,InPar.SpikesFilterWidth,InPar.SpikesTheta(Itheta),InPar.SpikesFilterGap); 
    % sum of spikes filters
    SpikesFilter = SpikesFilter + Ker(Itheta).Filt;
end
% Normalize SpikesFilter
SpikesFilter = SpikesFilter./sum(SpikesFilter(:));

if (~isfield_populated(Sim,BackField))
    Sim    = background(Sim,InPar.BackPar{:});
end

if (InPar.SubBack)
    SimB    = sub_background(Sim);
else
    SimB    = Sim;
end

Nsim = numel(Sim);
SpikesDB = Util.struct.struct_def({'Src'},Nsim,1);

for Isim=1:1:Nsim
    
    X = col_get(Sim(Isim),ColX);
    Y = col_get(Sim(Isim),ColY);

    SizeIm = imagesize(Sim(Isim));
    % Allocate the FlagImage of diffraction spikes
    FlagImage = false(fliplr(SizeIm));
    
    X = min(X,SizeIm(1));
    X = max(X,1);
    Y = min(Y,SizeIm(2));
    Y = max(Y,1);
    Ind = Util.array.sub2ind_fast(fliplr(SizeIm),round(Y),round(X));
    Flag = true(size(Ind));
    
    % Filter and normalize the image by std
    SimFilt = filter(SimB,SpikesFilter,'NormErr',true);
    SimFiltIndiv = SIM(Ntheta);
    for Itheta=1:1:Ntheta
        % Filter and normalize the image by std
        SimFiltIndiv(Itheta) = filter(SimB,Ker(Itheta).Filt,'NormErr',true);
        Flag = Flag & SimFiltIndiv(Itheta).(ImageField)(Ind)>InPar.SpikeFiltIndivThresh;
    end
    
    % Search for spike sources candidates
    Flag = Flag & SimFilt.(ImageField)(Ind)>InPar.SpikeFiltThresh & Sim.(ImageField)(Ind)>InPar.StarMinFlux; 

    %ds9.plot(X(Flag),Y(Flag),'rs')

    % go over all spikes
    Xs = X(Flag);
    Ys = Y(Flag);
    Ns = numel(Xs);
    SpikesDB(Isim).Src = Util.struct.struct_def({'XI','YI','X','Y','Length','Theta'},Ns,1);

    for Is=1:1:Ns
        XI = round(Xs(Is));
        YI = round(Ys(Is));

        %ArmLength = zeros(Narm,1);
        for Iarm=1:1:Narm
            Itheta = Iarm - Ntheta.*floor(Iarm./Ntheta - eps);

            SpikeRadVec = (1:1:InPar.MaxSpikeLength)';
            Xv{Iarm}           = XI + SpikeRadVec.*cosd(SpikesArmTheta(Iarm));
            Yv{Iarm}           = YI + SpikeRadVec.*sind(SpikesArmTheta(Iarm));
            % make sure X and Y are within image boundries
            Xv{Iarm} = max(Xv{Iarm},1);
            Xv{Iarm} = min(Xv{Iarm},SizeIm(1));
            Yv{Iarm} = max(Yv{Iarm},1);
            Yv{Iarm} = min(Yv{Iarm},SizeIm(2));
            
            MaxLength = find((abs(diff(Xv{Iarm}))~=0 | abs(diff(Yv{Iarm}))~=0)==0,1,'first');

            Ind = Util.array.sub2ind_fast(fliplr(SizeIm),round(Yv{Iarm}),round(Xv{Iarm}));
            Arm = SimFiltIndiv(Itheta).(ImageField)(Ind);
            %ArmLength(Iarm) = find(Arm<InPar.SpikeThreshold,1,'first');
            
            SpikesDB(Isim).Src(Is).XI = XI;   % rounded source position
            SpikesDB(Isim).Src(Is).YI = YI;
            SpikesDB(Isim).Src(Is).X  = Xs(Is);  % original source position
            SpikesDB(Isim).Src(Is).Y  = Ys(Is);
            
          
            % spike length
            Length = find(Arm<InPar.SpikeLengthThresh,1,'first');
            if (isempty(Length))
                Length = numel(Arm);
            end
            if (isempty(Length))
                % assume object is near edge
                % set to distance from edge
                Length = MaxLength;
            else
                Length = Length - 0.5.*InPar.SpikesFilterLength;
                Length = max(Length,0.5.*InPar.SpikesFilterLength);
            end
            
            SpikesDB(Isim).Src(Is).Length(Iarm) = Length;
            SpikesDB(Isim).Src(Is).Theta(Iarm)  = SpikesArmTheta(Iarm);  % [deg] from X axis
        end
        
        % set the length of spikes which are much longer (factor
        % InPar.MaxSpikeAssym) than the shortest spike arm length
        % to InPar.MaxSpikeAssym times the shortest arm length
        ShortestSpikeLength = min(SpikesDB(Isim).Src(Is).Length);
        SpikeLengthRatio    = SpikesDB(Isim).Src(Is).Length./ShortestSpikeLength;
        FlagLargeRatio      = SpikeLengthRatio>InPar.MaxSpikeAssym;
        SpikesDB(Isim).Src(Is).Length(FlagLargeRatio) = ShortestSpikeLength.*InPar.MaxSpikeAssym;
        
        
        for Iarm=1:1:Narm    
            % update the Flag image for the MASK object
            Length = SpikesDB(Isim).Src(Is).Length(Iarm);
           
            IndF = Util.array.sub2ind_fast(fliplr(SizeIm),round(Yv{Iarm}(InPar.SpikeGap:floor(Length))),round(Xv{Iarm}(InPar.SpikeGap:floor(Length))));
            FlagImage(IndF) = true;
            
            
            if (InPar.Plot)
                ds9.line_lt(SpikesDB(Isim).Src(Is).XI,SpikesDB(Isim).Src(Is).YI,SpikesDB(Isim).Src(Is).Length(Iarm),SpikesDB(Isim).Src(Is).Theta(Iarm));
            end
        end
    end
    % Update the MASK object
    Sim(Isim) = bitmask_set(Sim(Isim),FlagImage,InPar.BitName,InPar.BitType,InPar.CombFun);

end

