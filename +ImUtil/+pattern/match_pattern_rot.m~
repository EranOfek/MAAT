function Res=match_pattern_rot(Cat,Ref,varargin)
% Match two catalogs with unknown shift and rotation and find rotations
% Package: ImUtil.pattern
% Description: Given two catalogs with the same scale, but unknown shift,
%              rotation and axes flips. Find all the possible rotation and
%              flips that may result in matching the catalogs.
% Input  : - A matrix containing a catalog (e.g., [X, Y]).
%            If only one input argument is provided then work in simulation
%            mode and this parameter is the simulated rotation angle in
%            deg.
%          - A matrix containing a reference catalog (e.g., [X, Y]).
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Thresh'  - The minimum threshold in units of sigmas for
%                        selecting possible rotations that may match the
%                        two catalogs. Default is 5.
%            'MinThreshRatio' - Select only peaks which height is larger
%                        than the maximum peak hight multiple by 
%                        MinThreshRatio. Default is 0.5.
%            'Duplicate180'  - Duplicate solutions to solutions with
%                        +180 deg. Default is true. 
%            'Flip'    - Matrix of [X,Y] flips to try use in the matching
%                        process. Default is [1 1; 1 -1].
%                        This default with Duplicate180=true covers all
%                        possible flips.
%            'CatColX' - Column index for X axis in the Catalog.
%                        Default is 1.
%            'CatColY' - Column index for Y axis in the Catalog.
%                        Default is 2.
%            'RefColX' - Column index for X axis in the Reference.
%                        Default is 1.
%            'RefColY' - Column index for Y axis in the Reference.
%                        Default is 2. 
%            'HistDistEdges' - Vector of distamces to test.
%                        Default is (12:3:600)'.
%            'HistRotEdges'  - Vector of angular rotations to test.
%                        Default is (-90:0.2:90)'.
% Output : - A structure with the following fields:
%            'MatchedRot' - A matrix of all possible rotation and flips
%                        that may match the catalogs
%                        [Rot, Highet (sigma), Xflip, Yflip].
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Mar 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res=ImUtil.pattern.match_pattern_rot(56);
% Reliable: 2
%--------------------------------------------------------------------------

RAD = 180./pi;

if (nargin<2)
    if (nargin==0)
        Theta = 0;
    else
        Theta = Cat;
    end
    % simulation mode
    Ns = 5000;
    Nm = 2000;
    
    Ref = rand(Ns,2).*2048;
    Cat = Ref(Nm,:) + rand(Nm,2).*0.3;
    Cat = [Cat; rand(Nm,2).*2048];
    
    Ref = Ref.*[1 1];
    
    % Rotate Ref
    %Theta = -10;
    Rot = [cosd(Theta), -sind(Theta); sind(Theta), cosd(Theta)];
    
    Cat = (Rot*Cat.').';
end

DefV.Thresh               = 5;   % [sigma]
DefV.MinThreshRatio       = 0.5; % select only peaks > MaxPeak*MinThreshRatio
DefV.Duplicate180         = true;
DefV.Flip                 = [1 1; 1 -1;-1 -1; -1 1];
DefV.CatColX              = 1;
DefV.CatColY              = 2;
DefV.RefColX              = 1;
DefV.RefColY              = 2;
DefV.FitScale             = false;
DefV.HistDistEdges        = (12:3:300)';
DefV.HistRotEdges         = (-90:0.2:90)';
DefV.CutRefCat            = true;
DefV.SearchRangeX         = [-1000 1000];
DefV.SearchRangeY         = [-1000 1000];


InPar = InArg.populate_keyval(DefV,varargin,mfilename);

MaxDist = max(InPar.HistDistEdges);

ThetaShift = InPar.HistRotEdges - min(InPar.HistRotEdges);

Nflip = size(InPar.Flip,1);
AllListMax = zeros(0,4);
for Iflip=1:1:Nflip

    % X/Y columns
    Xc = Cat(:,InPar.CatColX);
    Yc = Cat(:,InPar.CatColY);
    Xr = Ref(:,InPar.RefColX);
    Yr = Ref(:,InPar.RefColY);

    Xc = Xc.*InPar.Flip(Iflip,1);
    Yc = Yc.*InPar.Flip(Iflip,2);

    %InPar.CutRefCat=false;

    if (InPar.CutRefCat)
        MinX = min(Xc) - max(abs(InPar.SearchRangeX));
        MaxX = max(Xc) + max(abs(InPar.SearchRangeX));
        MinY = min(Yc) - max(abs(InPar.SearchRangeY));
        MaxY = max(Yc) + max(abs(InPar.SearchRangeY));

        % use only Ref sources in SubCat region
        FinCut = Xr>MinX & ...
                 Xr<MaxX & ...
                 Yr>MinY & ...
                 Yr<MaxY;

        Xr = Xr(FinCut);
        Yr = Yr(FinCut);
    
    end

    
    % all possible differences between sources in Cat
    DXc = Xc - Xc.';
    DYc = Yc - Yc.';
    % all possible differences between sources in Ref
    DXr = Xr - Xr.';
    DYr = Yr - Yr.';
    
    Fc = abs(DXc(:))<MaxDist & abs(DYc(:))<MaxDist;
    Fr = abs(DXr(:))<MaxDist & abs(DYr(:))<MaxDist;
    DXc = DXc(Fc);
    DYc = DYc(Fc);
    DXr = DXr(Fr);
    DYr = DYr(Fr);
    
    
    % all possible distances/angle between sources in Cat
    Dc  = sqrt(DXc.^2 + DYc.^2);
    Tc  = atan(DYc./DXc);
    % all possible distances/angle between sources in Cat
    Dr  = sqrt(DXr.^2 + DYr.^2);
    Tr  = atan(DYr./DXr);
    

    % calc histograms
    Nc=histcounts2((Dc(:)),Tc(:).*RAD,InPar.HistDistEdges,InPar.HistRotEdges);
    Nr=histcounts2((Dr(:)),Tr(:).*RAD,InPar.HistDistEdges,InPar.HistRotEdges);


    %surface(Nc); shading interp; colorbar

    % subtract background
    Nc = Nc - mean(Nc(:));
    Nr = Nr - mean(Nr(:));
    %CC = (ifft2(fft2(Nr).*conj(fft2(Nc))));
    % cross correlated
    CC = (ifft2(fft2((Nc)).*conj(fft2((Nr)))));
    CC1 = CC(1,:);
    CC1 = (CC1 - mean(CC1))./std(CC1);

    % [ListExt] = Util.find.find_local_extramum(InPar.HistRotEdges(:),CC1(:))
    % I = find(ListExt(:,3)<0 & ListExt(:,2)>InPar.Thresh)
    % ListExt(I,:)

    Fmax = CC1>InPar.Thresh;
    CC1F = CC1(Fmax);
    ListMax = [ThetaShift(Fmax), CC1F(:)];

    if (~isempty(ListMax))
        if (InPar.Duplicate180)
            ListMax = [ListMax; ListMax(:,1)+180, ListMax(:,2)];
        end

        Fh = ListMax(:,2)>(max(ListMax(:,2)).*InPar.MinThreshRatio);
        ListMax = ListMax(Fh,:);
    else
        ListMax = zeros(0,2);
    end
    
    Nl = size(ListMax,1);
    AllListMax = [AllListMax; [ListMax, InPar.Flip(Iflip,:).*ones(Nl,1) ]];

end
Res.MatchedRot = AllListMax;