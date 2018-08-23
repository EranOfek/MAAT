function OutC=bin2d_cat(AstC,varargin)
%--------------------------------------------------------------------------
% bin2d_cat function                                         class/@AstCat
% Description: Bin an AstCat object by some 2-D coordinates
%              (e.g., X and Y coordinates), and calculate some mean
%              properties in each bin.
% Input  : - An AstCat object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ColBinX' - X-axis column by which to bin.
%                        Default is 'XWIN_IMAGE'.
%            'ColBinY' - Y-axis column by which to bin.
%                        Default is 'YWIN_IMAGE'.
%            'LimitX'  - The [Min, Max] X range in which to bin the data.
%                        If empty, then will attempt to read it from the
%                        SIM image size or the NAXIS header keywords.
%                        Default Min value is 0. Default is empty.
%            'LimitY'  - The [Min, Max] Y range in which to bin the data.
%                        If empty, then will attempt to read it from the
%                        SIM image size or the NAXIS header keywords.
%                        Default Min value is 0. Default is empty.
%            'NstepX'  - Number of bins in the X axis. Default is 10.
%            'NstepY'  - Number of bins in the Y axis. Default is 10.
%            'StepX'   - The X step size. If not empty then this overrides
%                        the NstepX parameter. Default is empty.
%            'StepY'   - The Y step size. If not empty then this overrides
%                        the NstepY parameter. Default is empty.
%            'OutCol'  - The default columns in the output catalog.
%                        This is a 3 column cell array - the columns are:
%                        the output column name, function to use in order
%                        to calculate the value, and the column name or
%                        index in the input catalog to use in the
%                        calculation.
%                        Function can be any valid function handle or
%                        'CenX'/'CenY' strings which will calculate the bin
%                        center. Default is:
%                           {'CenX',   'CenX',   [];...
%                            'CenY',   'CenY',   [];...
%                            'Nsrc',   @numel,  'XWIN_IMAGE';...
%                            'MeanX',  @mean,   'XWIN_IMAGE';...
%                            'MeanY',  @mean,   'YWIN_IMAGE'};
%            'AddOutCol' - Like 'OutCol', but here the user can add
%                          additional columns. Default is {}.
% Output : - An AstCat object of the binned data.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: OutC=bin2d_cat(S1,'StepX',1024,'StepY',1024)
% Reliable: 2
%--------------------------------------------------------------------------


ImageField      = 'Im';
CatField        = 'Cat';
ColCellField    = 'ColCell';

DefV.ColBinX            = 'XWIN_IMAGE';
DefV.ColBinY            = 'YWIN_IMAGE';
DefV.LimitX             = [];
DefV.LimitY             = [];
DefV.NstepX             = 10;
DefV.NstepY             = 10;
DefV.StepX              = [];
DefV.StepY              = [];
DefV.OutCol             = {'CenX',   'CenX',   [];...
                           'CenY',   'CenY',   [];...
                           'Nsrc',   @numel,  'XWIN_IMAGE';...
                           'MeanX',  @mean,   'XWIN_IMAGE';...
                           'MeanY',  @mean,   'YWIN_IMAGE'};
DefV.AddOutCol          = {}; %'MeanX2', @median,    'X2WIN_IMAGE'};
         
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

if ~(all(isfield_populated(AstC,ColCellField)))
    error('Some catalogs in AstCat are empty');
end


% sort input catalogs by Y axis:
AstC = sortrows(AstC,InPar.ColBinY);

% Number of lines in InPar.OutCol
InPar.OutCol = [InPar.OutCol; InPar.AddOutCol];
Noutcol = size(InPar.OutCol,1);

Ncat = numel(AstC);
OutC = AstCat(size(AstC));
for Icat=1:1:Ncat
    % for each AstCat object element
    
    OutCol = InPar.OutCol;
    % replace 3rd column in OutCol with column index:
    OutCol(:,3) = num2cell(colname2ind(AstC(Icat),OutCol(:,3)));
    
    
    % Get limits in which to bin
    if (isempty(InPar.LimitX) || isempty(InPar.LimitY))
        LimitX = [];
        LimitY = [];
        if (SIM.issim(AstC))
            % AstC is a SIM image attemp to use image size
            Size = size(AstC(Icat).(ImageField));
            LimitX = Size(2);
            LimitY = Size(1);
        else
            % attempt to read size from Header NAXIS keywords
            Naxis = naxis(AstC(Icat));
            if (isempty(Naxis))
                error('Can not find LimitX/LimitY from image size or NAXIS keywords');
            else
                LimitX = Naxis(1);
                LimitY = Naxis(2);
            end
        end
    else
        LimitX = InPar.LimitX;
    end
    
    if (numel(LimitX)==1)
        LimitX = [0 LimitX];
    end
    if (numel(LimitY)==1)
        LimitY = [0 LimitY];
    end
    
    % Get step size for X axis
    if (isempty(InPar.StepX))
        StepX = (LimitX(2)-LimitX(1))./InPar.NstepX;
    else
        StepX = InPar.StepX;
    end
    % Get step size for Y axis
    if (isempty(InPar.StepY))
        StepY = (LimitY(2)-LimitY(1))./InPar.NstepY;
    else
        StepY = InPar.StepY;
    end
    
    % X/Y columns content
    VecX      = col_get(AstC(Icat),InPar.ColBinX);
    VecY      = col_get(AstC(Icat),InPar.ColBinY);
    
    
    EdgeX  = (LimitX(1):StepX:LimitX(2));
    EdgeY  = (LimitY(1):StepY:LimitY(2));
    NbinX = numel(EdgeX) - 1;
    NbinY = numel(EdgeY) - 1;
    Ibin  = 0;
    for Iy=1:1:NbinY
        % sources are sorted by Y axis - so select range
        % lower index of Y range
        Iy1 = bin_sear(VecY,EdgeY(Iy));
        % upper index of Y range
        Iy2 = bin_sear(VecY,EdgeY(Iy+1));
        % sub catalog containing only the rows which Y value is in range:
        SubVecX = VecX(Iy1:Iy2);
        
        for Ix=1:1:NbinX
            Ibin = Ibin + 1;
            % search for X values in range within SubCat
            % IndFound contains the indices of the rows in Ix/Iy bin:
            IndFound = find(SubVecX>EdgeX(Ix) & SubVecX<=EdgeX(Ix+1)) + Iy1 - 1;
            % Sub cat containing the bin data:
            SubCat   = AstC(Icat).(CatField)(IndFound,:);
            % calculate statisics in Bin:
            % Number of entries
            
            for Ioutcol=1:1:Noutcol
                % for each requested output column
                
                if (isa(OutCol{Ioutcol,2},'function_handle'))
                    % The function handle to calculate in the bin
                    Fun = OutCol{Ioutcol,2};
                    OutC(Icat).(CatField)(Ibin,Ioutcol) = Fun(SubCat(:,OutCol{Ioutcol,3}));
                else
                    % special (pre-defined) string functions
                    % you can add here your favorite functions
                    switch lower(OutCol{Ioutcol,2})
                        case 'cenx'
                            % Calculate bin X center
                            OutC(Icat).(CatField)(Ibin,Ioutcol) = EdgeX(Ix)+StepX.*0.5;
                        case 'ceny'
                            % Calculate bin Y center
                            OutC(Icat).(CatField)(Ibin,Ioutcol) = EdgeY(Iy)+StepY.*0.5;
                        otherwise
                            error('Unknown OutCol second column string option');
                    end
                end
            end
            
        end
    end
    % append ColCell and Col to OutC
    OutC(Icat).(ColCellField) = OutCol(:,1)';
    OutC(Icat) = colcell2col(OutC(Icat));
    
end
