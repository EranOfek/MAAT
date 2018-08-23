function ztf_bias_flat(List,varargin)
% Prepare Bias and Flat images for ZTF
% Package: ImUtil.ZTF
% Description: Prepare Bias and Flat images for ZTF
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Nov 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: ImUtil.ZTF.ztf_bias_flat('ztf_*.fits.fz');
% Reliable: 
%--------------------------------------------------------------------------


DefV.StartDate            = -Inf;
DefV.EndDate              = Inf;
DefV.ListCCDID            = (1:1:16)';
DefV.Namp                 = 4;
DefV.ReadRawImagePar      = {};
DefV.BiasPar              = {};
DefV.FlatPar              = {};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% Create list of file names
[~,List] = Util.files.create_list(List,NaN);

% select images by CCDID, amplifier and filter
Data = VO.ZTF.ztf_imagename_prop(List);

if (numel(InPar.StartDate)==1)
    StartJD = InPar.StartDate;
else
    StartJD = celestial.time.julday(InPar.StartDate);
end
if (numel(InPar.EndDate)==1)
    EndJD = InPar.EndDate;
else
    EndJD = celestial.time.julday(InPar.EndDate);
end

FlagJD = [Data.JD]>StartJD & [Data.JD]<EndJD;
Data   = Data(FlagJD);
List   = List(FlagJD);



% for each CCDID
for IlistCCD=1:1:numel(InPar.ListCCDID)
    % index of CCD
    Iccd = InPar.ListCCDID(IlistCCD);
    
    % index of images in Data for the specific CCD
    IndCCD  = find([Data.CCD]==Iccd);
    % subset of Data for the specific CCD
    DataCCD = Data(IndCCD); 
    ListCCD = List(IndCCD);
    
    % find bias images
    Ibias = strcmp({DataCCD.ImExt},'b.fits.fz');
    
    % load bias images to SIM
    % and subtract overscan bias
    [SimBias,~,OverScan] = ImUtil.ZTF.read_raw_image(ListCCD,InPar.ReadRawImagePar{:});
    
    % construct bias image for CCD
    % for each amplifier
    for Iamp=1:1:InPar.Namp
        [Bias(Iamp),SummaryBias(Iamp),IsBias(Iamp).IsBias] = bias(SimBias(Iamp,:),InPar.BiasPar{:});
    end
    
    % search for all flat images
    Iflat = strcmp({DataCCD.ImExt},'f.fits.fz');
    DataCCDFlat = DataCCD(Iflat);
    ListCCDFlat = ListCCD(Iflat);
    
    % seach for all filters
    UniqueFilterID = unique([DataCCDFlat.FilterID]);
    % remove NaN filters (i.e., bias/dark)
    UniqueFilterID = UniqueFilterID(~isnan(UniqueFilterID));
    NuniqueFilter  = numel(UniqueFilterID);
    
    % for each filter
    for Ifilter=1:1:NuniqueFilter
        % search image of the same filter
        FilterID   = UniqueFilterID(Ifilter);  % current filter ID
        FlagFilter = [DataCCDFlat.FilterID] == FilterID;
        
        DataCCDFlatFilter = DataCCDFlat(FlagFilter);
        ListCCDFlatFilter = ListCCDFlat(FlagFilter);
        
        
        % load flat images to SIM
        [SimFlat,~,OverScanFlat] = ImUtil.ZTF.read_raw_image(ListCCDFlatFilter,InPar.ReadRawImagePar{:});
        % correct for gain
        SimFlat = gain_correct(SimFlat);
        
        for Iamp=1:1:InPar.Namp
            [Flat(Iamp),SummaryFlat(Iamp),IsFlat(Iamp).IsFlat]=flat(SimFlat(Iamp,:),InPar.FlatPar{:});
        end
        
        'a'
        
    end
        
end

    % read flat images
    
    % subtract bias
    
    % construct flat
    
    
    
    
    
    
