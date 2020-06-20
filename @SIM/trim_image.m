function [Sim,Shift,Sec,SecF]=trim_image(Sim,TrimSec,varargin)
% Trim a region from images in a SIM object
% Package: @SIM
% Description: Trim a region from images in a SIM object. If the trim 
%              region is outside the image bounderies, then the trimmed
%              image can be optionally padded to have the size of the
%              requested trim region.
%              Also update the WCS header.
%              This function can generate a trim section per image or
%              multiple trim sections from a single image.
% Input  : - A SIM object.
%          - A trim section: [Xmin, Xmax, Ymin, Ymax],
%            or [Xcenter, Ycenter, Xhalf_width, Yhalf_width],
%            or a string containing an header keyword name that contains
%            a CCDSEC like region (i.e., read using ccdsec.m).
%            If multipl lines are provided then assume that each line
%            corresponds to a SIM element.
%            If the number of trim section is larger than 1 and number of
%            SIM elements is 1, then retrive multiple trim section from a
%            single image.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'SectionMethod'- Trim section representation method:
%                       'section' - Trim section format is: [Xmin, Xmax, Ymin, Ymax].
%                       'center'  - Format: [Xcenter, Ycenter, Xhalf_width, Yhalf_width].
%                       Default is 'section'.
%            'SectionPad' - Pad value for out of bounderies regions.
%                       If empty then do not pad. Default is NaN.
%            'ExecField'  - A cell array (or a single string) of field
%                           names on which to execute the function.
%                           Default is {'Im','BackIm','ErrIm','WeightIm','Mask'}.
%            'ReplaceKey' - A three column cell array of {key,val,comment},
%                           or an HEAD object which keywords to replace
%                           or add to the SIM header.
%            'AddKey'     - Like 'ReplaceKey', but adding keywords,
%                           without replacment.
%            'MaskFun'    - Function that sets the bit mask:
%                           Fun(Sim1,MaskFunPar{:}).
%                           Default is empty.
%            'MaskFunPar' - Additional parameters to pass to MaskFun.
%                           Default is {}.
%            'UseNaN'     - A flag indicating what to do if CCDSEC value
%                           is NaN. If true then return NaN, if false then
%                           get CCDSEC using naxis or imagesize (of input
%                           is SIM). Default is false.
%            'GetKeyMethod'-Method by which to select the best keyword
%                           value. See getkey_fromlist.m for options.
%                           Default is 'first'.
%            'UpdateWCS'  - Update WCS. Default is true.
%            'UpdateCat'  - A logical flag indicating if to update the
%                           X,Y coordinates in the associated catalog.
%                           Default is true.
%            'TrimCat'    - Trim catalog. Default is true.
%            'ColX'       - Cell array of X coordinate column names to which
%                           to apply the X transformation.
%                           Default is
%                           {'XWIN_IMAGE','X','XPEAK_IMAGE','X_IMAGE'}.
%            'ColY'       - Cell array of Y coordinate column names to which
%                           to apply the Y transformation.
%                           Default is
%                           {'YWIN_IMAGE','Y','YPEAK_IMAGE','Y_IMAGE'}.
% Output : - A SIM object with the trimmed images and updated WCS header.
% See also: trim_image.m, SIM/stamp_coo.m, SIM/stamp_xy.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: SimT=trim_image(Sim,[-1 2 0 4]);
% Reliable: 2
%--------------------------------------------------------------------------



ImageField      = Sim.ImageField;
CatField        = AstCat.CatField;

DefV.SectionMethod        = 'section';
DefV.SectionPad           = NaN;
DefV.ExecField            = {'Im','BackIm','ErrIm','WeightIm','Mask'}; %{ImageField};
%DefV.CCDSEC               = [];
DefV.ReplaceKey           = {};
DefV.AddKey               = {};
DefV.MaskFun              = [];   % should work on SIM! and return a SIM
DefV.MaskFunPar           = {};   % e.g., MaskFun=@eq, MaskFunPar={0}
DefV.UseNaN               = false;
DefV.GetKeyMethod         = 'first';
DefV.UpdateWCS            = true;
DefV.UpdateCat            = true;
DefV.TrimCat              = true;
DefV.ColX                 = {'XWIN_IMAGE','X','XPEAK_IMAGE','X_IMAGE'};
DefV.ColY                 = {'YWIN_IMAGE','Y','YPEAK_IMAGE','Y_IMAGE'};
%DefV.MaskDic              = @MASK.def_bitmask_pipeline;
%DefV.MaskBit              = [];    % bit index or bit name - single or many
if (isempty(varargin))
    InPar = DefV;
else
    %InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
    InPar = InArg.populate_keyval(DefV,varargin,mfilename);
end

if (~iscell(InPar.ExecField))
    InPar.ExecField = {InPar.ExecField};
end

if (ischar(TrimSec))
    Ntrimsec = 1;
else
    Ntrimsec = size(TrimSec,1);
end
Nfn   = numel(InPar.ExecField);
Nsim  = numel(Sim);

if (Nsim==1 && Ntrimsec>1)
    %--- One input image and multiple trims ---
    Sim1 = Sim;
    Sim  = simdef(Ntrimsec,1);
    
    Shift = zeros(Ntrimsec,2);
    Sec   = zeros(Ntrimsec,4);
    SecF  = zeros(Ntrimsec,4);
   
    for Itrimsec=1:1:Ntrimsec
        % for each Trim section
        Isim = Itrimsec;    % for simplicity of copying code
        Sim(Itrimsec) = Sim1;
        
        % TrimSec is numeric
        Section = TrimSec(Itrimsec,:);
        
        for Ifn=1:1:Nfn
            % for each field 

            if (numel(Sim.(InPar.ExecField{Ifn}))>1)

                [Sim(Itrimsec).(InPar.ExecField{Ifn}),Shift(Itrimsec,:),Sec(Itrimsec,:),SecF(Itrimsec,:)] = ImUtil.Im.trim_image(Sim1.(InPar.ExecField{Ifn}),...
                                                              Section,...
                                                              InPar.SectionMethod,...
                                                              InPar.SectionPad);
            end

        end

        %--- Updtae WCS ---
        if (InPar.UpdateWCS)
            ValCRPIX  = mgetkey(Sim(Isim),{'CRPIX1','CRPIX2'});
            Sim(Isim) = update_key(Sim(Isim),{'CRPIX1',ValCRPIX{1}-Shift(Isim,1);...
                                              'CRPIX2',ValCRPIX{2}-Shift(Isim,2)});
            if (~isempty_wcs(Sim(Isim)))
                Sim(Isim) = populate_wcs(Sim(Isim));
            end
        end

        %--- Update mask ---
        if (~isempty(InPar.MaskFun))
            Sim(Isim) = InPar.MaskFun(Sim(Isim),InPar.MaskFunPar{:});
        end

        %--- Update header ---
        % replace
        if (~isempty(InPar.ReplaceKey))
            Sim(Isim) = replace_key(Sim(Isim),InPar.ReplaceKey);
        end
        % add
        if (~isempty(InPar.AddKey))
            Sim(Isim) = add_key(Sim(Isim),InPar.AddKey);
        end
    end 
else
    %--- Image per trim ---
    Shift = zeros(Nsim,2);
    Sec   = zeros(Nsim,4);
    SecF  = zeros(Nsim,4);

    for Isim=1:1:Nsim
        % for each SIM element
        %--- ccdsec ---
        if (ischar(TrimSec))
            Section = ccdsec(Sim(Isim),TrimSec,InPar.UseNaN,InPar.GetKeyMethod);
        else
            % if more than one line in TrimSec then select the Isim line:
            Section = TrimSec(min(Ntrimsec,Isim),:);
        end

        %--- call trim_image.m ---
        for Ifn=1:1:Nfn
            % for each field 
            if (numel(Sim(Isim).(InPar.ExecField{Ifn}))>1)
                [Sim(Isim).(InPar.ExecField{Ifn}),Shift(Isim,:),Sec(Isim,:),SecF(Isim,:)] = ImUtil.Im.trim_image(Sim(Isim).(InPar.ExecField{Ifn}),...
                                                              Section,...
                                                              InPar.SectionMethod,...
                                                              InPar.SectionPad);
            end
        end

        %--- Updtae WCS ---
        if (InPar.UpdateWCS)
            ValCRPIX  = mgetkey(Sim(Isim),{'CRPIX1','CRPIX2'});
            Sim(Isim) = update_key(Sim(Isim),{'CRPIX1',ValCRPIX{1}-Shift(Isim,1);...
                                              'CRPIX2',ValCRPIX{2}-Shift(Isim,2)});
            if (~isempty_wcs(Sim(Isim)))
                Sim(Isim) = populate_wcs(Sim(Isim));
            end
        end

        %--- Update mask ---
        if (~isempty(InPar.MaskFun))
            Sim(Isim) = InPar.MaskFun(Sim(Isim),InPar.MaskFunPar{:});
        end

        %--- Update header ---
        % replace
        if (~isempty(InPar.ReplaceKey))
            Sim(Isim) = replace_key(Sim(Isim),InPar.ReplaceKey);
        end
        % add
        if (~isempty(InPar.AddKey))
            Sim(Isim) = add_key(Sim(Isim),InPar.AddKey);
        end
        
        %--- trim catalog ---
        if (InPar.TrimCat)
            
            if (isfield_populated(Sim(Isim),CatField))
                
                NcolX = numel(InPar.ColX);
                NcolY = numel(InPar.ColY);
                Ncol  = max(NcolX,NcolY);
                for Icol=1:1:Ncol
                    IcolX = min(NcolX,Icol);
                    IcolY = min(NcolY,Icol);

                    ColX  = colname2ind(Sim(Isim),InPar.ColX{IcolX});
                    ColY  = colname2ind(Sim(Isim),InPar.ColY{IcolY});
                    if (~isnan(ColX) && ~isnan(ColY))
                        
                        % for X dimension
                        Flag = Sim(Isim).(CatField)(:,ColX) > Section(1) & ...
                               Sim(Isim).(CatField)(:,ColX) < Section(2) & ...
                               Sim(Isim).(CatField)(:,ColY) > Section(3) & ...
                               Sim(Isim).(CatField)(:,ColY) < Section(4);
                        Sim(Isim).(CatField) = Sim(Isim).(CatField)(Flag,:);
                        
                        
                        
                    end
                end
                
            end
            
            
        end
        
        %--- Update Catalog ---
        if (InPar.UpdateCat)
            
            if (isfield_populated(Sim(Isim),CatField))
                
                NcolX = numel(InPar.ColX);
                NcolY = numel(InPar.ColY);
                Ncol  = max(NcolX,NcolY);
                for Icol=1:1:Ncol
                    IcolX = min(NcolX,Icol);
                    IcolY = min(NcolY,Icol);

                    ColX  = colname2ind(Sim(Isim),InPar.ColX{IcolX});
                    ColY  = colname2ind(Sim(Isim),InPar.ColY{IcolY});
                    if (~isnan(ColX) && ~isnan(ColY))
                        
                        % for X dimension
                        Sim(Isim).(CatField)(:,ColX) = Sim(Isim).(CatField)(:,ColX) - Shift(Isim,1);
                        % for Y dimension
                        Sim(Isim).(CatField)(:,ColY) = Sim(Isim).(CatField)(:,ColY) - Shift(Isim,2);
                        
                    end
                end
                
            end
            
            
        end
        
        
        
    end

end
