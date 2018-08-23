function [Sim,SimDet,SimPeaks]=threshold(Sim,Threshold,varargin)
% Threshold a SIM image above a certain level.
% Package: @SIM
% Description: Threshold a SIM image above a certain level, and look for
%              peaks in the thresholded image. 
% Input  : - A SIM object.
%          - Threshold. This is one of the following:
%            A scalar that will be used as a threshold.
%            A matrix that will be used as a threshold.
%            If 'UseErrIm' argument is set to true (default), then the
%            threshold is the 'ErrIm' field in the SIM multiplied by the
%            threshold (i.e., this is the number of sigmas).
%            Default is 5.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField'- A string containing the field name in the SIM
%                        object on which to perform the thresholding.
%                        Default is 'Im'.
%            'UseErrIm'- A flag indicating if to use the 'ErrIm' field
%                        multiplied by the threshold parameter as the
%                        actual threshold (true) or to use the threshold
%                        parameter as is (false).
%                        Default is true.
%            'CombineFilter' - If tue then combine all the images in the
%                        SIM to produce a single threshold image (based on
%                        the maxima of all images).
%                        If false then threshold each image separatly.
%                        Default is false.
%            'MinArea' - Min area of source above threshold in the filtered
%                        image. Default is 1. Note that values above 1 will 
%                        change the meaning of the detection threshold.
%            'AreaOpenConn'- connectivity for bwareaopen.m Default is 8.
%            'RegionMaxConn'- connectivity for imregionalmax. Default is 8.
%            'ReplaceVal'- If not empty then replace NaN values in the
%                        MaxImage by this scalar. Default is [].
%            'ColXY'   - A cell array of the X and Y position columns in
%                        which to write
%                        the X and Y peak positions in the AstCat object.
%                        Default is {'XPIX_PEAK','YPIX_PEAK'}.
% Output : - The original SIM object, with the found threshold peaks in the
%            AstCat object.
%            Note that if 'CombineFilter' is true then the results will be
%            in the first element of the SIM.
%          - A SIM object with the logical thresholded image.
%          - A SIM object with the PeaksIm image.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S=threshold(S,5);
% Reliable: 2
%--------------------------------------------------------------------------

ImageField      = 'Im';
ErrField        = 'ErrIm';
CatField        = 'Cat';
ColCellField    = 'ColCell';
ColField        = 'Col';


Def.Threshold = 5;
if (nargin==1)
    Threshold  = Def.Threshold;
end
if (isempty(Threshold))
    Threshold  = Def.Threshold;
end

DefV.ExecField          = ImageField;   % a single field - not cell!
DefV.UseErrIm           = true;
DefV.CombineFilter      = false;
DefV.MinArea            = 1;
DefV.AreaOpenConn       = 8;  % [4 | 8]
DefV.RegionMaxConn      = 8;  % [4 | 8]
DefV.ReplaceVal         = []; %[-Inf 0];
DefV.ColXY              = {'XPIX_PEAK','YPIX_PEAK'};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% allocate output arguments
if (nargout>1)
    SizeSim = size(Sim);
    SimDet = SIM(SizeSim);
    if (nargout>2)
        SimPeaks = SIM(SizeSim);
    end
end


Nsim = numel(Sim);


if (InPar.CombineFilter)
    % CombineFilter == true
    % In this mode all the elemnts in the SIM will be combined to a single
    % PeaksIm based on the maximum of all individual images.
    % This is useful when you 
    
    PrevMaxImage = zeros(size(Sim(1).(InPar.ExecField)));
    for Isim=1:1:Nsim
        % for each SIM elelemt
        
        % threshold the image
        if (InPar.UseErrIm)
            % Thresholding is done relative to the ErrIm field
            DetSrc = Sim(Isim).(InPar.ExecField)>(Threshold.*Sim(Isim).(ErrField));
        else
            % thresholding is done relative to the Threshold parameter
            DetSrc = Sim(Isim).(InPar.ExecField)>Threshold;
        end
        
        % remove detection with number of pixels smaller than MinArea
        % For example if MinArea is 2 DetSrc with 1 pixel will be removed,
        % but those with 2 pixels will be left as they are
        if (InPar.MinArea>1)
            % In most cases the user will request for MinArea=1 so there is no
            % need to call this function
            DetSrc = bwareaopen(DetSrc,InPar.MinArea,InPar.AreaOpenConn);
        end

        if (nargout>1)
            SimDet(Isim).(InPar.ExecField) = DetSrc;
        end
        
        % locate local maxima in the filtered image
        %--- KNOWN PROBLEM if two identical value next to each other than
        %--- imregionalmax will find both
        MaxImage = DetSrc.*Sim(Isim).(InPar.ExecField);
        if (~isempty(InPar.ReplaceVal))
            MaxImage(isnan(MaxImage)) = InPar.ReplaceVal;
        end
        
        % combine all the MaxImage's
        PrevMaxImage = max(MaxImage,PrevMaxImage);
    
    end
    MaxImage = SIM;
    MaxImage.(InPar.ExecField) = PrevMaxImage;
    
    
    MaxImage = local_maxima(MaxImage,InPar.ColXY,InPar.ExecField,InPar.RegionMaxConn);
    
    % The catalog of oeaks will appear only in the first SIM element
    Sim(1).(CatField)     = MaxImage.(CatField);
    Sim(1).(ColCellField) = MaxImage.(ColCellField);
    Sim(1).(ColField)     = MaxImage.(ColField);
    
    
else
    % CombineFilter == false
    % In this model each element in the SIM object will be thresholded
    % seperatly.
    MaxImage = SIM;          % Define a single element SIM object
    for Isim=1:1:Nsim
        % for each SIM element

        % threshold the image
        if (InPar.UseErrIm)
            % Thresholding is done relative to the ErrIm field
            if (~isempty(Sim(Isim).(ErrField)))
                DetSrc = Sim(Isim).(InPar.ExecField)>(Threshold.*Sim(Isim).(ErrField));
            else
                error('UseErrIm option was selected but ErrIm doesnt exist');
            end
        else
            % thresholding is done relative to the Threshold parameter
            DetSrc = Sim(Isim).(InPar.ExecField)>Threshold;
        end

        % remove detection with number of pixels smaller than MinArea
        % For example if MinArea is 2 DetSrc with 1 pixel will be removed,
        % but those with 2 pixels will be left as they are
        if (InPar.MinArea>1)
            % In most cases the user will request for MinArea=1 so there is no
            % need to call this function
            DetSrc = bwareaopen(DetSrc,InPar.MinArea,InPar.AreaOpenConn);
        end

        if (nargout>1)
            SimDet(Isim).(InPar.ExecField) = DetSrc;
        end
        % locate local maxima in the filtered image
        %--- KNOWN PROBLEM if two identical value next to each other than
        %--- imregionalmax will find both
        
        MaxImage.(InPar.ExecField) = DetSrc.*Sim(Isim).(InPar.ExecField);
        if (~isempty(InPar.ReplaceVal))
            MaxImage.(InPar.ExecField)(isnan(MaxImage.(InPar.ExecField))) = InPar.ReplaceVal;
            %MaxImage(isnan(MaxImage)) = InPar.ReplaceVal;
        end

        % locate peaks
        MaxImage = local_maxima(MaxImage,InPar.ColXY,InPar.ExecField,InPar.RegionMaxConn);
        Sim(Isim).(CatField)     = MaxImage.(CatField);
        Sim(Isim).(ColCellField) = MaxImage.(ColCellField);
        Sim(Isim).(ColField)     = MaxImage.(ColField);
        
    end
end