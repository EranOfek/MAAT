function Sim=gradient2mask(Sim,varargin)
%--------------------------------------------------------------------------
% gradient2mask function                                        class/@SIM
% Description: Calculate the gradient of an image and/or a background
%              image in a SIM.
%              If the gradient is larger than some threshold value
%              then populate the mask image.
% Input  : - A SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Bit_BackGrad'
%            'Bit_ImGrad'
%            
% Output : - 
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jul 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S1=gradient2mask(S1)
% Reliable: 
%--------------------------------------------------------------------------

ImageField     = 'Im';
BackField      = 'BackIm';
ErrField       = 'ErrIm';


DefV.Bit_BackGrad       = 'Bit_BackGrad';
DefV.Bit_ImGrad         = 'Bit_ImGrad';
DefV.ImGradX            = [];
DefV.ImGradY            = [];
DefV.BackGradX          = [];
DefV.BackGradY          = [];
DefV.BackBufferSize     = [64 64];
DefV.ImThresh           = [40 40];   % [x,y]
DefV.ImThreshUnits      = 'sigma';   % 'sigma'|'count'
DefV.BackThresh         = [5 5];     % [x,y]
DefV.BackThreshUnits    = 'sigma';   % 'sigma'|'count'
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (numel(InPar.ImThresh)==1)
    InPar.ImThresh = InPar.ImThresh.*ones(2,1);
end
if (numel(InPar.BackThresh)==1)
    InPar.BackThresh = InPar.BackThresh.*ones(2,1);
end
if (numel(InPar.BackBufferSize)==1)
    InPar.BackBufferSize = InPar.BackBufferSize.*ones(2,1);
end

BufferSize = sqrt(prod(InPar.BackBufferSize));

Nsim = numel(Sim);
for Isim=1:1:Nsim
    %-----------------------
    %--- Image gradients ---
    %-----------------------
    if (~isempty(InPar.Bit_ImGrad))
        % User requested to populate the Bit_ImGrad bit
        if (isempty(InPar.ImGradX) || isempty(InPar.ImGradY))
            % Calculate the gradient of the Image
            [ImGx,ImGy] = gradient(Sim(Isim).(ImageField));
        else
            ImGx = InPar.ImGradX;
            ImGy = InPar.ImGradY;
        end

        switch lower(InPar.ImThreshUnits)
            case 'sigma'
                ErrIm  = Sim(Isim).(ErrField);
            case 'count'
                ErrIm  = 1;
            otherwise
                error('Unknown ImThreshUnits option');
        end
        ImFlag = abs(ImGx)./ErrIm>InPar.ImThresh(1) | ...
                 abs(ImGy)./ErrIm>InPar.ImThresh(2);

        % populate the mask image
        Sim(Isim) = bitmask_set(Sim(Isim),ImFlag,InPar.Bit_ImGrad);
    end
    
    %----------------------------
    %--- Background gradients ---
    %----------------------------
    if (~isempty(InPar.Bit_BackGrad))
        % User requested to populate the Bit_BackGrad bit
        if (isempty(InPar.BackGradX) || isempty(InPar.BackGradY))
            % Calculate the gradient of the Background
            [BackGx,BackGy] = gradient(Sim(Isim).(BackField));
        else
            BackGx = InPar.BackGradX;
            BackGy = InPar.BackGradY;
        end
        BackG = sqrt(BackGx.^2 + BackGy.^2);
        
        switch lower(InPar.BackThreshUnits)
            case 'sigma'
                [Mode,Std] = mode_fit(BackG);
                
            case 'count'
               Mode = 0;
               Std  = 1;
            otherwise
                error('Unknown BackThreshUnits option');
        end
        BackFlag = BackGx>(Mode+Std.*InPar.BackThresh(1)) | ...
                           BackGy>(Mode+Std.*InPar.BackThresh(2));
        
        % populate the mask image
        Sim(Isim) = bitmask_set(Sim(Isim),BackFlag,InPar.Bit_BackGrad);
    end
end
