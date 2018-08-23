function Sim=image_shift_fft(Sim,Shift,varargin)
% Shift image pixels using the Fourier transform shift theorem.
% Package: @SIM
% Description: Shift image pixels using the Fourier transform
%              shift theorem. The function also update the WCS and the
%              associated catalog.
%              This function uses ImUtil.Im.image_shift_fft.m.
%              This works well when the image does not contain sharp
%              artifacts. Sharp artifacts will produce ringing.
%              Note that the shift is defined on the content of the image,
%              rather than the image boundries - e.g., the stars will be
%              shifted in the requested direction.
% Input  : - A SIM object.
%          - A two column matrix of shifts to apply to each image.
%            Each row corresponds to image in the SIM array. If a single
%            row, then apply the same shift for all images.
%            First column is for the X shift, and second column is for the
%            Y shift. The shift is defined on the content of the image,
%            rather than the image boundries - e.g., the stars will be
%            shifted in the requested direction.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField'  - A cell array (or a single string) of field
%                           names on which to execute the function.
%                           Default is {'Im'}.
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
%            'UpdateWCS'  - Update WCS. Default is true.
%            'UpdateCat'  - Update catalog. Default is true.
%            'PosColX'    - Names of X position columns to update in the
%                           catalog. Default is {'XWIN_IMAGE'}.
%            'PosColY'    - Names of Y position columns to update in the
%                           catalog. Default is {'YWIN_IMAGE'}.
% Output : - Shifted SIM images.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=image_shift_fft(S,[1.6 11.2]);
% Reliable: 2
%--------------------------------------------------------------------------

ImageField      = 'Im';
ColCellField    = 'ColCell';


DefV.ExecField            = {ImageField};
DefV.ReplaceKey           = {};
DefV.AddKey               = {};
DefV.MaskFun              = [];   % should work on SIM! and return a SIM
DefV.MaskFunPar           = {};   % e.g., MaskFun=@eq, MaskFunPar={0}
DefV.UpdateWCS            = true;
DefV.UpdateCat            = true;
DefV.PosColX              = {'XWIN_IMAGE'};
DefV.PosColY              = {'YWIN_IMAGE'};
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (~iscell(InPar.ExecField))
    InPar.ExecField = {InPar.ExecField};
end

Nsim   = numel(Sim);
Nshift = size(Shift,1);
Nf     = numel(InPar.ExecField);
for Isim=1:1:Nsim
    % for each SIM element
    
    % define shift for SIM:
    Ishift = min(Isim,Nshift);
    DX = Shift(Ishift,1);
    DY = Shift(Ishift,2);
    
    for If=1:1:Nf
        % for each SIM field
        Sim(Isim).(InPar.ExecField{If}) = ImUtil.Im.image_shift_fft(Sim(Isim).(InPar.ExecField{If}),DX,DY);
        
    end
            
    %--- Updtae WCS ---
    if (InPar.UpdateWCS)
        ValCRPIX  = mgetkey(Sim(Isim),{'CRPIX1','CRPIX2'});
        warning('update header in image_shift_fft wasnt tested - check sign!');
        Sim(Isim) = update_key(Sim(Isim),{'CRPIX1',ValCRPIX{1}-DX;...
                                          'CRPIX2',ValCRPIX{2}-DY});
        if (~isempty_wcs(Sim(Isim)))
            Sim(Isim) = populate_wcs(Sim(Isim));
        end
    end
        
    %--- Update catalog ---
    if (InPar.UpdateCat)
        if (~isempty(Sim(Isim).(ColCellField)))
            % Update the coordinates in the SIM AstCat object
            % Shift by DX, DY:
            Sim(Isim) = col_fun(Sim(Isim),@plus,InPar.PosColX,[],{DX});
            Sim(Isim) = col_fun(Sim(Isim),@plus,InPar.PosColY,[],{DY});
            
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

        
        