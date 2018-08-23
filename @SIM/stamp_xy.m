function [Out]=stamp_xy(Sim,XY,varargin)
% Get image stamps around requested positions
% Package: @SIM
% Description: Given a SIM object, and a single list of multiple X/Y
%              positions, get image stamps around requested positions,
%              including sub pixel shifts.
%              Regions outside image boundries are always padded.
% Input  : - A  SIM object. Number of elements should be either 1 or equal
%            to the number of requested positions.
%            If a single element than retrive multiple stamps from the same
%            image.
%          - Two column matrix of [X,Y] positions to extract.
%            Alternativly, this can be a single AstCat object with the
%            X, Y positions, where X Y column names are indicated by
%            'PosCatCol'.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField'  - A cell array (or a single string) of field
%                           names on which to execute the function.
%                           Default is {'Im'}.
%            'StampSize' - Stamp size is 1 + twice this parameter.
%                          Default is 10.
%            'Align'     - Align stamps to sub pixel position using
%                          requested coordinates. If false than will get
%                          stamp around rounded value of coordinates.
%                          Alignment is done using ImUtil.Im.image_shift_fft.m.
%                          Default is true.
%            'HalfPixShift'- Shift the stamp by additional half pixel.
%                          This controls if the source center is in the
%                          pixel center or on a pixel edge.
%                          Default is true.
%            'OutType'   - Output type:
%                          'sim' - A SIM object.
%                          'cube'- A cube in which the image index is the
%                                  first dimension.
%                          Default is 'cube'.
%            'PadValue'  - Pad value. Default is NaN.
%            'UpdateWCS' - Update WCS of SIM output. Default is true.
%            'PosCatCol' - The column names in the AstCat object indicating
%                          the X/Y positions.
%                          Default is {'XWIN_IMAGE','YWIN_IMAGE'}.
% Output : - Either a SIM object with the multiple stamps, or a cube of
%            stamps in which the image index is the first dimension.
% See als: SIM/trim_image.m, SIM/stamp_coo.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Out=stamp_xy(S1,XY);
% Reliable: 2
%--------------------------------------------------------------------------

ImageField      = 'Im';
%BackField       = 'BackIm';
CatField        = 'Cat';

DefV.ExecField          = {ImageField};
DefV.StampSize          = 10;
DefV.Align              = true;
DefV.HalfPixShift       = true;
DefV.OutType            = 'cube';   % 'cube'|'sim'
DefV.PadValue           = NaN;
DefV.UpdateWCS          = true;
DefV.PosCatCol          = {'XWIN_IMAGE','YWIN_IMAGE'};
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (~iscell(InPar.ExecField))
    InPar.ExecField = {InPar.ExecField};
end
Nf = numel(InPar.ExecField);

if (strcmp(InPar.OutType,'cube') && Nf>1)
    error('For cube output number of ExecField must be 1');
end

% Convert AstCat input into a two column matrix of [X,Y] positions
if (AstCat.isastcat(XY))
    if (numel(XY)>1)
        error('AstCat object must contain a single element');
    end
    ColInd = colname2ind(XY,InPar.PosCatCol);
    XY = XY.(CatField)(:,ColInd);
end

% HalfPixShift responsible to put the star center on a pixel edge rather
% than on pixel center.
if (InPar.HalfPixShift)
    HalfPixShift = -0.5;
else
    HalfPixShift = 0;
end

Nxy   = size(XY,1);

if (Nxy==0)
    error('Empty list of coordinates');
end
Nsim  = numel(Sim);
Nobj  = max(Nxy,Nsim);
Shift = XY - InPar.StampSize - HalfPixShift;   % Note the sign of HalfShiftPix is negative!

for Iobj=1:1:Nobj
   
    % for each object:
    % either SIM or X/Y position
    Ixy  = min(Iobj,Nxy);
    Isim = min(Iobj,Nsim);
    for If=1:1:Nf
        XY_I = XY(Ixy,:);
        % trim by whole pixel
        XY_Trim  = round(XY_I);
        % Shift by sub pixel
        % This is the Shift of the source rather than the shift of the
        % image
        XY_Shift = XY_Trim - XY_I + HalfPixShift; % the sign here was checked
       
        [Image,~,~,~] = ImUtil.Im.trim_image(Sim(Isim).(InPar.ExecField{If}),[XY_Trim, InPar.StampSize, InPar.StampSize],'center',InPar.PadValue);
        if (InPar.Align)
            Image                      = ImUtil.Im.image_shift_fft(Image,XY_Shift(1),XY_Shift(2));
        end
        switch lower(InPar.OutType)
            case 'sim'
                if (Iobj==1)
                    Out = simdef(Nobj,1);
                end
                Out(Iobj) = Sim(Isim);
                
                Out(Iobj).(InPar.ExecField{If}) = Image;
                
                %--- Updtae WCS ---
                if (InPar.UpdateWCS)
                    ValCRPIX  = mgetkey(Out(Isim),{'CRPIX1','CRPIX2'});
                    Out(Iobj) = update_key(Out(Iobj),{'CRPIX1',ValCRPIX{1}-Shift(Ixy,1);...
                                                      'CRPIX2',ValCRPIX{2}-Shift(Ixy,2)});
                    if (~isempty_wcs(Out(Isim)))
                        Out(Iobj) = populate_wcs(Out(Iobj));
                    end
                end
                
            case 'cube'
                if (Iobj==1)
                    % initiate Cube
                    Out = zeros(Nobj,2.*InPar.StampSize+1,2.*InPar.StampSize+1);
                end
                Out(Iobj,:,:) = Image;
                
                
            otherwise
                error('Unknown OutType option');
        end
    end
    
    
end
                
                
    


