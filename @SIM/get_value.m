function [ValPeak,ValInterp,ValMask]=get_value(Sim,XY,varargin)
% Get value of images and mask in a SIM object.
% Package: @SIM
% Description: Get value of images and mask in a SIM object.
%              This includes the pixel peak value, interpolated value, and
%              or operator applied to all pixels within some radius from
%              the requested position.
% Input  : - A SIM object.
%          - List of [X,Y] pixel coordinates.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecFields'   - SIM class fields on which to execute the
%                             operations. Default is {SIM.ImageField}
%                             (i.e., the image field).
%            'InterpMethod' - Interpolation method. Default is 'cubic'.
%            'MaskRadius'   - Radius around position from which to read the
%                             MASK image values. Default is 3 pixels.
%            'Operation'    - Operation on MASK image value.
%                             Default is @Util.array.bitor_array.
%            'CommonSize'   - Is all images have the same size (faster).
%                             Default is true.
% Output : - A cube withe the pixel value at the rounded position.
%            The cube coordinates are [Sim index, field index, coordinate]
%          - A cube with the pixel interpolated values.
%          - A cube with the MASK value.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [ValPeak,ValInterp,ValMask]=get_value(Sim,[100 100]);
% Reliable: 
%--------------------------------------------------------------------------



DefV.ExecFields           = {SIM.ImageField};
DefV.InterpMethod         = 'cubic';
DefV.MaskRadius           = 3;   % radius if NaN than use 'nearest'
DefV.Operation            = @Util.array.bitor_array;
DefV.CommonSize           = true;    % assume all images have common size
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (~iscell(InPar.ExecFields))
    InPar.ExecFields = {InPar.ExecFields};
end

X = XY(:,1);
Y = XY(:,2);
RoundX = round(X);
RoundY = round(Y);

if (InPar.CommonSize)
    Size    = size(Sim(1).(InPar.ExecFields{1}));
    CellInd = ImUtil.Im.find_within_radius_cell(Size,X,Y,InPar.MaskRadius,true);
end

Nxy       = size(XY,1);
Nf        = numel(InPar.ExecFields);
Nsim      = numel(Sim);
ValPeak   = zeros(Nsim,Nf,Nxy);
ValInterp = zeros(Nsim,Nf,Nxy);
ValMask   = zeros(Nsim,Nf,Nxy);
for Isim=1:1:Nsim
    for If=1:1:Nf
        Size    = size(Sim(Isim).(InPar.ExecFields{If}));
        Ind     = Util.array.sub2ind_fast(Size,RoundY,RoundX);
        ValPeak(Isim,If,:) = Sim(Isim).(InPar.ExecFields{If})(Ind);
        if (nargout>1)
            ValInterp(Isim,If,:) = interp2(Sim(Isim).(InPar.ExecFields{If}),Y,X,InPar.InterpMethod);
            if (nargout>2)
                if (~InPar.CommonSize)
                    %Size    = size(Sim(Isim).(InPar.ExecFields{If}));
                    CellInd = ImUtil.Im.find_within_radius_cell(Size,X,Y,InPar.MaskRadius,true);
                end
                ValMask(Isim,If,:) = mask4src(Sim(Isim),CellInd,InPar.Operation);
            end
        end
    end
end
        
    
    