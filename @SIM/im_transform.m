function [Sim]=im_transform(Sim,AstT,varargin)
%--------------------------------------------------------------------------
% im_transform function                                         class/@SIM
% Description: Transformation of SIM images object.
% Input  : - A SIM object.
%          - An AstTran object containing the transformation to apply to
%            the images.
%            Alternatively this is a cell array of AstTran objects,
%            each cell corresponds to a SIM element.
%            Or this is a structure array (element in the array for each
%            element in the SIM) which the field designated in
%            'TranField' (see below) contains the AstTran object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField' - The SIM fields on which to execute the operator.
%                          Default is {'Im'}.
%            'InterpMethod'- Transformation interpolation method using
%                          interp2.m. Default is 'linear'.
%            'TranField' - If the transformation is a structure than this
%                          indicate the field name in the structure array
%                          containing the AstTran object.
%                          Default is 'AstT'.
%            'DeleteWCS' - Delete WCS from header. Default is false.
%                          Note that after the transformation the header
%                          WCS will be incorrect.
% Output : - A SIM object with the transformed images.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstT = cell2asttran({'x_shift',1;'y_shift',10});
%          [Sim]=im_transform(S(2),AstT);
%          AstT = cell2asttran({'x_shift',1;'y_shift',10;'xy_rot',[1 0]});
% Reliable: 
%--------------------------------------------------------------------------

ImageField   = 'Im';


DefV.ExecField            = {ImageField};
DefV.InterpMethod         = 'linear';
DefV.TranField            = 'AstT';
DefV.DeleteWCS            = false;
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


ImSize = imagesize(Sim);   % [X,Y] and not [I,J]

Nsim = numel(Sim);
Nf   = numel(InPar.ExecField);

for Isim=1:1:Nsim,
    % for each SIM element
    
    % decide if a new mesh grid is needed
    % a new mesh grid is needed if the first image
    % or image size changed.
    if (Isim==1),
        CalcMesh = true;
    else
        if all(ImSize(Isim,:)==ImSize(Isim-1,:)),
            CalcMesh = false;
        else
            CalcMesh = true;
        end
    end
    
    % calculate mesh grid
    if (CalcMesh),
        [MatX,MatY] = meshgrid((1:1:ImSize(Isim,1)),(1:1:ImSize(Isim,2)));
    end
       
    if (iscell(AstT)),
        Tran = AstT{Isim};
    elseif (isstruct(AstT)),
        Tran = AstT(Isim).(TranField);
    elseif (isasttran(AstT)),
        Tran = AstT;
    else
        error('Unknown AstT input type');
    end
    
    not working...
        
    T = maketform('custom', 2, 2, [], INVERSE_FCN, TDATA)
    
    for If=1:1:Nf,
        % for each field in the SIM element
        
        % transform pixel coordinates
        TransXY = transform(Tran,[MatX(:),MatY(:)],'OutType','mat');
        % transform image using interp2
        Sim(Isim).(InPar.ExecField{If}) = reshape(interp2(MatX,MatY,...
                                                          Sim(Isim).(InPar.ExecField{If}),...
                                                          TransXY(:,1),TransXY(:,2),InPar.InterpMethod),...
                                                  ImSize(Isim,2),ImSize(Isim,1));
    end
    
    % delete WCS
    if (InPar.DeleteWCS),
        % delete WCS from header
        Sim(Isim) = delete_wcs(Sim(Isim));
    end
    
end
    
        
        
        

