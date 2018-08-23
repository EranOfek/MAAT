function [Sim]=interp2(Sim,X,Y,varargin)
% Apply interp2.m for a SIM object.
% Package: @SIM
% Description: Apply interp2.m for a SIM object.
% Input  : - A SIM object.
%          - A matrix of X coordinates over which to interepolate the
%            images. If this, and the next input arguments are vector,
%            then will construct the X and Y matrices using meshgrid.m.
%          - Like the second argument, but for the Y axis.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField'  - A cell array (or a single string) of field
%                           names on which to execute the function.
%                           Default is {'Im'}.
%            'InterpMethod'-Interpolation method. See interp2.m for
%                           options. Default is 'cubic'.
%            'ExtrapVal'  - Extrapolation value. See interp2.m for
%                           options. Default is {}.
%            'MaskInterpMethod'- Interpolation method for the .Mask field.
%                           See interp2.m for options.
%                           Default is 'nearest'.
%            'ReplaceKey' - A three column cell array of {key,val,comment},
%                           or an HEAD object which keywords to replace
%                           or add to the SIM header.
%            'AddKey'     - Like 'ReplaceKey', but adding keywords,
%                           without replacment.
%            'UseFast'    - Use interp2fast.m {true|false}. This is usefuk
%                           when the interpolated region is much smaller
%                           than the whole image. Default is false.
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - A SIM object with the requested interpolated objects.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Sim]=interp2(S,(1:0.5:100),(1:0.5:100));
% Reliable: 2
%--------------------------------------------------------------------------

ImageField      = 'Im';
MaskField       = 'Mask';

DefV.ExecField          = {ImageField};
DefV.InterpMethod       = 'cubic';
DefV.ExtrapVal          = {};
DefV.MaskInterpMethod   = 'nearest';
DefV.ReplaceKey         = {};
DefV.AddKey             = {};
DefV.UseFast            = false;
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


if (~iscell(InPar.ExecField)),
    InPar.ExecField = {InPar.ExecField};
end



if (isvector(X) && isvector(Y)),
    % construct X/Y matrices:
    [MatX,MatY] = meshgrid(X,Y);
else
    % Assume X and Y are already in matrix form
    MatX = X;
    MatY = Y;
end



Nfn  = numel(InPar.ExecField);
Nsim = numel(Sim);
for Isim=1:1:Nsim,
    % for each SIM element
  
    for Ifn=1:1:Nfn,
        % for each field 
        if (~strcmp(InPar.ExecField{Ifn},MaskField)),
            if (InPar.UseFast),
                Size = size(Sim(Isim).(InPar.ExecField{Ifn}));
                [InX,InY] = meshgrid((1:1:Size(2)),(1:1:Size(1)));
                Sim(Isim).(InPar.ExecField{Ifn}) = interp2fast(InX(1,:),InY(:,1),Sim(Isim).(InPar.ExecField{Ifn}),MatX,MatY,InPar.InterpMethod,InPar.ExtrapVal{:});
            else
                Sim(Isim).(InPar.ExecField{Ifn}) = interp2(Sim(Isim).(InPar.ExecField{Ifn}),MatX,MatY,InPar.InterpMethod,InPar.ExtrapVal{:});
            end
        else
            if (InPar.UseFast),
                Size = size(Sim(Isim).(InPar.ExecField{Ifn}));
                [InX,InY] = meshgrid((1:1:Size(2)),(1:1:Size(1)));
                Sim(Isim).(InPar.ExecField{Ifn}) = interp2fast(InX(1,:),InY(:,1),Sim(Isim).(InPar.ExecField{Ifn}),MatX,MatY,InPar.MaskInterpMethod,InPar.ExtrapVal{:});
            else
                Sim(Isim).(InPar.ExecField{Ifn}) = interp2(Sim(Isim).(InPar.ExecField{Ifn}),MatX,MatY,InPar.MaskInterpMethod,InPar.ExtrapVal{:});
            end
        end
        
    end
    
        %--- Update header ---
    % replace
    if (~isempty(InPar.ReplaceKey)),
        Sim(Isim) = replace_key(Sim(Isim),InPar.ReplaceKey);
       
    end
    % add
    if (~isempty(InPar.AddKey)),
        Sim(Isim) = add_key(Sim(Isim),InPar.AddKey);
      
    end
end

