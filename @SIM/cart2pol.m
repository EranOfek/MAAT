function Sim=cart2pol(Sim,Center,OutTheta,OutR,varargin)
% Convert SIM image from cartezian coordinates to polar coordinates.
% Package: @SIM
% Description: Convert SIM image from cartezian coordinates to polar
%              coordinates.
% Input  : - A SIM object.
%          - Center [X, Y] to use of the origin on the polar coordinates.
%          - Vector of the requested output polar coordinates theta
%            (radians).
%          - Vector of the requested output polar coordinates radius
%            (pixels). 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'X'          - A vector of X coordinates represing the images
%                           X coordinates. If empty then use (1:1:SizeX).
%                           Default is empty.
%            'Y'          - A vector of Y coordinates represing the images
%                           Y coordinates. If empty then use (1:1:SizeY).
%                           Default is empty.
%            'InterpMethod'- Interpolation method. See interp2.m for
%                           options. Default is 'linear'.
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
% Output : - A SIM object containing the images in polar coordinates.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sp=cart2pol(Sim,[512 512],(1:1:360),(1:1:512));
% Reliable: 
%--------------------------------------------------------------------------


DefV.X                  = [];
DefV.Y                  = [];
DefV.InterpMethod       = 'cubic';
DefV.ExecField          = {'Im'};
DefV.ReplaceKey         = {};
DefV.AddKey             = {};
DefV.MaskFun            = [];   % should work on SIM! and return a SIM
DefV.MaskFunPar         = {};   % e.g., MaskFun=@eq, MaskFunPar={0}
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

[MatOutTheta,MatOutR] = meshgrid(OutTheta,OutR);

if (~iscell(InPar.ExecField))
    InPar.ExecField = {InPar.ExecField};
end
Nf = numel(InPar.ExecField);

Nsim = numel(Sim);
OldVecX = [];
OldVecY = [];
for Isim=1:1:Nsim
    % for each SIM element
    for If=1:1:Nf
        % for each SIM field
        
        [Sy,Sx] = size(Sim(Isim).(InPar.ExecField{If}));
        if (isempty(InPar.X))
            VecX = (1:1:Sx);
        else
            VecX = InPar.X;
        end
        if (isempty(InPar.Y))
            VecY = (1:1:Sy);
        else
            VecY = InPar.Y;
        end
        if (numel(OldVecX)==numel(VecX) && ...
            numel(OldVecY)==numel(VecY))
            % no need to regenerated MatX / MatY
        else
            [MatX,MatY] = meshgrid(VecX-Center(1),VecY-Center(2));
            OldVecX = VecX;
            OldVecY = VecY;
        end
        [OutX,OutY] = pol2cart(MatOutTheta,MatOutR);
        Sim(Isim).(InPar.ExecField{If}) = interp2(MatX,MatY,Sim(Isim).(InPar.ExecField{If}),OutX,OutY,InPar.InterpMethod);
      
    end
end

            
    
