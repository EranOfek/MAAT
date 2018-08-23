function varargout=sim2cube(Sim,varargin)
% Convert the images in a SIM object array into a cube.
% Package: @SIM
% Description: Convert the images in a SIM object array into a cube.
% Input  : - A SIM object array.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField'  - A cell array (or a single string) of field
%                           names which contains images that should be
%                           converted into a cube.
%                           Default is {'Im'}.
%            'ImDim'      - The dimension in which to store the image
%                           Index (1|3). Default is 1.
% Output : * The number of output arguments is equal to the number of
%            strings in the ExecField input. Each output is the image cube
%            corresponding to the SIM object field in the ExecField
%            input.
% See also: cube2sim.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: C=sim2cube(S,'ExecField',{'Im','BackIm'});
%          C=sim2cube(S);
% Reliable: 2
%--------------------------------------------------------------------------



DefV.ExecField          = {'Im'};
DefV.CCDSEC             = [];
DefV.ImDim              = 1;
if (numel(varargin)>0)
    %InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
    InPar = InArg.populate_keyval(DefV,varargin,mfilename);
else
    InPar = DefV;
end

if (~iscell(InPar.ExecField))
    InPar.ExecField = {InPar.ExecField};
end
Nf = numel(InPar.ExecField);

Nsim = numel(Sim);
varargout = cell(1,Nf);
for If=1:1:Nf
    
    Size1 = size(Sim(1).(InPar.ExecField{If}));
    % init
    varargout{If} = zeros(Size1(1),Size1(2),Nsim);
    % copy SIM into cube
    for Isim=1:1:Nsim
        varargout{If}(:,:,Isim) = Sim(Isim).(InPar.ExecField{If});
    end
        
    if (InPar.ImDim==3)
        % do nothing
        
    elseif (InPar.ImDim==1)
        varargout{If} = shiftdim(varargout{If},2);
        
%         varargout{If} = zeros(Nsim,Size1(1),Size1(2));
%         for Isim=1:1:Nsim,
%             varargout{If}(Isim,:,:) = Sim(Isim).(InPar.ExecField{If});
%         end
        
    else
        error('Unknown ImDim option');
    end
end


