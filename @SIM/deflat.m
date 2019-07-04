function [Sim,FlatSim,NotIsFlat,Summary]=deflat(Sim,FlatSim,varargin)
% Correct images by flat field. Construct flat if needed.
% Package: @SIM
% Description: Correct images by flat field. Construct flat if needed.
% Input  : - A SIM object containing images to deflat.
%          - A SIM object containing flat images. SIM element per filter
%            group. If empty, then construct the Flat images using
%            SIM/flat.m. Default is empty.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'FilterKey' - A filter keyword string or a cell array of
%                          filter keyword strings by which to group the
%                          input images. A flat field image will be
%                          constructed for each group.
%                          Default is 'FILTER'.
%            'FlatPar'   - A cell array of parameters to pass to
%                          SIM/flat.m. Default is {}.
%            'IsFlatPar' - A cell array of parameters to pass to
%                          SIM/isflat.m. Default is {}.
%            'NotIsFlat' - A vector of logicals (or indices) indicating
%                          which input SIM images should be flat
%                          divided. If empty, then will attempt to
%                          identify the non-flat images and divide the
%                          flat from these images. If true, then divide
%                          the flat from all the images.
%                          Default is true.
% Output : - A SIM object with the flat-divided images.
%          - A SIM object with the flat images.
%          - vector of logicals for non-flat images.
%          - A summary structure (see SIM/flat.m for details).
%            Empty, if FlatSIm is provided by user.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Aug 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: S=deflat(S,FlatSim);
% Reliable: 2
%--------------------------------------------------------------------------


if (nargin==1)
    FlatSim = [];
end

DefV.FilterKey          = 'FILTER';
DefV.FlatPar            = {};
DefV.IsFlatPar          = {};
DefV.NotIsFlat          = true;    % list of images from which to subtract the bias
%DefV.ExecField          = {SIM.ImageField};

DefV.Verbose            = false;
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% construct super flat per filter group
IsFlat = [];
if (isempty(FlatSim))
    % User did not supply FlatSim
    % construct flat images
    IsFlat = isflat(Sim,InPar.IsFlatPar{:});
    [FlatSim,Summary] = flat(Sim,'IsFlat',IsFlat,InPar.FlatPar{:},'FilterKey',InPar.FilterKey);
else
    Summary = [];
end

% identify images from which to divide the superflat
if (isempty(InPar.NotIsFlat))
    if (isempty(IsFlat))
        % IsFlat was not populated
        IsFlat = isflat(Sim,InPar.IsFlatPar{:});
    end
    NotIsFlat = ~IsFlat;
elseif (numel(InPar.NotIsFlat)==1 && all(InPar.NotIsFlat))
    % divide flat from all the images
    NotIsFlat = true(size(Sim));
else
    % Assume NotIsFlat is in the correct form
    % containing: vector of logicals or indices
    NotIsFlat = InPar.NotIsFlat;
end


% for each Flat
Nfilter = numel(FlatSim);
FilterKeyVal = mgetkey(FlatSim,InPar.FilterKey);

for Ifilter=1:1:Nfilter
    KeyVal    = [InPar.FilterKey(:).',FilterKeyVal(Ifilter,1)];
    %GroupFlag = iskeyval(Sim(:),KeyVal(:)) & NotIsFlat(:);
    GroupFlag = iskeyval(Sim(:),KeyVal{:}) & NotIsFlat(:); % Na'ama, 2018-06-06
    
    if any(GroupFlag) % Na'ama, 2018-06-06
        % divide super flat
        % and combine the mask of the image with that of the flat
        Sim(GroupFlag) = bfun2sim(Sim(GroupFlag),FlatSim(Ifilter),@rdivide,'MaskFun',@mask_add,'MaskFunPar',{@bitor}); %,'ExecField',InPar.ExecField);

        if InPar.Verbose
            fprintf('%s: There are %d images in %s Filter deflat.\n',datestr(now,31),sum(GroupFlag),KeyVal{2})
        end
    end
end