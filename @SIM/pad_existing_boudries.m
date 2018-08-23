function Sim=pad_existing_boudries(Sim,Width,PadVal,ExecField)
%--------------------------------------------------------------------------
% pad_existing_boudries function                                class/@SIM
% Description: Set the edges of images in a SIM object to some value.
% Input  : - A SIM object.
%          - Distance from boundries to pad. Default is 5.
%          - Pad value. Default is 0.
%          - SIM fields to pad. Default is {'Im'}.
% Output : - A padded SIM object.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=pad_existing_boudries(Sim,5,NaN);
% Reliable: 2
%--------------------------------------------------------------------------

ImageField = 'Im';

Def.Width     = 5;
Def.PadVal    = 0;
Def.ExecField = {ImageField};

if (nargin==1),
    Width     = Def.Width;
    PadVal    = Def.PadVal;
    ExecField = Def.ExecField;
elseif (nargin==2),
    PadVal    = Def.PadVal;
    ExecField = Def.ExecField;
elseif (nargin==3),
    ExecField = Def.ExecField;
elseif (nargin==4),
    % do nothing
else
    error('Illigal number of input arguments: Sim, [Width, PadVal, Fields]');
end


Nsim = numel(Sim);
if (~iscell(ExecField)),
    ExecField = {ExecField};
end
Nf   = numel(ExecField);

for Isim=1:1:Nsim,
    for If=1:1:Nf,
        Sim(Isim).(ExecField{If})(1:Width,:)         = PadVal;
        Sim(Isim).(ExecField{If})(end-Width+1:end,:) = PadVal;
        Sim(Isim).(ExecField{If})(:,1:Width)         = PadVal;
        Sim(Isim).(ExecField{If})(:,end-Width+1:end) = PadVal;
    end
end
