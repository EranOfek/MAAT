function Sim=sub_background(Sim,varargin)
% Subtract background from SIM images (see aslo background).
%--------------------------------------------------------------------------
% sub_background function                                       class/@SIM
% Description: Given a SIM object array of images in which the background
%              field ('BackIm') is populated, subtract the background from
%              the image in each SIM element.
%              For background estimation and subtraction see also
%              background.m.
% Input  : - A SIM object in which the background image field is populated.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ReplaceKey'- A three column cell array of {key,val,comment},
%                          or an HEAD object which keywords to replace
%                          or add to the SIM header.
%                          Default is {}.
%            'AddKey'    - Like 'ReplaceKey', but adding keywords,
%                          without replacment.
%                          Default is {}.
% Output : - A background subtracted SIM object.
% See also: SIM/background.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sb=sub_background(S);
% Reliable: 2
%--------------------------------------------------------------------------

ImageField    = 'Im';
BackField     = 'BackIm';


DefV.AddKey             = {};
DefV.ReplaceKey         = {};
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Nsim = numel(Sim);
for Isim=1:1:Nsim
    % for each SIM element
    %Sim(Isim).(ImageField) = Sim(Isim).(ImageField) - Sim(Isim).(BackField);
    Sim(Isim).(ImageField) = double(Sim(Isim).(ImageField)) - Sim(Isim).(BackField); % Na'ama, 20180828
    
    % Add/modify keywords in header
    if (~isempty(InPar.ReplaceKey))
        Sim = replace_key(Sim,InPar.ReplaceKey);
    end
    if (~isempty(InPar.AddKey))
        Sim = add_key(Sim,InPar.AddKey);
    end
    
end
