function [Sim,Gain]=gain_correct(Sim,varargin)
% Correct images in SIM for gain (i.e., set gain to 1).
% Package: @SIM
% Description: Correct images in SIM for gain (i.e., set gain to 1).
% Input  : - A SIM object with images and catalog.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField'     - SIM fields on which to execute the gain
%                              correction. Default is {'Im'}.
%            'ExecFieldSqrt' - SIM fields on which to execute the
%                              sqrt(gain)correction. Default is {}.
%            'OrigGainKey'   - An header keyword name in which to store
%                              the original gain value.
%                              Default is 'ORIGGAIN'.
%            'CatColUpdtae'  - A cell array of catalog column names for
%                              which to apply the gain factor.
%                              Default is {}.
%            'GainKeys'      - Cell array of GAIN keywords in which to look
%                              for the gain . Alternatively this can be a
%                              vector of numbers (gain per image).
%                              Defauly is {'GAIN'}.
%            'SelectionMethod'- Gain key val selection method.
%                              See getkey_gain.m for options.
%                              Defauly is 'first'.
% Output : - A SIM object with the images and catalog corrected for the
%            gain factor.
%          - Vector of the original gain factors.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jun 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=gain_correct(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

ImageField         = 'Im';


DefV.ExecField          = {ImageField};
DefV.ExecFieldSqrt      = {};
DefV.OrigGainKey        = 'ORIGGAIN';
DefV.CatColUpdate       = {};
% getkey_gain.m arguments
DefV.GainKeys           = {'GAIN'};
DefV.SelectionMethod    = 'first';
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% Get gain of all images
if (isnumeric(InPar.GainKeys))
    Gain = InPar.GainKeys;
    InPar.GainKeys = {'GAIN'};
else
    Gain = getkey_gain(Sim,'GainKeys',InPar.GainKeys,'SelectionMethod',InPar.SelectionMethod);
end

% Deal with Gain==NaN
if (any(isnan(Gain)))
    warning('Gain==NaN identified - set gain to 1');
    Gain(isnan(Gain)) = 1;
end

% Deal with Gain==0
if (any(Gain==0))
    warning('Gain==0 identified - set gain to 1');
    Gain(Gain==0) = 1;
end


% Correct for gain factor
if (any(Gain~=1))
    % No need to run if Gain~=1
    Sim = bfun2sim(Sim,Gain,@times,'ExecField',InPar.ExecField);
    % Apply sqrt(Gain) for specific fields
    if (~isempty(InPar.ExecFieldSqrt))
        Sim = bfun2sim(Sim,sqrt(Gain),@times,'ExecField',InPar.ExecFieldSqrt);
    end
end


% Update header keywords + Catalog
Nsim = numel(Sim);
for Isim=1:1:Nsim
    % Original gain
    Sim(Isim) = replace_key(Sim(Isim),InPar.OrigGainKey,Gain(Isim),'Original GAIN [e-/ADU]');
    % New gain
    Sim(Isim) = replace_key(Sim(Isim),InPar.GainKeys{1},1,'GAIN (after setting to 1) [e-/ADU]');
    
    % Update catalog
    if (~isempty(InPar.CatColUpdate))
        Sim = col_fun(Sim(Isim),@times,InPar.CatColUpdtae,[],Gain(Isim));
    end
end



