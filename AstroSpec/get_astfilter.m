function AstF=get_astfilter(varargin)
%--------------------------------------------------------------------------
% get_astfilter function                                         AstroSpec
% Description: Get astronomical filters in AstFilter class format using
%              exact name search on family or band names.
%              This function replaces get_filter.m.
%              OBSOLETE: Ues AstFilter.get instead.
% Input  : - Family name or a cell array of family names.
%            If empty then search only by band name.
%            Note that this function have a different behavior if
%            the first argument is of AstFilter class.
%          - Band name or a cell array of band names.
%            If empty then search only by family name.
%          - Behaviour options:
%            'i' - case insenstive search (default).
%            'c' - case sensitive search.
% Output : - Astronomical filter class structure array of selected
%            filters.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstF=get_astfilter('SDSS','r');
%          AstF=get_astfilter([],'r');
%          AstF=get_astfilter;
%          AstF=get_filter({'PTF','sdss'});
% Reliable: 2
%--------------------------------------------------------------------------

LoadCat = false;
if (numel(varargin)>0)
    if (isastfilter(varargin{1})),
        AstFilterCat = varargin{1};
        if (numel(varargin)>1),
            varargin = varargin(2:end);
        else
            varargin = {};
        end
    else
        LoadCat = true;
    end
else
    LoadCat = true;
end

if (LoadCat),
    AstFilterCat = load_check('AstFilterCat.mat');
end

Def.Family = [];
Def.Band   = [];
Def.Behav  = 'i'; 

Narg = numel(varargin);
if (Narg==0),
    Family = Def.Family;
    Band   = Def.Band;
    Behav  = Def.Behav;
elseif (Narg==1),
    Family = varargin{1};
    Band   = Def.Band;
    Behav  = Def.Behav;
elseif (Narg==2),
    Family = varargin{1};
    Band   = varargin{2};
    Behav  = Def.Behav;
elseif (Narg==3),
    Family = varargin{1};
    Band   = varargin{2};
    Behav  = varargin{3};
    
else
    error('Illegal number of input arguments');
end

switch lower(Behav)
    case 'i'
        % case insensetive
        BehavFun = @strcmpi;
    case 'c'
        BehavFun = @strcmp;
    otherwise
        error('Unknown Behav option');
     
end

Nf           = numel(AstFilterCat);

if (isempty(Family) && isempty(Band)),
    % get all filters
    IndF = (1:1:Nf);    
else
    if (isempty(Band)),
        % Family is empty
        % Band is given
        if (iscell(Family)),
            Nfam = numel(Family);
            IndF = [];
            for Ifam=1:1:Nfam,
                IndF = [IndF, find(BehavFun({AstFilterCat.family},Family{Ifam}))];
            end
        else
            IndF = find(BehavFun({AstFilterCat.family},Family));
        end
    elseif (isempty(Family)),
        % Band is empty 
        % Family is given
        if (iscell(Band)),
            Nband = numel(Band);
            IndF = [];
            for Iband=1:1:Nband,
                IndF = [IndF, find(BehavFun({AstFilterCat.band},Band{Iband}))];
            end
        else
            IndF = find(BehavFun({AstFilterCat.band},Band));
        end
    else
        % both band and family are provided
        IndF = find(BehavFun({AstFilterCat.band},Band) & ...
                    BehavFun({AstFilterCat.family},Family));
    end
end
AstF = AstFilterCat(IndF);
        
            
        
        
