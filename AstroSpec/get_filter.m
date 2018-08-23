function OutF=get_filter(FamilyNameStr,BandNameStr,varargin)
%--------------------------------------------------------------------------
% get_filter function                                            AstroSpec
% Description: Search and get astronomical Filter information and
%              transmission curve.
% Input  : - String (or cell array of strings) containing
%            filter family name.
%            Alternatively, this can be a two column matrix containing
%            the transmission curve.
%            If empty (i.e. []), donot search the family keyword.
%            Default is [].
%          - String (or cell array of strings) containing
%            band name.
%            If empty (i.e. []), donot search the band keyword.
%            If [] and FamilyNameStr is empty, then return
%            a list of all available families.
%            Default is [].
%          * Arbitrary number of pairs of arguments ...,key,val,...
%            The following keywords are available:
%            'Case ' - Case sensitive search {'y' | 'n'}, default is 'n'.
%            'Out'   - Output type: 'a' - (AstFilter) Astronomical filter
%                                         class.
%                                   's' - structure array.
%                                   'c' - cell array (default).
% Output : - Filter structure containing all the filters
%            found in database.
%            Output filter structure containing:
%                      F.family      - Family name
%                      F.band        - Band name
%                      F.T           - [Wave[Ang], Transm]
%                      F.nT          - [Wave[Ang], Norm. Transm]
%                      F.min_wl      - Min wavelength [Ang]
%                      F.max_wl      - Max wavelength [Ang]
%                      F.eff_wl      - Effective wavelength [Ang]
%                      F.half_width  - Filter half width [Ang]
% Tested : Matlab 7.0
%     By : Eran O. Ofek / Dovi Poznanski   Jan 2006
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Needed : Require Filters.mat containing the filters database.
% Example: OutF=get_filter('2MASS',[]);  % get all 2MASS filters
%           OutF=get_filter([],'u');      % get all u-band filters
%           OutF=get_filter('HST-ACS','F435W');  % get ACS F435W filter.
%           OutF=get_filter([],NaN);      % get all families
%           OutF=get_filter;              % get all families
%           OutF=get_filter([5000 0;5001 1;6000 1;6001 0],[],'out','a');
% Reliable: 1
%--------------------------------------------------------------------------

if (nargin==0),
   FamilyNameStr = [];
   BandNameStr   = [];
elseif (nargin==1),
   BandNameStr   = [];
else
   % do nothing
end


DefV.Case = 'n';
DefV.Out  = 'c';  % 'a' - for array; 'c' - for cell. 
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});


if (isnumeric(FamilyNameStr)),
    % Filter transmission curve is provided
    OutF.family{1}      = 'user';
    OutF.band{1}        = 'user';
    OutF.T{1}           = FamilyNameStr;
    OutF.nT{1}          = [FamilyNameStr(:,1), FamilyNameStr(:,2)./trapz(FamilyNameStr(:,1),FamilyNameStr(:,2))];
    OutF.min_wl{1}      = min(FamilyNameStr(:,1));
    OutF.max_wl{1}      = max(FamilyNameStr(:,1));
    OutF.eff_wl{1}      = sum(OutF.nT{1}(:,1).*OutF.nT{1}(:,2))./sum(OutF.nT{1}(:,2));
    li=min(find(FamilyNameStr(:,2)>0.2*max(FamilyNameStr(:,2))));
    ui=max(find(FamilyNameStr(:,2)>0.2*max(FamilyNameStr(:,2))));
    OutF.half_width{1}  =(FamilyNameStr(ui,1)-FamilyNameStr(li,1))./2;
    OutF.comments{1}    = '';
    OutF.source{1}      = '';
    OutF.Tunits{1}      = '';
else
    % Filter family is a char array
    
    load Filters.mat;

    if (isempty(FamilyNameStr)==1 && isnan(BandNameStr)==1),
       % retrieve all famalies
       Fam = 1;
       BandNameStr = [];
    else
       Fam = 0;
    end

    if (iscell(FamilyNameStr)==0),
       FamilyName{1} = FamilyNameStr;
    else
       FamilyName = FamilyNameStr;
    end
    if (iscell(BandNameStr)==0),
       BandName{1} = BandNameStr;
    else
       BandName = BandNameStr;
    end

    Nf = length(FamilyName);
    Nb = length(BandName);

    AllI = [];
    for If=1:1:Nf,
       for Ib=1:1:Nb,
          % search family
          if (isempty(FamilyName{If})==0),
             switch lower(InPar.Case)
                  case 'n'
                     FlagF = strcmpi(FamilyName{If},F.family);
                  case 'y'
                     FlagF = strcmp(FamilyName{If},F.family);
                  otherwise
                     error('Unkown caseSens Option');
             end
          else
             % all are flaged
             FlagF = ones(size(F.family));
          end
          % search band
          if (isempty(BandName{If})==0),
             switch lower(InPar.Case)
                  case 'n'
                     FlagB = strcmpi(BandName{If},F.band);
                  case 'y'
                     FlagB = strcmp(BandName{If},F.band);
                  otherwise
                     error('Unkown caseSens Option');
             end
          else
             % all are flaged
             FlagB = ones(size(F.family));
          end

          FoundI = find(FlagF==1 & FlagB==1);

          % indices of found filters
          AllI = [AllI; FoundI];
       end
    end

    switch Fam
     case 0
        OutF.family      = F.family(AllI);
        OutF.band        = F.band(AllI);
        OutF.T           = F.T(AllI);
        OutF.nT          = F.nT(AllI);
        OutF.min_wl      = F.min_wl(AllI);
        OutF.max_wl      = F.max_wl(AllI);
        OutF.eff_wl      = F.eff_wl(AllI);
        OutF.half_width  = F.half_width(AllI);
        OutF.comments    = F.comments(AllI);
        OutF.source      = F.source(AllI);
        OutF.Tunits      = F.Tunits(AllI);
     case 1
        [GroupStr]=group_cellstr(F.family,'y');
        OutF = GroupStr;
     otherwise
        error('Unknown Fam Option');
    end
end

%--- set output ---
switch lower(InPar.Out)
 case 'c'
    % do nothing
 case {'s','a'}
    N = length(OutF.T);
    switch lower(InPar.Out)
        case 's'
            OutA = struct('family',cell(N,1),'band',cell(N,1),...
                  'T',cell(N,1),'nT',cell(N,1),...
                  'min_wl',cell(N,1),'max_wl',cell(N,1),...
                  'eff_wl',cell(N,1),'half_width',cell(N,1));
        case 'a'
            OutA = astfilterdef(N,1);
    end
    
    for In=1:1:N,
       OutA(In).family     = OutF.family{In};
       OutA(In).band       = OutF.band{In};
       OutA(In).T          = OutF.T{In};
       OutA(In).nT         = OutF.nT{In};
       OutA(In).min_wl     = OutF.min_wl{In};
       OutA(In).max_wl     = OutF.max_wl{In};
       OutA(In).eff_wl     = OutF.eff_wl{In};
       OutA(In).half_width = OutF.half_width{In};
%OutF.comments{In}
       OutA(In).comments   = OutF.comments{In};
       OutA(In).source     = OutF.source{In};
       OutA(In).Tunits     = OutF.Tunits{In};
    end
    OutF = OutA;
    
 otherwise
    error('Unknown Out option');
end
