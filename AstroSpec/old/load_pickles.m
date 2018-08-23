function AstS=load_pickles(SpC,SpL)
%--------------------------------------------------------------------------
% load_pickles function                                          AstroSpec
% Description: Load the Pickles stellar spectra library into a AstSpec
%              class object. By default will load all the normal stars.
% Input  : - Optional spectral type to load (e.g., 'G'). If empty, load
%            all. Default is empty.
%          - Optional luminosity class to load (e.g., 'V'). If empty, load
%            all. Default is empty.
% Output : - AstSpec class object will all the spectra.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Feb 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstS=load_pickles;
%          AstS=load_pickles('g','v');
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==0),
    SpC = [];
    SpL = [];
elseif (nargin==1),
    SpL = [];
else
    % do nothing
end
    

DirPic = '~/matlab/data/PicklesStellarSpec';
PWD    = pwd;
cd(DirPic);

[~,List] = create_list('uk*.mat',NaN);
RE  = regexp(List,'uk(?<SpClass>[obafgkm]+)(?<SpNum>\d+)(?<SpLum>\w+)','names');
Nre = numel(RE);

Ind = 0;
for Ire=1:1:Nre,
    if (~isempty(RE{Ire})),
        Ind = Ind + 1;
        SRE(Ind) = RE{Ire};
        ListN{Ind} = List{Ire};
    end
end

Nre = numel(SRE);
if (isempty(SpC)),
    FlagC = true(1,Nre);
else
    FlagC = strcmpi({SRE.SpClass},SpC);
end
if (isempty(SpL)),
    FlagL = true(1,Nre);
else
    FlagL = strcmpi({SRE.SpLum},SpL);
end

Flag = FlagC & FlagL;


%AstS = spec_read_mat('uk*.mat');
AstS = spec_read_mat(ListN(Flag));
SRE  = SRE(Flag);

Ns   = numel(AstS);
for Is=1:1:Ns,
    AstS(Is).z = 0;
    AstS(Is).source = 'Pickles 1998';
    
    if (numel(SRE(Is).SpNum)>1),
        SpNum = str2double(SRE(Is).SpNum(1))+0.5;
    else
        SpNum = str2double(SRE(Is).SpNum);
    end
    AstS(Is).ObjName = sprintf('%s %3.1f %s',upper(SRE(Is).SpClass),SpNum,upper(SRE(Is).SpLum));
    
end

cd(PWD);