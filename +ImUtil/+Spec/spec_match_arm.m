function [Match,Out]=spec_match_arm(InfoS,varargin)
%--------------------------------------------------------------------------
% spec_match_arm function                                           ImSpec
% Description: Given the InfoS structure returned by spec_classify_images.m
%              try to group the observations into objects (e.g., images
%              taken at the same coordinates).
%              Next, let the use edit the groups and the imnages to use.
% Input  : - InfoS structure returned by spec_classify_images.m
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'MatchDist'  Matching distance [arcsec]. Default is 30.
% Output : - 
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: InfoS = spec_classify_images('*.fits');
%          Match=spec_match_arm(InfoS);
% Reliable: 
%--------------------------------------------------------------------------
RAD     = 180./pi;
SEC_DAY = 86400;

%global SI

%DefV.AddFields = {};
DefV.MatchDist = 30;   % [arcsec]
InPar          = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar.MatchDist = InPar.MatchDist./(3600.*RAD);   % [radians]

UniqueArm = unique([InfoS.ArmIndex].');
Nua       = numel(UniqueArm);

for Iua=1:1:Nua,
    IndImArm    = find([InfoS.ArmIndex].'==UniqueArm(Iua));
    IndOtherArm = find([InfoS.ArmIndex].'~=UniqueArm(Iua));
    OtherStartTime = [InfoS(IndOtherArm).JD].';
    OtherEndTime   = OtherStartTime + [InfoS(IndOtherArm).ExpTime].'./SEC_DAY;
    Nim      = numel(IndImArm);
    for Iim=1:1:Nim,
        % for each image in the current arm configuration
        % look for images taken at the same time window
        TimeWin = [[InfoS(IndImArm(Iim)).JD].', [InfoS(IndImArm(Iim)).JD].'+[InfoS(IndImArm(Iim)).ExpTime./SEC_DAY].'];
        
        % look for time intersecting observations
        Flag = (OtherStartTime>TimeWin(1) & OtherEndTime<TimeWin(2)) | ...
               (OtherStartTime<TimeWin(1) & OtherEndTime>TimeWin(1)) | ...
               (OtherStartTime<TimeWin(2) & OtherEndTime>TimeWin(2));
        Match(IndImArm(Iim)).IndByTime = IndOtherArm(Flag);
        
        % look for observations with similar coordinates 
        % in all arms
        Dist = sphere_dist(InfoS(IndImArm(Iim)).RA, InfoS(IndImArm(Iim)).Dec, ...
                           [InfoS.RA].', [InfoS.Dec].','deg');
        
        Match(IndImArm(Iim)).IndByCoo = find(Dist<InPar.MatchDist);
    end
end
% group by coordinates and no gaps in time
[~,SI] = sort([InfoS.JD].');
Group = zeros(1,Nim).*NaN;
GroupInd = 0;
for Iim=1:1:Nim,
    switch lower(InfoS(Iim).ImType)
        case 'object'
           if (~isempty(Match(Iim).IndByCoo)),
               if (isnan(Group(Iim))),

                  GroupInd = GroupInd + 1;
                  %Group(Iim) = GroupInd;
                  Group(Match(Iim).IndByCoo) = GroupInd;
               end
           end
        otherwise
            % do nothing
    end
end

       
   



% edit log of observations
Nim    = numel(InfoS);
UseIm  = num2cell(ones(1,Nim));
GroupIm= num2cell(Group);
OrigInd= num2cell((1:1:Nim));
Columns= {'File','Object','ImType','Date','ExpTime','ArmInd','RA','Dec','Dichroic','SlitW','Turret','OrigInd','UseIm','GroupIm'};
InfoS1 = struct('ImageFileName',{InfoS.ImageFileName},...
                'Object',{InfoS.Object},...
                'ImType',{InfoS.ImType},...
                'Date',{InfoS.Date},...
                'ExpTime',{InfoS.ExpTime},...
                'ArmIndex',{InfoS.ArmIndex},...
                'RA',{InfoS.RA},...
                'Dec',{InfoS.Dec},...
                'Dichroic',{InfoS.Dichroic},...
                'SlitWidth',{InfoS.SlitWidth},...
                'Turret',{InfoS.Turret},...
                'OrigInd',OrigInd,...
                'UseIm',UseIm,...
                'GroupIm',GroupIm);
            
C=squeeze(struct2cell(InfoS1)).';
% sort by time
[~,SI] = sort([InfoS.JD].');

TmpC = C(SI,:);
TmpC = md_edit(TmpC,'label',{'123',Columns});

Out.ColCell = Columns;
Out.Col     = cell2struct(num2cell(1:1:length(Columns)),Out.ColCell,2);
Out.C       = sortrows(TmpC,Out.Col.OrigInd);  %TmpC; %(SI,:);

close;

%assignin('base','TmpC',TmpC)
%assignin('base','SI',SI)
%openvar('TmpC');

% drawnow;
% pause(0.5); % wait for client to become visible
% 
% % Get handle of variables client in the Variable Editor
% %http://blogs.mathworks.com/community/2008/04/21/variable-editor/
% jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
% jClient = jDesktop.getClient('TmpC');
% hjClient = handle(jClient,'CallbackProperties');
% 
% % Instrument the client to fire callback when it's closed
% set(hjClient,'ComponentRemovedCallback',{@my_callback_epar,'TmpC'});
% 
% %C=evalin('base','TmpC');
% %R = input('Press any key when finished editing','s');
% % need to wait...




% function my_callback_epar(varEditorObj,eventData,varname)
% % do nothing
% %global FunParsFullName
% %global SI
% 
% 
% C = TmpC;
% OrigInd = [C{:,end-2}];
% IsUse = [C{:,end-1}];
% %IsUse = IsUse(SI);
% UserGroup = [C{:,end}];
% %UserGroup = UserGroup(SI);
% CellOrigInd = num2cell(OrigInd);
% [InfoS(OrigInd).OrigInd] = deal(CellOrigInd{:});
% CellIsUse   = num2cell(IsUse);
% [InfoS(OrigInd).IsUse] = deal(CellIsUse{:});
% CellUserGroup   = num2cell(UserGroup);
% [InfoS(OrigInd).UserGroup] = deal(CellUserGroup{:});
% 
% 
% % an example for all the spectra of one source
% %InfoS([InfoS.UserGroup]==4).Object
% 
% %save(FunParsFullName,'IsUse','UserGroup');
% %TmpC = evalin('base','TmpC');
% %IsUse = evalin('base','IsUse');
% %UserGroup = evalin('base','UserGroup');
% 
% 
