function add_specphot_stand(FileName,Factor,Name,RA,Dec,MagV,SpecType,Comment)
%------------------------------------------------------------------------------
% add_specphot_stand function                                        AstroSpec
% Description: Add a spectrophotometric standard to list of standard stars.
% Input  : - File name containing the standard star spectrum, with the
%            following columns: 
%            Wavelength [A]
%            Flux [erg/cm/cm/s/A * 10^16].
%          - Factor to multiply flux in order to get the correct units.
%          - String or cell array of strings containing the star names.
%            The star may have more than one name.
%          - J2000.0 R.A., in radians, [H M S] or sexagesimal string.
%          - J2000.0 Dec., in radians, [Sign D M S] or sexagesimal string.
%          - V-band magnitude.
%          - String containing spectral type of star.
%          - Comments. Default is ''.
% Output : null
% Result : Updating the SpecPhot_Stand.mat file
% Tested : Matlab 7.16
%     By : Eran O. Ofek                    Sep 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: search_specphot_stand.m
% Example: add_specphot_stand('fhr4554.dat',1e-16,{'HR4554','HD103287'},[11 53 49.84732],[+1 53 41 41.1350],2.44,'A0Ve');
% Reliable: 2
%------------------------------------------------------------------------------
Dir = Util.files.which_dir('SpecPhot_Stand.mat');
PWD = pwd;
IndExist = [];

if (nargin==7),
   Comment = '';
end


if (isempty(Dir)),
   Dir = PWD;
   I   = 1;
else
   SpecPhot_Stand = Util.IO.load2('SpecPhot_Stand.mat');
   I = length(SpecPhot_Stand);
   I = I + 1;

   for Is=1:1:I-1,
       if (sum((strcmp(lower(SpecPhot_Stand(Is).Name),lower(Name{1})))==0)>0),
 	  IndExist = Is;
      end
   end   


end

if (ischar(Name)),
   Name = {Name};
end

Spec = load2(FileName);
Spec(:,2) = Spec(:,2).*Factor;


if (isempty(IndExist)==0),
   % star name exist in list
   Ans = input('Star name is already exist - do you want to update entry? [Y/N]','s');

   switch lower(Ans)
    case 'y'
       I = IndExist;
    case 'n'
       error('Star name is already exist');
    otherwise
       error('Unknown Ans option');
   end
end


SpecPhot_Stand(I).Name     = Name;
SpecPhot_Stand(I).RA       = celestial.coo.convertdms(RA,'gH','r');
SpecPhot_Stand(I).Dec      = celestial.coo.convertdms(Dec,'gD','r');
SpecPhot_Stand(I).MagV     = MagV;
SpecPhot_Stand(I).SpecType = SpecType;
SpecPhot_Stand(I).Spec     = Spec;
SpecPhot_Stand(I).Comment  = Comment;

% setup grid points
R = input('Do you want to select grid points? [Y/N]');
if (strcmpi(R,'y')),
    load Telluric.mat;
    clf;
    semilogy(SpecPhot_Stand(I).Spec(:,1),SpecPhot_Stand(I).Spec(:,2),'k-');
    hold on;
    
    FlagT = Util.array.find_ranges_flag(SpecPhot_Stand(I).Spec(:,1),Telluric.Cat(:,1:2));
   
    TelS = SpecPhot_Stand(I).Spec(:,2).*(FlagT>1);
    TelS(TelS==0) = NaN;
    plot(SpecPhot_Stand(I).Spec(:,1),TelS,'r-','LineWidth',2);
    
    Points=plot_select_points(SpecPhot_Stand(I).Spec(:,1),SpecPhot_Stand(I).Spec(:,2));
    Selected=[[Points.X]',[Points.Y]'];
    [~,UI] = unique(Selected(:,1));
    Selected = Selected(UI,:);
    %plot(Selected(:,1),Selected(:,2),'bx')
    
    %Res=plot_int({SpecPhot_Stand(I).Spec(:,1),SpecPhot_Stand(I).Spec(:,2),'-',{},'X','Y'});
    %waitfor(gcf,'KeyPressFcn','');

    SpecPhot_Stand(I).GridWave=Selected(:,1);  %SpecPhot_Stand(I).Spec(Res.IndRM,1);
end


cd(Dir);
save SpecPhot_Stand.mat SpecPhot_Stand
cd(PWD);
