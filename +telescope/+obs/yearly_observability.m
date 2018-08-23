function [DeltaHA_Hour]=yearly_observability(Year,ObjCoo,ObsCoo,TimeZone,AirMass,DeltaT)
% Plot yearly observability chart
% Pacakge: telescope.obs
% Description: Plot a yearly observability chart for an object.
% Input  : - Year (scalar).
%          - Matrix of objects coordinates: [[RA], [Dec]],
%            where [RA] is a column vector of R.A. in radians,
%            where [Dec] is a column vector of Dec. in radians.
%            If this matrix has 7 columns, then the 3 first
%            columns are taken as the R.A. [H M S],
%            and the next columns taken as the Dec. [Sign D M S].
%            One object per line.
%          - Observatory coordinates [East Long, North Lat] in radians,
%            or observatory name. For options of observatory names.:
%            see observatory_coo.m
%            Default is 'Wise'.
%          - East Time Zone [hours], default is 2.
%          - Minimum airmass lines [Airmass].
%            cuve for the times the object will cross (upper and lower),
%            the specified airmass.
%            If NaN, then no airmass line is plotted.
%            Default is NaN.
%          - Delta T (=TDT-UT1) in fraction of day.
%            Default is [1/1440].
% Output : - The maximum (abs. value) of the Hour Angle [hours], for which
%            the object is found in the needed AirMass.
% Plot   : Rise/set (red solid lines); 6/12/18/Twilight (red dotted lines);
%          midnight (green dotted); object crosses airmass (blue dashed)
% See also: daily_observability.m, obspl.m
% Tested : Matlab 5.3
%     By : Eran O. Ofek                    Sep 2001
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: telescope.obs.yearly_observability(2001,[10 0 0 -1 20 0 0],'Wise',2,2.5);
% Bug    : telescope.obs.yearly_observability(2001,[10 0 0 +1 20 0 0],'SAO',0,2.5);
% Reliable: 2
%------------------------------------------------------------------------------
RAD = 180./pi;

if (nargin==2)
   ObsCoo   = 'Wise';
   TimeZone = 2;
   AirMass  = NaN;
   DeltaT   = 1./1440;   
elseif (nargin==3)
   TimeZone = 2;
   AirMass  = NaN;
   DeltaT   = 1./1440;   
elseif (nargin==4)
   AirMass  = NaN;
   DeltaT   = 1./1440;   
elseif (nargin==5)
   DeltaT   = 1./1440;   
elseif (nargin==6)
   % do nothing
else
   error('Illegal number of input arguments');
end

if (ischar(ObsCoo)==1)
   [ObsCoo] = telescope.obs.observatory_coo(ObsCoo);
end

ColN = length(ObjCoo(1,:));
if (ColN==2)
   % already in radians
elseif (ColN==7)
   RA  = celestial.coo.convertdms(ObjCoo(:,1:3),'H','r');
   Dec = celestial.coo.convertdms(ObjCoo(:,4:7),'D','R');
   ObjCoo = [RA, Dec];
else
   error('ObjCoo should have 2 or 7 columns');
end


TZ_Offset = -1.*(ObsCoo(1).*RAD - TimeZone.*15)./360;  % fraction of day

Step   = 10;
Length = 366;
Margin = 0.02;

Nobj = length(ObjCoo(:,1));

% find H.A. observability zone
MinAlt = (pi./2 - celestial.coo.hardie_inv(AirMass));
CosH   = (sin(MinAlt) - sin(ObjCoo(:,2)).*sin(ObsCoo(2)))./(cos(ObjCoo(:,2)).*cos(ObsCoo(2)) );
DeltaH = acos(CosH);
DeltaHA_Hour = DeltaH.*RAD./15;

JD1    = celestial.time.julday([1 1 Year 0]);
JD     = (JD1:Step:JD1+Length).';
DayNum = (0:Step:Length-1).';
[Time,Ang]=celestial.SolarSys.sun_rise_set(JD,ObsCoo,TimeZone,DeltaT);

% calculate Equation of time
[RA,Dec,R,SL,EquationTime]=celestial.SolarSys.suncoo(JD,'a');
EquationTime = EquationTime./1440;  % convert to fraction of day

% calculate LST
LST = celestial.time.lst(JD,ObsCoo(1),'a');

% plot Sun rise/set/twilight
Offset = -1;
Time(:,6) = Offset + Time(:,6);
Time(:,7) = Offset + Time(:,7);
Time(:,8) = Offset + Time(:,8);
Time(:,9) = Offset + Time(:,9);

if ((Time(1,4)-Time(1,6))>1)
   Time(1,4) = Time(1,4) - 1;
end
if ((Time(1,3)-Time(1,7))>1)
   Time(1,3) = Time(1,3) - 1;
end
if ((Time(1,2)-Time(1,8))>1)
   Time(1,2) = Time(1,2) - 1;
end
if ((Time(1,1)-Time(1,9))>1)
   Time(1,1) = Time(1,1) - 1;
end


for Col=1:1:9
   Offset      = round(Time(2:end,Col)-Time(1,Col));
   Offset      = [0; Offset];
   Time(:,Col) = Time(:,Col) - Offset;
end

for Col=7:1:9
   I = find(abs(Time(:,6)-Time(:,Col))>0.5);
   Time(I,Col) = Time(I,Col) + sign(Time(I,6)-Time(I,Col));
end
for Col=1:1:3
   I = find(abs(Time(:,4)-Time(:,Col))>0.5);
   Time(I,Col) = Time(I,Col) + sign(Time(I,4)-Time(I,Col));
end

plot(Time(:,4),DayNum,'r-');
hold on;
plot(Time(:,1),DayNum,'r:');
plot(Time(:,2),DayNum,'r:');
plot(Time(:,3),DayNum,'r:');

plot(Time(:,6),DayNum,'r-');
plot(Time(:,7),DayNum,'r:');
plot(Time(:,8),DayNum,'r:');
plot(Time(:,9),DayNum,'r:');

%plot equation of time
plot(TZ_Offset+EquationTime,DayNum,'g:');

MinX = min(Time(:,6));
MaxX = max(Time(:,4));

Axis = axis;
axis([MinX-Margin, MaxX+Margin, 0, Length]);

PlotPosition = get(gcf,'Position');
set(gcf,'Position',[100 100 600 800]);

MinH = floor((MinX-Margin).*24);
MaxH = floor((MaxX+Margin).*24);

HourList    = [MinH:1:MaxH].';
FracList    = HourList./24;
I           = find(HourList<0);
HourList(I) = 24+HourList(I);

set(gca,'XTick',FracList);
TickN = length(HourList);
for I=1:1:TickN
   HourLabel(I,1:2) = sprintf('%2d',HourList(I));
end
set(gca,'XTickLabel',HourLabel);
set(gca,'FontSize',14);
if (TimeZone==0)
   h = xlabel('Universal Time');
else
   h = xlabel('Local Time');
end
set(h,'FontSize',18);

% delete right y-labels tick
set(gca,'YTickLabel','');
Pos = get(gca,'YTick');
for I=1:1:length(Pos)
   plot([MinX-Margin;MinX-Margin.*0.5],[Pos(I);Pos(I)],'w-');
   plot([MaxX+Margin.*0.5;MaxX+Margin],[Pos(I);Pos(I)],'w-');
end

% plot month ticks/labels
MonthName = strvcat('Jan','Feb','Mar','Apr',...
                    'May','Jun','Jul','Aug',...
                    'Sep','Oct','Nov','Dec');
for Month=1:1:12
   JDm   = celestial.time.julday([1 Month Year 0]);
   DayNm = JDm - JD1;
   plot([MaxX + Margin.*0.3; MaxX + Margin],[DayNm;DayNm],'k-');
   plot([MaxX + Margin.*0.6; MaxX + Margin],[DayNm+10;DayNm+10],'k-');
   plot([MaxX + Margin.*0.6; MaxX + Margin],[DayNm+20;DayNm+20],'k-');
   h=text(MaxX + Margin.*1.2, DayNm, MonthName(Month,1:3));
   set(h,'FontSize',12);
end


% plot objects transit

for I=1:1:Nobj
   %--- transit local time
   TLT = ObjCoo(I,1)./(2.*pi) - LST + TimeZone./24;      
   SortedTLT = sortrows([TLT, DayNum]);
   SortedTLT = [SortedTLT; [SortedTLT(:,1)+1, SortedTLT(:,2)]];
   Jdiff = find(diff(SortedTLT(:,2))>20);

   SortedTLT(Jdiff,1:2) = zeros(length(Jdiff),2).*NaN;
   plot(SortedTLT(:,1),SortedTLT(:,2),'b--');

   %--- Airmass lines
   if (isnan(AirMass)==1)
      % do nothing
   else
      %--- min H.A. local time
      TLT = ObjCoo(I,1)./(2.*pi) - LST + TimeZone./24 + DeltaH(I)./(2.*pi);      
      SortedTLT = sortrows([TLT, DayNum]);
      SortedTLT = [[SortedTLT(:,1)-1, SortedTLT(:,2)]; SortedTLT; [SortedTLT(:,1)+1, SortedTLT(:,2)]];
      Jdiff = find(diff(SortedTLT(:,2))>20);
      SortedTLT(Jdiff,1:2) = zeros(length(Jdiff),2).*NaN;
      plot(SortedTLT(:,1),SortedTLT(:,2),'b:');

      %--- max H.A. local time
      TLT = ObjCoo(I,1)./(2.*pi) - LST + TimeZone./24 - DeltaH(I)./(2.*pi);      
      SortedTLT = sortrows([TLT, DayNum]);
      SortedTLT = [[SortedTLT(:,1)-1, SortedTLT(:,2)]; SortedTLT; [SortedTLT(:,1)+1, SortedTLT(:,2)]];
      Jdiff = find(diff(SortedTLT(:,2))>20);
      SortedTLT(Jdiff,1:2) = zeros(length(Jdiff),2).*NaN;
      plot(SortedTLT(:,1),SortedTLT(:,2),'b:');
   end
end




