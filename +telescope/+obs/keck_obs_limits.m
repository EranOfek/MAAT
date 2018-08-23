function [RTS,AM]=keck_obs_limits(Obs,Date,RA,Dec)
% Rise/set time for object in Keck observatory given telescope limits.
% Package: telescope.obs
% Description: Given a date and object celestial positions, calculate the
%              rise and set time of an object in Keck observatory, given
%              the Nasmyth mount limits.
% Input  : - Observatory: 
%            {'KeckI' | 'KeckII'}
%          - Date [D M Y] or [JD].
%          - List of target's RA [H M S], ['HH:MM:SS.S'] or [rad].
%            see convertdms.m for more details.
%          - List of target's Dec [Sign D M S], ['+DD:MM:SS.S'] or [rad].
%            see convertdms.m for more details.
% Output : - Crossing the visibility limits [Rise, Transit, Set] in UT
%            fraction of day. 
%          - Airmass at [Rise, Transit, Set].
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Oct 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [RTS,AM]=keck_obs_limits('KeckI',[1 7 2011],[18 00 00],[-1 30 0 0])
% Reliable: 2
%------------------------------------------------------------------------------

RAD         = 180./pi;
ColAz       = 1;
ColAlt      = 2; 
KeckI_Ind   = 1;
KeckII_Ind  = 2;
GeoPos(KeckI_Ind,:)    = [-155.478 19.828]./RAD;
GeoPos(KeckII_Ind,:)   = [-155.478 19.828]./RAD;
Limits{KeckI_Ind}    = [  0     18
                          5     18
                          5.001 33
                        145     33
                        145.001 18
                        360     18]./RAD;
Limits{KeckII_Ind}   = [  0     18
                        185     18
                        185.001 33
                        335     33
                        335.001 18
                        360     18]./RAD;

switch Obs
 case 'KeckI'
    ObsInd = KeckI_Ind;
 case 'KeckII'
    ObsInd = KeckII_Ind;
 otherwise
    error('Unknown Obs Option');
end


 
% Date  
if (size(Date,2)==3),
   JD = julday([Date, 0]);
else
   JD = Date;
end

% RA
if (isstr(RA)==1),
   RA = convertdms(RA,'SH','r');
else
   if (size(RA,2)==3),
      RA = convertdms(RA,'H','r');
   elseif (size(RA,2)==1),
      RA = RA;
   else
      error('Illegal RA format');
   end
end

% Dec
if (isstr(Dec)==1),
   Dec = convertdms(Dec,'SD','R');
else
   if (size(Dec,2)==4),
      Dec = convertdms(Dec,'D','R');
   elseif (size(Dec,2)==1),
      Dec = Dec;
   else
      error('Illegal Dec format');
   end
end


Nobj = length(RA);
FracDay = [0:1./1440:1].';


for Iobj=1:1:Nobj,
   HorCoo = horiz_coo([RA(Iobj), Dec(Iobj)],JD+FracDay,GeoPos(ObsInd,:),'h');

   AltLimit = interp1(Limits{ObsInd}(:,ColAz), Limits{ObsInd}(:,ColAlt), HorCoo(:,ColAz),'linear');
   FlagOn   = HorCoo(:,ColAlt) > AltLimit;
%[HorCoo(:,ColAz).*RAD, HorCoo(:,ColAlt).*RAD, AltLimit.*RAD]
   RiseInd  = find(diff(FlagOn)>0)+1;
   SetInd   = find(diff(FlagOn)<0)+1;
   [MaxAlt,MaxAltInd] = max(HorCoo(:,ColAlt));
   RiseAlt  = HorCoo(RiseInd,ColAlt);
   TranAlt  = HorCoo(MaxAltInd,ColAlt);
   SetAlt   = HorCoo(SetInd,ColAlt);
   RTS(Iobj,:) = [FracDay(RiseInd), FracDay(MaxAltInd), FracDay(SetInd)];

   AM(Iobj,:)  = [hardie(pi./2-RiseAlt), hardie(pi./2-TranAlt), hardie(pi./2-SetAlt)];

end

