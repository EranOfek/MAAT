function [RA,Dec]=get_simbad_coo(ObjectName);
%-------------------------------------------------------------------------
% get_simbad_coo function                                       Catalogue
% Description: Convert object name to coordinates using
%              the SIMBAD name server.
% Input  : - String or cell vector of object names.
% Output : - RA [radians].
%          - Dec [radians].
% Tested : Matlab 7.0
%     By : Eran O. Ofek                       May 2006
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% Example: [RA,Dec]=get_simbad_coo({'Z Cam';'Sirius';'beta UMa'});
%-------------------------------------------------------------------------
URL_Part1 = 'http://simbad.u-strasbg.fr/sim-id.pl?protocol=html&Ident=';

if (iscell(ObjectName)==0),
   ObjectList{1} = ObjectName;
else
   ObjectList    = ObjectName;
end

N   = length(ObjectList);
RA  = zeros(N,1).*NaN;
Dec = zeros(N,1).*NaN;
for I=1:1:N,
   % replace spaces in name by '+'
   IndSpace = strfind(ObjectList{I},' ');
   ObjectList{I}(IndSpace) = '+';

   StrSIMBAD = urlread(sprintf('%s%s',URL_Part1,ObjectList{I}));
   J = strfind(StrSIMBAD,'ICRS 2000.0 coordinates');
   SubStr = StrSIMBAD(J+52:J+78);
   if (isempty(SubStr)==0),
      [RA_H,RA_M,RA_S,Dec_D,Dec_M,Dec_S,Junk1,Junk2,Junk3]=strread(SubStr,'%d %d %f %s %d %f %s %s %s');
      Sign = Dec_D{1}(1:1);
      Dec_D = str2num(Dec_D{1}(2:end));
      switch Sign
       case '+'
          DecSign = +1;
       case '-'
          DecSign = -1;
       otherwise
          error('unable to read declination sign');
      end
      RA(I)  = convertdms([RA_H RA_M RA_S],'H','r');
      Dec(I) = convertdms([DecSign Dec_D Dec_M Dec_S],'D','R');
   end
end
