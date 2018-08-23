function SN=update_sn_cat;
%---------------------------------------------------------------------------
% update_sn_cat function                                          Catalogue
% Description: Update the SN catalog from
%              http://cfa-www.harvard.edu/iau/lists/Supernovae.html
% Input  : - NULL
% Output : - Structure containing SNe catalog.
% Example: SN=update_sn_cat;
%          SN.update        - date at which structure was produced
%          SN.Name
%          SN.Host
%          SN.DiscJD
%          SN.Mag
%          SN.RA
%          SN.Dec
%          SN.DiscRef
%          SN.PosRef
%          SN.Type
%          SN.Discoverer
% Tested : Matlab 7.0
%     By : Eran O. Ofek
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%---------------------------------------------------------------------------
RAD   = 180./pi;

Nhead = 27;

TempFile = 'Tmp.SN.url';
%system(sprintf('lynx -dump -dont_wrap_pre -nolist http://cfa-www.harvard.edu/iau/lists/Supernovae.html > %s',TempFile));

Cell = file2str(TempFile,'cell');

%delete(TempFile);

Nc   = length(Cell);
Isn  = 0;
for I=Nhead+1:1:Nc-7,
I
   if (length(Cell{I}>1)),
      Isn = Isn + 1;
      SN.Name{Isn}    = Cell{I}(1:8);
      SN.Host{Isn}    = Cell{I}(9:25);

      Year   = sscanf(Cell{I}(26:29),'%d');
      Month  = sscanf(Cell{I}(31:32),'%d');
      Day    = sscanf(Cell{I}(34:35),'%d');
      if (isstr(Year)==1),
         % not a SN
 	 SN.DiscJD(Isn) = NaN;
         SN.Mag(Isn)    = NaN;
      else
         Mag = sscanf(Cell{I}(64:70),'%f');
         if (strcmp(Mag,'')==1),
 	    SN.Mag(Isn)  = NaN;
         else
	    SN.Mag(Isn)  = Mag;
         end
         if (isstr(Month)==1),
 	    Month = 0;
         end
         if (isstr(Day)==1),
 	    Day   = 0;
         end
         SN.DiscJD(Isn) = julday([Day Month Year]);
      end   
      RA     = sscanf(Cell{I}(88:99),'%f %f %f');
      if (strcmp(RA,'')==1),
 	 % alternate calculation of RA/Dec
	 [RA, Dec] = get_simbad_coo(SN.Host{Isn});

         DelRA  = sscanf(Cell{I}(54:56),'%f');
         DelDec = sscanf(Cell{I}(59:61),'%f');
         if (strcmp(DelRA,'')==1 | strcmp(DelDec,'')==1),
            DelRA  = 0;
            DelDec = 0;
         else
            switch Cell{I}(57)
             case 'W'
                DelRA = -DelRA;
             case 'E'
                DelRA = DelRA;
             otherwise
                % do nothing: error('Unknonw DelRA sign');
            end
            switch Cell{I}(62)
             case 'N'
                DelDec = DelDec;
             case 'S'
                DelDec = -DelDec;
             otherwise
                % do nothing: error('Unknonw DelDec sign');
            end
         end

	 SN.RA(Isn) = RA + (DelRA./cos(Dec))./(3600.*RAD);
         SN.Dec(Isn) = Dec + DelDec./(3600.*RAD);


      else
         Dec    = sscanf(Cell{I}(101:112),'%f %f %f');
         DecSign = Cell{I}(100);
         switch DecSign
          case '+'
             Sign = 1;
          case '-'
             Sign = -1;
          otherwise
             error('Bad sign in SN declination');
         end
         if (length(RA)==2),
            SN.RA(Isn) = convertdms([RA(1), RA(2)],'HM','r');
         elseif (length(RA)==3),
            SN.RA(Isn) = convertdms([RA(1), RA(2) RA(3)],'H','r');
         else
            SN.RA(Isn) = NaN;
         end
         if (length(Dec)==1),
            SN.Dec(Isn) = Sign.*abs(Dec(1))./RAD;
         elseif (length(Dec)==2),
            SN.Dec(Isn) = convertdms([Sign abs(Dec(1)), Dec(2)],'DM','R');
         elseif (length(RA)==3),
            SN.Dec(Isn) = convertdms([Sign abs(Dec(1)), Dec(2) Dec(3)],'D','R');
         else
            SN.Dec(Isn) = NaN;
         end
      end      
      SN.DiscRef{Isn} = Cell{I}(72:87);
      SN.PosRef{Isn}  = Cell{I}(114:130);
      if (length(Cell{I})<137),
	 SN.Type{Isn}    = NaN;
         SN.Discoverer   = NaN;
      else
         SN.Type{Isn}    = Cell{I}(131:136);
         SN.Discoverer   = Cell{I}(145:end);
      end
   end   
end

SN.update = date;

