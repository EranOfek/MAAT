function Out=get_lines1(In,MaxDist);
%------------------------------------------------------------------------------
% get_lines1 function                                                AstroSpec
% Description: Search for a spectral line, by name or by eavelength,
%              among a list of selected spectral lines.
% Input  : - Line name (e.g., 'H \alpha') or
%            column vector of lines restframe wavelength [Ang]
%            or two column matrix of ranges of lines restframe
%            wavelength [Ang].
%            If String or cell array of strings then search for matching
%            line names in the Spectral Line database.
%            If a single column vector then  return the nearest line.
%            If two column matrix then return all the lines within the
%            region defined by each row.
%          - Max dist [Ang].
%            In case that the single column input is used then
%            return only lines found within this distance.
%            If negative number then returns all the lines and
%            not only the nearest line.
%            Default is 10 Ang.
% Output : - Cell array of lines wavelength and name, or
%            column vector of the wavelength of all the lines found.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                      July 2005
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
% See Also: get_lines.m
% Example: get_lines1('H \alpha');
%          get_lines1(6563);
% Needed: LinesDB_Selected.mat   (database of selected spectral line)
% Reliable: 1
%-----------------------------------------------------------------------
MaxDist_Def = 10;

if (nargin==1),
   MaxDist  = MaxDist_Def;
elseif (nargin==2),
   % do nothing
else
   error('Illegal number of input arguments');
end


if (MaxDist<0),
   MaxDist   = abs(MaxDist);
   FindLines = 'all';
else
   FindLines = 'nearest';
end

if (isstr(In)),
   In1{1} = In;
   clear In;
   In     = In1;
end

load LinesDB_Selected.mat;
Ndb     = length(SpLines.Line);

Out     = zeros(0,1);
if (iscell(In)),
   % search by line name
   Nin = length(In);

   for Iin=1:1:Nin,
      Out(Iin) = NaN;
      for Idb=1:1:Ndb,
         if (strcmpi(In{Iin},SpLines.Name{Idb})),
            % line found - return wavelength
            Out(Iin) = SpLines.Line{Idb};
         end
      end
   end

else
   % search by wavelength
   [Nin,Nc] = size(In);
   
   if (Nc==2),
      %--- Convert to single column case ---
      FindLines = 'all';
      MaxDist   = (max(In,[],2) - min(In,[],2)).*0.5;
      In        = mean(In,2);
   end
   if (length(MaxDist)==1),
      MaxDist = MaxDist.*size(In);
   end

   Out  = cell(0,2);
   Iout = 0;
   for Idb=1:1:Ndb,
      Diff = SpLines.Line{Idb} - In;

      switch FindLines
       case 'all'
          Ifound = find(abs(Diff)<=MaxDist);
          if (isempty(Ifound)),
             % not found
          else
             % line found
             Iout        = Iout + 1;
             Out{Iout,1} = SpLines.Line{Idb};
             Out{Iout,2} = SpLines.Name{Idb};
          end
       case 'nearest'
          [Min,MinI] = min(abs(Diff));
          if (Min>MaxDist),
             % not found
          else
             % line found
             Iout        = Iout + 1;
             Out{Iout,1} = SpLines.Line{Idb};
             Out{Iout,2} = SpLines.Name{Idb};
          end
       otherwise
          error('Unknown FindLines Option');
      end
   end
end
