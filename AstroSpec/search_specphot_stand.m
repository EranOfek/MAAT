function Res=search_specphot_stand(RA,Dec,Size,OutType)
%--------------------------------------------------------------------------
% search_specphot_stand function                                 AstroSpec
% Description: Search a spectroscopic standard star, by coordinates or by
%              name in the SpecPhot_Stand.mat structure.
%              The search by name option is case insensitive and blanks
%              insensitive. Moreover, the name search can be either exact
%              or by substring.
% Input  : - Standard star name.
%            Alternatively this can be a J2000.0 R.A. In this case the
%            second input argument is J2000.0 Dec (see convertdms.m
%            for options).
%          - J2000.0 Dec. see convertdms.m for details).
%            Alternatively, this can be one of the following strings:
%            'exact' - to perform an exact name search (default).
%            'substr'- to search for substring within the names.
%          - Search radius [arcsec]. Default is 60.
%          - Spectrum output type {'mat'|astspec'}.
%            Default is 'astspec'.
% Output : - Structure array of all stnadtad stars found.
% Tested : Matlab 2011b
%     By : Eran O. Ofek                    Apr 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: add_specphot_stand.m
% Example: Res=search_specphot_stand('Hz44');
%          Res=search_specphot_stand('Hz','substr');
%          Res=search_specphot_stand('21:51:11.021',[1 28 51 50.36]);
% Reliable: 2
%--------------------------------------------------------------------------
import Util.string.*

RAD  = 180./pi;

Def.Size = 60;
Def.NameType = 'exact';
Def.OutType  = 'astspec';
if (nargin==1),
   % Search by name
   SearchMethod = 'name';
   NameType     = Def.NameType;
   OutType      = Def.OutType;
elseif (nargin>1),
   if (ischar(Dec)),
      switch lower(Dec)
       case {'exact','substr'}
          SearchMethod = 'name';
          NameType     = Dec;
          OutType      = Def.OutType;
       otherwise
          % search by coordinates
          SearchMethod = 'coo';
          if (nargin==2),
             Size    = Def.Size;
             OutType = Def.OutType;
          elseif (nargin==3),
              OutType = Def.OutType;
          else
              % do nothing
          end
      end
   else
      % Search by coordinates
      SearchMethod = 'coo';
      if (nargin==2),
         Size    = Def.Size;
         OutType = Def.OutType;
      elseif (nargin==3),
          OutType = Def.OutType;
      else
          % do nothing
      end
   end
else
   error('Illegal number of input arguments');
end

load SpecPhot_Stand.mat


IndF = 0;
switch SearchMethod
 case 'name'
    % search by name
    Name = RA;
    N = length(SpecPhot_Stand);

    Res = struct('Name',{},'RA',{},'Dec',{},'MagV',{},'SpecType',{},'Spec',{},'Comment',{},'GridWave',{});
    for In=1:1:N,
       switch lower(NameType)
        case 'exact'
           C = strcmp(lower(spacedel(SpecPhot_Stand(In).Name)),...
   		      lower(spacedel(Name)));
        case 'substr'

	       C = ~Util.cell.isempty_cell(strfind(lower(spacedel(SpecPhot_Stand(In).Name)),...
			      lower(spacedel(Name)) ));


        otherwise
	        error('Unknown NameType option');
       end
       if (~isempty(find(C)==1)),
          % found
	      IndF = IndF + 1;
          Res(IndF) = SpecPhot_Stand(In);  
          
       end
    end
    
    
 case 'coo'
    % search by coordinates
    D = sphere_dist([SpecPhot_Stand.RA].',[SpecPhot_Stand.Dec].',RA,Dec);
    I = find(D<(Size./(RAD.*3600)));

    if (~isempty(I))
       Res = SpecPhot_Stand(I);
    else
       % not found
       Res = [];
    end

 otherwise
    error('Illegal number of input arguments');
end

switch lower(OutType)
    case 'astspec'
        Res.Spec         = AstSpec.mat2astspec(Res.Spec);
        Res.Spec.z       = 0;
        Res.Spec.ObjName = Res.Name;
    case 'mat'
        % do nothing
    otherwise
        error('Unknown OutType option');
end
