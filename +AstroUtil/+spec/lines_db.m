function Ans=lines_db(Line,varargin)
% Search spectral line by name or wavelength and add ionization potential
% Package: AstroUtil.spec
% Description: Search spectral line by name or wavelength and add
%              ionization potential information.
%              The DB list of lines contains 46663 spectral lines
%              (Reader et al. 1980; 1981) for 99 atomic species.
%              Neutral through quadruply ionized atoms are tabulated. 
%              The Main list of lines is a smaller list of prominant lines.
% Input  : - Line wavelength, or name (e.g., 'H I').
%            If empty then upload entire LinesDB.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'SearchWidth' - Wavelength search semi width.
%                            Default is 10 Ang.
%            'SortBy'      - Sort output by one of the following columns:
%                            'Z'|'Intensity'|'Wave'|'Name'|'IonPotential'.
%                            Default is 'Z'.
%            'AddIonPot'   - Add ionization potential to table.
%                            Default is true.
%            'q'           - Apply a query to catalog.
%                            E.g., 'IonPotential>10 & Z<10'.
%                            Default is empty.
%            'DB'          - Lines database to use. Options are:
%                            'DB' - Main database of over 46,000 lines.
%                            'Main' - List of main prominant lines.
%                            Default is 'DB'.
% Output : - An AstCat object containing the requested lines.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Apr 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Ans=AstroUtil.spec.lines_db('H I');
%          Ans=AstroUtil.spec.lines_db(6564,'SearchWidth',2);
%          Ans=AstroUtil.spec.lines_db;
%          Ans=AstroUtil.spec.lines_db(6564,'q','IonPotential>10 & Z<10','SortBy','IonPotential');
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==0)
    Line = [];
end

DefV.SearchWidth          = 10;
DefV.SortBy               = 'Z';  % 'Z'|'Intensity'|'Wave'|'Name'|'IonPotential'
DefV.AddIonPot            = true;
DefV.q                    = [];
DefV.DB                   = 'DB';  % 'DB'|'Main'
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% load spectral lines DB
switch lower(InPar.DB)
    case 'db'
        LinesDB = cats.spec.lines.LinesDB;
    case 'main'
        LinesDB = cats.spec.lines.MainLines;
    otherwise
        error('Unknown DB option');
end


if (isempty(Line))
    Ans = LinesDB;
else
    if (ischar(Line))
        % search by line name
        Ifl = find(strcmp(Line,LinesDB.Cat.Name));
    else
        % search by line wavelength
        Ifl = find(abs(LinesDB.Cat.Wave-Line)<InPar.SearchWidth);
    end

    Ans = LinesDB;
    Ans.Cat = LinesDB.Cat(Ifl,:); 



    % add ionization potential
    if (InPar.AddIonPot)
        [P]=AstroUtil.spec.ionization_potential(Ans.Cat.Z,Ans.Cat.Ion);
        Ans.Cat = [Ans.Cat,table(P)];
        Ans.ColCell  = [Ans.ColCell, 'IonPotential'];
        Ans.ColUnits = [Ans.ColUnits, 'eV'];
        Ans = colcell2col(Ans);
        Ans.Cat.Properties.VariableNames = Ans.ColCell;
        Ans.Cat.Properties.VariableUnits = Ans.ColUnits;
    end

    
end

% apply query
if (~isempty(InPar.q))
    Ans = query(Ans,InPar.q);
end

% sort
if (~isempty(InPar.SortBy))
    Ans.Cat = sortrows(Ans.Cat,InPar.SortBy);
end