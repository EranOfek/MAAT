function [List,TableCell]=prep_function_list(List)
% Prepare list of all functions in the Astro Toolbox
% Package: Util.code
% Description: Prepare list of all functions in the Astro Toolbox under
%              current directory.
% Input  : - Input argument for recursive calls.
% Output : - List of all functions.
%          - Cell array of selected data for each function.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Feb 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [List,Table]=Util.code.prep_function_list
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==0)
    List = {};
end

AddDir = {};

DirP = dir('+*');
DirC = dir('@*');
DirM = dir('*.m');
%DirD = Util.files.dir_cell({'ImSpec','ImPhot','ImAstrom','plotting','Swift'});
NdirP = numel(DirP);
NdirC = numel(DirC);
NdirM = numel(DirM);
%NdirD = numel(DirD);

% go over all packages:
for IdirP=1:1:NdirP
    cd(DirP(IdirP).name)
    List = Util.code.prep_function_list(List);
    cd ..
end
% go over all classes
for IdirC=1:1:NdirC
    cd(DirC(IdirC).name)
    List = Util.code.prep_function_list(List);
    cd ..
end
% % go over all directories
% for IdirD=1:1:NdirD
%     cd(DirD(IdirD).name)
%     List = Util.code.prep_function_list(List);
%     cd ..
% end

% go over all functions
Nl = numel(List);
for IdirM=1:1:NdirM
    %Nl = Nl+1;
    %List=init_List(List,Nl);
   
    
    %List(Nl).Path = DirM(IdirM).folder;
    
    FunName = DirM(IdirM).name;
    FunName
    
    if (~strcmp(FunName,'startup.m'))
    
    
        FileCell = Util.files.file2str(FunName,'cell');
        %List(Nl).Nlines = numel(FileCell);

        if any(~Util.cell.isempty_cell(strfind(FileCell,'classdef '))) && ~strcmp(FunName,'prep_function_list.m')
            % class definition and methods


            Imethod  = find(~Util.cell.isempty_cell(strfind(FileCell,'methods')));
            ImethodS = find(~Util.cell.isempty_cell(strfind(FileCell,'methods (Static)')));
            Ifun     = find(~Util.cell.isempty_cell(strfind(FileCell,'function')) & ...
                             Util.cell.isempty_cell(strfind(FileCell,'%')) & ...
                             Util.cell.isempty_cell(strfind(FileCell,'function_handle')));
            Nfun     = numel(Ifun);
            for If=1:1:Nfun
                Nl = Nl + 1;
                List=init_List(List,Nl);
                % add FunName name...
                Lfun = FileCell{Ifun(If)};
                I1=strfind(Lfun,'function');
                I2=strfind(Lfun,'=');
                I3=strfind(Lfun,'(');
                if (isempty(I2))
                    Istart = length('function ')+I1;
                else
                    Istart = I2+1;
                end
                if (isempty(I3))
                    Iend=length(Lfun);
                else
                    Iend=I3-1;
                end
                List(Nl).FunName = Lfun(Istart:Iend);

                LineDesc = FileCell{Ifun(If)+1};
                Icom = strfind(LineDesc,'%');
                List(Nl).ShortDesc = LineDesc(Icom+1:end);
                List(Nl).IsPackage = false;
                List(Nl).IsClass   = true;
                Diff = Ifun(If)-Imethod;
                Diff = min(Diff(Diff>0));
                DiffS = Ifun(If)-ImethodS;
                DiffS = min(DiffS(DiffS>0));
                if (Diff==DiffS)
                    List(Nl).IsStatic = true;
                else
                    List(Nl).IsStatic = false;
                end
                List(Nl).LastUpdate = DirM(IdirM).date(1:11);

                List(Nl).Path = DirM(IdirM).folder;
                if (If==1)
                    List(Nl).Nlines = numel(FileCell);
                end
            end


        else
            % function or methods
            Nl = Nl + 1;
            List=init_List(List,Nl);

            List(Nl).Nlines = numel(FileCell);
            List(Nl).Path = DirM(IdirM).folder;
            List(Nl).FunName = FunName;


            PWD = pwd;
            RE = regexp(PWD,filesep,'split');
            Idir = find(strcmp(RE,'fun'));
            Pc = RE(Idir+1:end);
            List(Nl).Package = Pc;

            fprintf('%4d',Nl);
            fprintf('%40s     ',FunName);
            Npc = numel(Pc);
            for Ipc=1:1:Npc
                fprintf('%20s',Pc{Ipc});
            end
            fprintf('\n');

            List(Nl).ShortDesc = FileCell{2}(2:end);

            %List(Nl).ShortDesc
            FunName
            List(Nl).ShortDesc
            if numel(List(Nl).ShortDesc)<2
                List(Nl).Formatted = false;
            else
                if ~isempty(strfind(List(Nl).ShortDesc(1:2),'-'))
                    List(Nl).Formatted = false;
                else
                    List(Nl).Formatted = true;
                end
            end

            if (strcmp(Pc(1),'+'))
                List(Nl).IsPackage = true;
            else
                List(Nl).IsPackage = false;
            end
            if (strcmp(Pc(1),'@'))
                List(Nl).IsClass = true;
            else
                List(Nl).IsClass = false;
            end


            List(Nl).LastUpdate = DirM(IdirM).date(1:11);

        end
    end
end

if (nargout>1)
    % prep web page with list of functions
    Nlist = numel(List);
    TableCell = cell(Nlist,5);
    for Ilist=1:1:Nlist
        FunPath = regexprep(List(Ilist).Path,'/home/eran/','./');
        FunPath = sprintf('%s%s%s',FunPath,'/',List(Ilist).FunName);
        TableCell{Ilist,1} = sprintf('<a href="%s">%s</a>',FunPath,List(Ilist).FunName);
        TableCell{Ilist,2} = List(Ilist).Path(23:end);
        TableCell{Ilist,3} = List(Ilist).IsStatic;
        TableCell{Ilist,4} = List(Ilist).LastUpdate;
        TableCell{Ilist,5} = List(Ilist).ShortDesc;
    end
    % generate HTML table:
    %www.html_table('Table.html','TableCell',TableHead);
end


end

function List=init_List(List,Nl)
    List(Nl).FunName    = NaN;
    List(Nl).Path       = NaN;
    List(Nl).Package    = NaN;
    List(Nl).ShortDesc  = NaN;
    List(Nl).Formatted  = false;
    List(Nl).IsPackage  = false;
    List(Nl).IsClass    = false;
    List(Nl).LastUpdate = NaN;
    List(Nl).IsStatic   = false;
    List(Nl).Nlines     = 0;
end
