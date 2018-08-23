%--------------------------------------------------------------------------
% InArg class                                                        class
% Description: A static class to deal with key,val input arguments, set
%              their default and edit their values.
%              Type "InArg." followed by <tab> to see the full list of
%              functions.
% Input  : null
% Output : null
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jul 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef InArg 

    %----------------------
    %--- Static methods ---
    %----------------------
    methods (Static)
        
        % Return the file path in which the default arguments are stored.
        function DefArgPath = default_arg_path
            % file path in which the default arguments are stored
            % Description: Return the file path in which the default
            %              arguments are stored.
            %              By default this function assumes that the
            %              default arguments are stored in:
            %              '
            % Input  : null
            % Output : - The path of the directory containing the default
            %            argument files.
            
            
            % Get the user name:
            UserHome   = Util.OS.get_userhome;
            % UserName = char(java.lang.System.getProperty('user.name'))
            % char(java.net.InetAddress.getLocalHost.getHostName);
            DefFileLocalPath = 'matlab/.FunPars';
            DefArgPath    = sprintf('%s%s%s%s',UserHome,filesep,DefFileLocalPath,filesep);
            
        end
        
        % Name that contains the file name of the default arguments
        function DefArgFile = default_arg_file(CallerFun)
            % file name that contains the file name of the default arguments
            % Description: Return the file name that contains (without full
            %              path) the file name of the default arguments for
            %              the function that called this function.
            % Input  : - A string of the caller function name, or a the
            %            relative depth of the caller function in the
            %            dbstack.
            %            Default i2 s.
            % Output : - The file name of the default arguments for
            %            the function that called this function.
            DefFileSuffix   = '.par.mat';
            
            Def.CallerFun = 2;
            if (nargin==0)
                CallerFun = Def.CallerFun;
            end
            if (isempty(CallerFun))
                CallerFun = Def.CallerFun;
            end
            
            if (ischar(CallerFun))
                % do nothing
            else
                % CallerFun is the depth of call
                [ST] = dbstack;                % get the name of caller functions
                if (length(ST)==1)
                    CallerFun = 'session';
                else
                    CallerFun = ST(CallerFun).name;        % Name of caller function 
                end
            end
            DefArgFile = sprintf('%s%s',CallerFun,DefFileSuffix);
            
            
        end
            
        % Save the default arguments file
        function save_arg_file(CallerFun,InPar)
            % Save the default arguments file
            % Description: Save the default arguments file
            % Input  : - The caller function name.
            %          - The InPar structure containing the keyword and
            %            values to save.
            %            Default is to parse it from the function.
            % Output : null.
           
            if (nargin==1)
                % parse default arguments from function
                InPar = InArg.default_arg(CallerFun);
            end
            
            Path     = InArg.default_arg_path;
            FileName = InArg.default_arg_file(CallerFun);
            FullFileName = sprintf('%s%s',Path,FileName);
            save(FullFileName,'InPar');
            
        end
        
        % Load default arguments file
        function InPar=load_arg_file(CallerFun)
            % Load default arguments file into the InPar structure
            % Description: Load default arguments file into the InPar
            %              structure.
            % Input  : - Caller function name.
            % Output : - A structure containing the keyword and values
            %            Return empty if file not found.
            Path     = InArg.default_arg_path;
            FileName = InArg.default_arg_file(CallerFun);
            FullFileName = sprintf('%s%s',Path,FileName);
            if (java.io.File(FullFileName).exists)
                InPar = Util.IO.load2(FullFileName);
            else
                InPar = [];
            end
                
        end
        
        % Delete a default arguments file
        function delete_arg_file(CallerFun)
            % Delete a default arguments file
            % Description: Delete a default arguments file.
            % Input  : - Caller function name.
            % Output : null.
            Path     = InArg.default_arg_path;
            FileName = InArg.default_arg_file(CallerFun);
            FullFileName = sprintf('%s%s',Path,FileName);
            delete(FullFileName);
        end
        
        % Populate a structure array from key,val pairs
        function InPar=populate_keyval(DefV,VarArgIn,CallerFun,ArgExistInDefV,CheckDefFile)
            %A structure array containing the default arguments and their values
            % Input  : - A structure array containing the default
            %            arguments and their values. The field name
            %            corresponds to the argument name and the field
            %            value to the argument value.
            %          - A cell array of the key,val pairs of input
            %            arguments (i.e., varargin).
            %          - The caller function name. If empty then use
            %            default_arg_file.m. Default is empty.
            %          - A flag indicating if to check if the arguments in
            %            keyword name appears in the DefV field names.
            %            Default is true.
            %          - A flag indicating if to check if the default
            %            arguments file exist and if so to use it.
            %            Default is true
            
            
            Def.CallerFun      = [];
            Def.ArgExistInDefV = true;
            Def.CheckDefFile   = true;
            if (nargin==2)
                CallerFun      = Def.CallerFun;
                ArgExistInDefV = Def.ArgExistInDefV;
                CheckDefFile   = Def.CheckDefFile;
            elseif (nargin==3)
                ArgExistInDefV = Def.ArgExistInDefV;
                CheckDefFile   = Def.CheckDefFile;
            elseif (nargin==4)
                CheckDefFile   = Def.CheckDefFile;
            elseif (nargin==5)
                % do nothing
            else
                error('Illegal number of input arguments: populate_keyval(DefV,VarArgIn,[CallerFun,ArgExistInDefV,CheckDefFile])');
            end
            
            
                
            if (CheckDefFile)
                DefArgFile = InArg.default_arg_file(CallerFun);
                DefArgPath = InArg.default_arg_path;
                DefArgFilePath = sprintf('%s%s',DefArgPath,DefArgFile);
                
                % Note that using the java method is much faster than
                % the exist function
                if (java.io.File(DefArgFilePath).exists)
                    % User default parameter file was found - use it
                    % instead of the DefV
                    DefV = Util.IO.load2(DefArgFilePath);
                end
%                 if (exist(DefParFile,'file')>0),
%                     % User default parameter file was found - use it
%                     % instead of the DefV
%                     DefV = load(DefParFile);
%                 end
            end
            
            Narg  = numel(VarArgIn);
            if (Narg==0)
                % Return the Default arguments as is
                InPar = DefV;
            else
                % Modify the default argument according to the user request
                
                if (Narg.*0.5~=round(Narg.*0.5))
                    error('Must supply an even number of key,val arguments');
                end

                InPar = DefV;
                
                FieldNames = fieldnames(DefV);
                Nf         = numel(FieldNames);
                for Iarg=1:2:Narg-1
                    FoundFlag = strcmpi(VarArgIn{Iarg},FieldNames);
                    if all(~FoundFlag)
                        % User supplied argument is not in DefV structure
                        if (ArgExistInDefV)
                            error('User supplied argument (%s) is not in default arguments structure',VarArgIn{Iarg});
                        else
                            % Ignore arguments which are not in DefV
                            
                        end
                    else
                        % User supplied argument found in DefV structure
                        % set value to user supplied input
                        CurrField = FieldNames{FoundFlag};
                        %InPar.(VarArgIn{Iarg}) = VarArgIn{Iarg+1};
                        InPar.(CurrField) = VarArgIn{Iarg+1};
                    end
                end
            end
            
        end
        
        % Parse default arguments of a matlab function
        function DefV=default_arg(FunName)
            % Parse default arguments of a matlab function
            % Description: Parse default arguments of a matlab function
            %              that uses the DefV, InArg.populate_keyval syntax.
            % Input  : - Function name.
            % Output : - A structure with the function default keyword and
            %            their values.
            
            FunName_m = sprintf('%s.m',FunName);
            FunName   = which(FunName);

            Line  = Util.files.file2str(FunName_m,'cellremove3dots');
            Idefv = strfind(Line,'DefV.');
            Iline = find(~Util.cell.isempty_cell(Idefv));
            N     = numel(Iline);
            for I=1:1:N
               evalc(Line{Iline(I)});
            end

%            FID = fopen(FunName,'r');
%            Cont = true;
%            while (~feof(FID) || Cont)
%                Line = fgetl(FID);
%                if (isnumeric(Line))
%                    Cont = false;
%                end
%                if (numel(Line)>4)
%                    if (~isempty(strfind(Line,'InArg.populate_keyval')))
%                        Cont = false;
%                    end
%                    if (strcmp(Line(1:5),'DefV.'))
%                        % Need to take care of cont line (...)
%                        if (strcmp(Util.string.spacedel(Line(end-2:end)),'...'))
%                            error('default_arg does not work on cont lines');
%                        end
%                            
%                        evalc(Line);
%                    end
%                end
%            end
%            fclose(FID);
            
        
        end
        
        % edit arguments
        function epar(FunName)
            % Edit the default arguments file for a function
            % Description: Edit the default arguments file for a function
            % Input  : - Function name.
            % Output : - A structure array with the editted parameters.
            
            % Get user default arguments from file
            DefVar = InArg.load_arg_file(FunName);
            if (isempty(DefVar))
                % arguments file doesn't exist create the file
                InArg.save_arg_file(FunName);
                % Get user default arguments from file
                DefVar = InArg.load_arg_file(FunName);
            end
            
            
            assignin('base','DefVar',DefVar);

            %openvar('DefVar');
            evalin('base','openvar(''DefVar'');');
            drawnow;
            pause(0.5); % wait for client to become visible

            % Get handle of variables client in the Variable Editor
	        %http://blogs.mathworks.com/community/2008/04/21/variable-editor/
            jDesktop = com.mathworks.mde.desk.MLDesktop.getInstance;
            jClient = jDesktop.getClient('DefVar');
            hjClient = handle(jClient,'CallbackProperties');

            % Instrument the client to fire callback when it's closed
            set(hjClient,'ComponentRemovedCallback',{@InArg.my_callback_epar,FunName,DefVar});

            % reload changes
            %DefVar = load2(FunParsFullName);
            %InPar = DefVar;
            
            
        end
        
        % Aux function for InArg.epar
        function my_callback_epar(varEditorObj,eventData,FunName,DefVar)
            % Description: Aux function for InArg.epar
            
            % user closed openvar editor - get the DefVar from the session
            DefVar = evalin('base','DefVar');
            % Save it in the user default arguments file
            InArg.save_arg_file(FunName,DefVar)
            
        end
        
        % unlearn arguments
        function InPar=unlearn(FunName)
            % Set the user default arguments file of a function to its default state
            % Description: Set the user default arguments file of a
            %              function to its default state (i.e., the default
            %              parameters in the function).
            % Input  : - Function name.
            % Output : - A structure array with the arguments.
            
            % Delete the user arguments file
            InArg.delete_arg_file(FunName);

            % Get the function default parameters
            InPar = InArg.default_arg(FunName);
            
            % Save the default arguments
            InArg.save_arg_file(FunName,InPar);
            
        end
    end
end
       
