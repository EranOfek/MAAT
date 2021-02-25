
% A logFile class to handle log files
% Package: @logFile
% Description: 
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Fen 2020
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
% Example:
%   % create a logFile object, and write messages
%     L = logFile
%     L.writeLog('executed command x');
%     L.writeLog('executed command y');
%   % to see the default log file name and directory in which the info is
%   % stored:
%     L.FileName, L.Dir
%   % the file name is still open - so in order to see its content
%   % either delete the object
%     delete(L)
%   % or close the file
%     L.closeFile
%


classdef logFile < handle
    properties (SetAccess = public)
        FileNameTemplate     = 'logFile_%s_%s.log';   % default is logFile_<logOwner>_<YYYYMMDD>.log | logFile_<logOwner>.log
        FileName     = '';
        Dir          = pwd;                   % default is pwd
        logOwner     = 'unknown';             % e.g., the process/instrument that uses the logFile
        Verbose      = false                  % {true,[false]}
        FID          = [];
        Counter      = 0;
        
    end
  
    % constructor/destructor
    methods

        %-------------------------
        %--- Class constructor ---
        %-------------------------
        
        function H=logFile(varargin)
            % HEAD constructor
            % Package: @HEAD
            
            
            %FullPath = sprintf('%s%s%s',H.Dir,filesep,H.FileName);
            
        end
        
        function delete(H)
            % delete logFile class and close log file
            % Package: @logFile
            
            if ~isempty(H.FID)
                fclose(H.FID);
                H.FID = [];
            end
            
            delete(H);
        end
 
        function closeFile(H)
            % close FileName
            % Package: @logFile
            
            if ~isempty(H.FID)
                fclose(H.FID);
                H.FID = [];
            end
            
            
        end
        
        
    end
    

    % private utilities
    methods 
        function update_FileName(H)
            % update FileName according to date
            % package: @logFile
            
            switch numel(strfind(H.FileNameTemplate,'%s'))
                case 0
                    % file name doesn't include date template
                    FileName = H.FileNameTemplate;
                case 1
                    % FileName includes logOwner template
                    FileName = sprintf(H.FileNameTemplate,H.logOwner);
                case 2
                    % FileName includes logOwner and date template
                    DateStr = datestr(now,'yyyymmdd');
                    FileName = sprintf(H.FileNameTemplate,H.logOwner,DateStr);
                otherwise
                    error('Unknown FileName template');
            end

                
            if strcmp(FileName,H.FileName)
                % new file name = old file name
                % do not do anything
                
                if isempty(H.FID)
                    H.FID = fopen(H.FileName,'a+');
                end
                
            else
                % FileName was updated
                % close old file
                if ~isempty(H.FID)
                    fclose(H.FID)
                    H.FID = [];
                end
                % open new file ID
                H.FID = fopen(FileName,'a+');
                H.FileName = FileName;
            end
        end
    end
    
    % getters/setters
    methods
        function H=set.logOwner(H,logOwnerString)
            % setter for logOwner - update FileName
            % Package: @logFile
            % Input  : - logOwnerString
        
            if ~ischar(logOwnerString)
                error('logOwner must be a string');
            end
            H.logOwner = logOwnerString;
            
            % update FileName
            H.update_FileName;
        
        end
        
         function H=set.FileName(H,FileNameString)
            % setter for FileName - update FileName
            % Package: @logFile
            % Input  : - FileNameString | NaN - NaN will restore default
            
        
            
            if ~ischar(FileNameString)
                if isnan(FileNameString)
                    H.FileName = 'logFile_%s_%s.log';
                else
                    error('FileName must be a string or NaN');
                end
            else
                H.FileName = FileNameString;
            end
            
            
            % update FileName
            H.update_FileName;
        
         end
        
    end
    
    % write logFile message
    methods
        function writeLog(H,Message)
            % writeing a message into the logFile
            % Package: @logFile
            % Input  : - 
            
            H.update_FileName;  % update file name
            
%             if isempty(H.FID)
%                 % fileID doesn't exist - open a new file
%                 H.FID = fopen(sprintf('%s%s%s',H.Dir,filesep,H.FileName),'a+');
%             end
            
            % write
            DateStr = datestr(now,'yyyy-mm-dd HH:MM:SS.FFF');
            [SI,I] = dbstack;
            fprintf(H.FID,'%s %s %s\n',DateStr,H.logOwner,Message);
            H.Counter = H.Counter + 1;
            
            if (H.Verbose)
                % print message to screen
                fprintf('%s %s %s\n',DateStr,H.logOwner,Message);
            end
            
        end
    end % end methods
end

            
