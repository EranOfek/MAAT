%--------------------------------------------------------------------------
% ds9 class                                                          class
% Description: A static class for intearction with the ds9 display.
%              This include functions to load images, change their
%              properties, create and plot region files, printing, image
%              examination, interaction with SIM content and more.
%              Type "ds9." followed by <tab> to see the full list of
%              functions.
%              Full manual is available in manual_ds9.pdf
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Jul 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef ds9
    % A static class
    
    % Constructor method (display)
    methods
        function Obj=ds9(Image,varargin)
            % Display images in ds9
            % Package: @ds9
            % Description: Display images in ds9 (constructor).
            %              See also ds9.disp.
            % Input  : - Images to display. One of the following:
            %            1. String containing a single FITS image name.
            %            2. String with wild cards and ranges for multiple
            %               FITS images.
            %            3. Cell array of FITS image names.
            %            4. A matrix.
            %            5. A cube (display a cube).
            %            6. Cell array of matrices.
            %            7. A SIM object.
            %          * A vector of frames (see 'Frame') in which to load
            %            the images into followed by ...,key,val,... pairs,
            %            or just the ...,key,val,... pairs.
            %            'Frame' - Frame indices or frame parameters
            %                      (e.g., 'next'). Default is (1:N), where
            %                      N is the number of images to display.
            %            'Scale' - See ds9.scale. Default is 'mode zscale'.
            %            'CMap'  - See ds9.cmap. Default is [].
            %            'Colorbar' - See ds9.colorbar. Default is 'no'.
            %            'Orient'- See ds9.orient. Default is [].
            %            'Rotate'- See ds9.rotate. Default is [].
            %            'Zoom'  - See ds9.zoom. Default is 1.
            % Output : null
            % Example: ds9(rand(100,100),1,'Zoom',2)
            %          ds9(rand(100,100),'Zoom',2);
            %          ds9(Sim(1:3));
            %          ds9('MyFitsImage.fits');
            %          ds9('Images*.fits');
            % Reliable: 2
            
            % ds9_disp(Image,varargin{:});
            ds9.disp(Image,varargin{:});
            clear Obj; % Output Object Not Assigned, therefore need to clear it
        end
        
    end
    
    
    
    % Static methods for ds9
    % Basic methods (system, command construction, open, exit,...)
    methods (Static)
        
        % execute xpa command
        function [Answer,Status]=system(String,varargin)
            % Construct and execute an xpa command
            % Package: @ds9
            % Input  : - String to execute. The string may include any
            %            special characters (e.g., %s, %d) that will be
            %            populated by the additional parameters.
            %          * Arbitrary number of additional parameters that
            %            will populate the %s %d like format.
            % Output : - Answer
            % Example: ds9.system('xpaset -p ds9 frame frameno %d',FrameNumber);
            % Reliable: 2
            String = sprintf(String,varargin{:});
            if ismac
               String = strcat('set DYLD_LIBRARY_PATH "";', String);
               [Status,Answer]=system(String);
            else
               [Status,Answer]=system(String);
            end
            if (Status~=0)
                if (strfind(Answer,'not found'))
                    % It is possible that xpa is not installed
                    fprintf('\n It seems that xpa and or ds9 are not installed\n');
                    fprintf('ds9 instellation: http://ds9.si.edu/site/Download.html\n');
                    fprintf('xpa instellation: http://ds9.si.edu/site/XPA.html\n');
                end
                error('Command: %s failed - Answer: %s',String,Answer);
                
            end
            
        end

        % construct xpa command
        function String=construct_command(varargin)
            % Construct an arguments string for ds9 command
            % Package: @ds9
            % Input  : * Arbitrary number of arguments from which to
            %            construct a string.
            %            Default is ''.
            % Output : The constructed string
            % Example: ds9.construct_command('scale linear');
            %          ds9.construct_command scale linear
            %          ds9.construct_command('log',100); % return 'log 100'
            %          ds9.construct_command('log 100');
            % Reliable: 2
            String = '';
            Narg = numel(varargin);
            for Iarg=1:1:Narg
                if (ischar(varargin{Iarg}))
                    String = sprintf('%s %s',String,varargin{Iarg});
                else
                    if (mod(varargin{Iarg},1)==0)
                        % argument is an integer
                        String = sprintf('%s %d',String,varargin{Iarg});
                    else
                        % argument is a double
                        String = sprintf('%s %f',String,varargin{Iarg});
                    end
                end
            end
            
        end
        
        % xpaset
        function xpaset(Command)
            % Execute an xpaset command
            % Package: @ds9
            % Input  : - The command to follow the 'xpaset -p ds9'.
            % Output : null
            % Reliable: 2
            ds9.system('xpaset -p ds9 %s',Command);
        end
        
        % xpaget
        function Ans=xpaget(Command)
            % Execute an xpaget command
            % Package: @ds9
            % Input  : - The command to follow the 'xpaget ds9'.
            % Output : - Return output.
            % Reliable: 2
            Ans=ds9.system('xpaget ds9 %s',Command);
        end
        
        % xpa help
        function xpahelp
            % Open the XPA command help web page
            % Package: @ds9
            % Input  : null
            % Output : null
            % Example: ds9.xpahelp
            % Reliable: 2
           
            web('http://ds9.si.edu/doc/ref/xpa.html');
            
        end
        
        % check if ds9 is open
        function IsOpen=isopen
            % Check if ds9 is open (linux/mac only)
            % Package: @ds9
            % Input  : null
            % Output : A logical flag indicating if ds9 is open for the
            %          user.
            % Example: ds9.isopen
            % Reliable: 2
            if (ismac)
                Status = system('ps -A | grep ds9 | grep -v grep');
            else
                Status = system('ps -xug | grep ds9 | grep -v grep');
            end
            if (Status==1)
                % not open
                IsOpen = false;
            else
                % open
                IsOpen = true;
            end
            
        end
        
        % open ds9
        function open
            % Open ds9 dispaly window if not exist already and set mode to region
            % Package: @ds9
            % Input  : null
            % Output : null
            % Problems: The command ds9 is not recognized by the bash interpreter, try adding an alias
            %           to the bash profile using ‘vim .bash_profile’ adding the following line to the
            %           profile ‘alias ds9="open -a /Applications/SAOImageDS9.app/Contents/MacOS/ds9"
            % Example: ds9.open
            % Reliable: 2
            Wait = 3;  % s
            IsOpen = ds9.isopen;
            if (IsOpen)
                % do nothing - already open
                fprintf('ds9 is already open\n');
            else
                Status=system('ds9&');
                if (Status~=0)
                    warning('Can not open ds9');
                else
                    pause(Wait);
                end
            
            end
            ds9.mode('region');
            
        end
        
        % mode
        function mode(Mode)
            % Set ds9 mode
            % Package: @ds9
            % Input  : - ds9 mode. 
            %            none|region|crosshair|colorbar|pan|zoom|rotate|catalog|examine.
            %            Default is 'region'.
            % Output : null
            
            if (nargin==0)
                Mode = 'region';
            end
            
            ds9.system('xpaset -p ds9 mode %s',Mode);
        end
        
        % exit ds9
        function exit
            % exit ds9
            % Package: @ds9
            % Input  : null
            % Output : null
            % Reliable: 2
            ds9.system('xpaset -p ds9 exit');
        end
        
    end
    
    % Frame methods
    % (frame, delete_frame, clear_frame, disp, load, load1, disp, write2fits, write2sim)
    methods (Static)
        
        % set/get frame
        function Answer=frame(FrameNumber)
            % Set ds9 frame
            % Package: @ds9
            % Input  : - If not given than only get the current frame 
            %            number. If numeric then set frame number.
            %            Alternatively, if a string than set to frame name.
            %            Possible names: 'first'|'prev'|'next'|'last'
            %            Additional possibilities include:
            %            'hide'|'move first'|'move last'|'move back'|
            %            'move forward'|'match wcs'|'lock wcs'|...
            % Output : - Current frame number.
            % Reliable: 2
            if (nargin==0)
                % get frame number
                Answer = ds9.system('xpaget ds9 frame frameno');
            else
                % set frame number
                if ischar(FrameNumber)
                    % string frame
                    ds9.system('xpaset -p ds9 frame %s',FrameNumber);
                else
                    % numeric frame
                    ds9.system('xpaset -p ds9 frame frameno %d',FrameNumber);
                end
                
                % get current frame number if needed
                if (nargout>0)
                    Answer = ds9.system('xpaget ds9 frame frameno');
                end
            end
            pause(0.2);
        end
        
        % delete frame
        function delete_frame(FrameNumber)
            % Delete frames from ds9
            % Package: @ds9
            % Description: Delete frame number, some frames or all frames
            %              from ds9.
            % Input  : - Vector of frames to delete.
            %            If empty then delete current frame. Default.
            %            If Inf then delete all frames.
            % Output : null
            % Example: ds9.delete_frame(Inf)
            % Reliable: 2
           

             if (nargin==0)
                FrameNumber = [];
             end
            
            if (isempty(FrameNumber))
                % delete current frame
                ds9.system('xpaset -p ds9 frame delete');
            else
                if (isinf(FrameNumber))
                    % get all frame numbers
                    Ans = ds9.system('xpaget ds9 frame all');
                    Frames = regexp(Ans,' ','split');
                    Frames = Frames(1:end-1);
                    FrameNumber = cellfun(@str2double,Frames);
                end
                
                % delete requested frames
                Nf = numel(FrameNumber);
                for If=1:1:Nf
                    ds9.system('xpaset -p ds9 frame frameno %d',FrameNumber(If));
                    ds9.system('xpaset -p ds9 frame delete');
                end
                
            end
        end
              
        % clear frame
        function clear_frame(FrameNumber)
            % Clear frames from ds9
            % Package: @ds9
            % Description: Clear frame number, some frames or all frames
            %              from ds9.
            % Input  : - Vector of frames to clear.
            %            If empty then clear current frame. Default.
            %            If Inf then clear all frames.
            % Output : null
            % Example: ds9.clear_frame(Inf)
            % Reliable: 2
           
            if (nargin==0)
                FrameNumber = [];
            end
                
            
            if (isempty(FrameNumber))
                % clear current frame
                ds9.system('xpaset -p ds9 frame clear');
            else
                if (isinf(FrameNumber))
                    % get all frame numbers
                    Ans = ds9.system('xpaget ds9 frame all');
                    Frames = regexp(Ans,' ','split');
                    Frames = Frames(1:end-1);
                    FrameNumber = cellfun(@str2double,Frames);
                end
                
                % delete requested frames
                Nf = numel(FrameNumber);
                for If=1:1:Nf
                    ds9.system('xpaset -p ds9 frame frameno %d',FrameNumber(If));
                    ds9.system('xpaset -p ds9 frame clear');
                end
                
            end
        end
        
        % Load a FITS image into ds9 frame
        function [Frame,CurFile]=load(FitsFile,FrameNumber)
            % Load a FITS image into ds9 frame
            % Package: @ds9
            % Description: Load a FITS image into ds9 frame
            % Input  : - FITS file name. This is either a fits file name
            %            in the current directory or a full path to the
            %            file name.
            %          - Frame number.
            %            Default is 'next'.
            % Output : - Current frame number
            %          - Current fits file name
            % Reliable: 2
            if (nargin==1)
                FrameNumber = 'next';
            end
            
            Frame = ds9.frame(FrameNumber);
            
            % check if file name contains directory seperator
            if (~isempty(strfind(FitsFile,filesep)))
                % file name contains directory seperator
                % Assume that FitsFile is a full path
                % do nothing
            else
                % Assume FitsFile is not a full path
                FitsFile = sprintf('%s%s%s',pwd,filesep,FitsFile);
            end
            
            % load FITS file 
            ds9.system('xpaset -p ds9 file %s',FitsFile);
            
            if (nargout>1)
                CurFile = ds9.system('xpaget ds9 file');
            end

        end
        
        % Load a FITS image into ds9 frame number 1
        function [Frame,CurFile]=load1(FitsFile)
            % Load a FITS image into ds9 frame number 1
            % Package: @ds9
            % Description: Load a FITS image into ds9 frame number 1
            % Input  : - Fits file name
            % Output : - Current frame number
            %          - Current fits file name
            % Reliable: 2
            [Frame,CurFile] = ds9.load(FitsFile,1);
            
        end
        
        function Frame=url(URL,FrameNumber)
            % Load FITS file from a URL
            % Package: @ds9
            % Description: Load FITS file from a URL.
            % Input  : - URL string.
            %          - Frame number.
            %            Default is 'next'.
            % Output : - Frame number.
            % Example: ds9.url(URL);
            % Reliable: 2
            
            if (nargin==1)
                FrameNumber = 'next';
            end
            
            Frame = ds9.frame(FrameNumber);
            ds9.system('xpaset -p ds9 url %s',URL);
        end
        
        % display images in all formats
        function disp(Images,varargin)
            % Display images in ds9 (use ds9 for short cut)
            % Package: @ds9
            % Description: Display images in ds9 
            %              For shortcut use the ds9 constructor
            %              (see examples).
            % Input  : - Images to display. One of the following:
            %            1. String containing a single FITS image name.
            %            2. String with wild cards and ranges for multiple
            %               FITS images.
            %            3. Cell array of FITS image names.
            %            4. A matrix.
            %            5. A cube (display a cube).
            %            6. Cell array of matrices.
            %            7. A SIM object.
            %          * A vector of frames (see 'Frame') in which to load
            %            the images into followed by ...,key,val,... pairs,
            %            or just the ...,key,val,... pairs.
            %            'Frame' - Frame indices or frame parameters
            %                      (e.g., 'next'). Default is (1:N), where
            %                      N is the number of images to display.
            %            'Scale' - See ds9.scale. Default is 'mode zscale'.
            %            'CMap'  - See ds9.cmap. Default is [].
            %            'Colorbar' - See ds9.colorbar. Default is 'no'.
            %            'Orient'- See ds9.orient. Default is [].
            %            'Rotate'- See ds9.rotate. Default is [].
            %            'Zoom'  - See ds9.zoom. Default is 1.
            % Output : null
            % Example: ds9.disp(rand(100,100),1,'Zoom',2)
            %          ds9.disp(rand(100,100),'Zoom',2);
            %          ds9.disp(Sim(1:3));
            %          ds9.disp('MyFitsImage.fits');
            %          ds9.disp('Images*.fits');
            % Reliable: 2
            
            Narg = numel(varargin);
            if (Narg.*0.5)~=floor(Narg.*0.5)
                % odd number of arguments
                % assume 1st argument is frame
                Frame = varargin{1};
                varargin = varargin(2:end);
            else
                Frame = [];
            end
            
            DefV.Frame                 = [];
            DefV.Scale                 = 'mode zscale';
            DefV.CMap                  = [];
            DefV.Colorbar              = 'no';
            DefV.Orient                = [];
            DefV.Rotate                = [];
            DefV.Zoom                  = 1;
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);
            
            % open ds9
            ds9.open;
            
            if (isempty(Frame))
                % use default instead
                Frame = InPar.Frame;
            end
                       

            ImageField = SIM.ImageField;
            IsFits = false;
            
            % prepare cell array of string (FITS images)
            if (ischar(Images))
                [~,List] = Util.files.create_list(Images,NaN);
                IsFits = true;
                Nim = numel(List);
            end
            
            % prepare cell array of matrices
            if (isnumeric(Images))
                List   = {Images};
                IsFits = false;
                Nim = numel(List);
            end
            
            if (SIM.issim(Images))
                %List   = Images;  %.(ImageField);
                IsFits = false;
                Nim = numel(Images);
            end
            
%             if (imCl.isimCl(Images))
%                 %List   = Images;  %.(ImageField);
%                 IsFits = false;
%                 Nim = numel(Images);
%             end
            
            
            % create tmp file names
            IsTmp   = false;
            
            if (~IsFits)
                % create tmp image
                IsTmp     = true;
                for Iim=1:1:Nim
                    TmpName = sprintf('%s.fits',tempname);
                    % create FITS files
                    %fitswrite_my(List{Iim},TmpName);
                    if (SIM.issim(Images))
                        % to fix a problem - delete the PCOUNT and GCOUNT
                        % header keywords
                        %Images(Iim) = delete_key(Images(Iim),{'PCOUNT','GCOUNT','PSCALET1','PSCALET2'});
                        Images(Iim) = delete_key(Images(Iim),{'EXTVER','PCOUNT','GCOUNT','PSCALET1','PSCALET2',...
                            'XTENSION','TTYPE1','TFORM1','TTYPE2','TFORM2','TTYPE3','TFORM3'}); % Na'ama, 20180831, to display fits.fz (fpacked) images
                        
                        
                        Nkey = size(Images(Iim).Header,1);
                        Bl   = cell(Nkey,1);
                        [Bl{1:Nkey}] = deal(' ');
                        Images(Iim).Header = [Images(Iim).Header(:,1:2),  Bl];
                        FITS.write(Images(Iim).(ImageField),TmpName,'Header',Images(Iim).Header);
                    else
                        FITS.write(List{Iim},TmpName);
                    end
                    List{Iim} = TmpName;
                end
            end
           
            if (isempty(Frame))
                Frame = (1:1:Nim).';
            end
            if (~ischar(Frame))
                Frame = Frame(:).*ones(Nim,1);
            end
            
            for Iim=1:1:Nim
                % for each image                
               
                if (ischar(Frame))
                    ds9.load(List{Iim},Frame);
                else
                    ds9.load(List{Iim},Frame(Iim));
                end
                
                % set image properties
                if (~isempty(InPar.Scale))
                    ds9.scale(InPar.Scale);
                end
                if (~isempty(InPar.CMap))
                    ds9.cmap(InPar.CMap);
                end
                if (~isempty(InPar.Colorbar))
                    ds9.colorbar(InPar.Colorbar);
                end
                if (~isempty(InPar.Orient))
                    ds9.orient(InPar.Orient);
                end
                if (~isempty(InPar.Rotate))
                    ds9.rotate(InPar.Rotate);
                end
                if (~isempty(InPar.Zoom))
                    ds9.zoom(InPar.Zoom);
                end
                   
               
            end
            
            if (IsTmp)
                Util.files.delete_cell(List);
            end
            pause(0.2);
            
        end
        
        % save ds9 frame to FITS file
        function read2fits(FileName)
            % Save a ds9 frame to FITS image
            % Package: @ds9
            % Description: Save a ds9 frame to FITS image
            % Input  : - FITS file name to save.
            % Output : null
            % Example: ds9.tofits('try.fits');
            % Reliable: 
            
            % Read entire image from ds9
            Answer  = ds9.system('xpaset -p ds9 save %s',FileName);
        end
        
        % save ds9 frame to SIM image
        function Sim=read2sim
            % Save a ds9 frame in a SIM object
            % Package: @ds9
            % Description: Save a ds9 frame in a SIM object
            % Input  : null.
            % Output : - A SIM object.
            % Eaxmple: S=ds9.tosim;
            % Reliable: 2
            
            % Read entire image from ds9
            TmpName = sprintf('%s.fits',tempname);
            Answer  = ds9.system('xpaset -p ds9 save %s',TmpName);
            %Sim     = image2sim(TmpName);
            Sim     = FITS.read2sim(TmpName);
            delete(TmpName);
        end
    end
    
    % Frame properties methods
    % (scale, cmap, colorbar, orient, pan, rotate, zoom, header)
    methods (Static)
        % Scale image intensity
        function scale(varargin)
            % Set the intensity scale of an image in ds9
            % Package: @ds9
            % Description: Set the intensity scale of an image in ds9
            % Input  : * Arbitrary number of arguments to pass to the ds9
            %            scale command.
            %            Default is 'mode zscale'.
            % Output : null
            % Example: ds9.scale('linear')
            %          ds9.scale('log',100)
            %          ds9.scale('log 100');
            %          ds9.scale('limits',1520,1900)
           % Reliable: 2
            if (nargin==0)
                Mode = 'mode zscale';
            else
                Mode = '';
            end
            
            Mode = ds9.construct_command(Mode,varargin{:});
            
            ds9.system('xpaset -p ds9 scale %s',Mode);
            pause(0.2);
            
        end
        
        % Set image color map
        function cmap(varargin)
            % Set the color map of an image in ds9
            % Package: @ds9
            % Description: Set the color map of an image in ds9
            % Input  : * Arbitrary number of arguments to pass to the ds9
            %            cmap command.
            %            Default is 'invert yes'.
            % Output : null
            % Example: ds9.cmap('Heat')
            %          ds9.cmap Heat
            %          ds9.cmap('value',5)
            %          ds9.cmap('cmap','match');
            % Reliable: 2
            if (nargin==0)
                Mode = 'invert yes';
            else
                Mode = '';
            end
            
            Mode = ds9.construct_command(Mode,varargin{:});
            
            ds9.system('xpaset -p ds9 cmap %s',Mode);
            pause(0.2);

        end
        
        % Set image colorbar
        function colorbar(varargin)
            % set the colorbar of an ds9 image
            % Package: @ds9
            % description: set the colorbar of an ds9 image
            % Input  : * Arbitrary number of arguments to pass to the ds9
            %            colorbar command.
            %            Default is 'no'.
            % Output : null
            % Example: ds9.colorbar('vertical')
            %          ds9.colorbar('size',20)
            % Reliable: 2
            if (nargin==0)
                Mode = 'no';
            else
                Mode = '';
            end
            
            Mode = ds9.construct_command(Mode,varargin{:});
            
            ds9.system('xpaset -p ds9 colorbar %s',Mode);
            pause(0.2);
            
        end
        
        % Set image orientation
        function orient(varargin)
            % Set the x/y orientation of an image in ds9
            % Package: @ds9
            % Description: Set the x/y orientation of an image in ds9
            % Input  : * Arbitrary number of arguments to pass to the ds9
            %            orient command.
            %            Default is 'none'.
            % Output : null
            % Example: ds9.orient('x')
            %          ds9.orient('xy')
            % Reliable: 2
            if (nargin==0)
                Mode = 'none';
            else
                Mode = '';
            end
            
            Mode = ds9.construct_command(Mode,varargin{:});
            
            ds9.system('xpaset -p ds9 orient %s',Mode);
            pause(0.2);
            
        end
        
        % Set image pan
        function pan(varargin)
            % Set the pan (cursor location) of an image in ds9
            % Package: @ds9
            % Description: Set the pan (cursor location) of an image in ds9
            % Input  : * Arbitrary number of arguments to pass to the ds9
            %            pan command.
            %            Default is 'to 0 0 image'.
            % Output : null
            % Example: ds9.pan('200 200 image')
            % Reliable: 2
            if (nargin==0)
                Mode = 'to 0 0 image';
            else
                Mode = '';
            end
            
            Mode = ds9.construct_command(Mode,varargin{:});
            
            ds9.system('xpaset -p ds9 pan %s',Mode);
            pause(0.2);
            
        end
        
        % Set image rotate
        function rotate(varargin)
            % Set the rotation of an image in ds9
            % Package: @ds9
            % Description: Set the rotation of an image in ds9
            % Input  : * Arbitrary number of arguments to pass to the ds9
            %            rotate command.
            %            Default is ''.
            % Output : null
            % Example: ds9.rotate('to',45)
            % Reliable: 2
            if (nargin==0)
                Mode = '';
            else
                Mode = '';
            end
            
            Mode = ds9.construct_command(Mode,varargin{:});
            
            ds9.system('xpaset -p ds9 rotate %s',Mode);
            pause(0.2);
            
        end
        
        % Set image zoom
        function zoom(varargin)
            % Set the zoom of an image in ds9
            % Package: @ds9
            % Description: Set the zoom of an image in ds9
            % Input  : * Arbitrary number of arguments to pass to the ds9
            %            zoom command. If the first argument is a number
            %            than the word "to" will be added (absolute zoom).
            %            Default is 'to to 1'.
            % Output : null
            % Example: ds9.zoom(2,4)
            %          ds9.zoom('in')
            % Reliable: 2
            if (nargin==0)
                Mode = 'to 1';
            else
                if (isnumeric(varargin{1}))
                    Mode = 'to ';
                else
                    Mode = '';
                end
            end
            
            
            Mode = ds9.construct_command(Mode,varargin{:});
            ds9.system('xpaset -p ds9 zoom %s',Mode);
            pause(0.2);
            
        end
        
        % header display
        function header(varargin)
            % Description: Display FITS header dialaog
            % Input  : - Parameters to fits header dialaog.
            %            Default is ''.
            % Output : null
            % Example: ds9.header(1)
            %          ds9.header('save',1,'foo.txt')
            %          ds9.header('close')
            % Reliable: 2
            Mode = '';
            Mode = ds9.construct_command(Mode,varargin{:});
            ds9.system('xpaset -p ds9 header %s',Mode);
            
        end
        
    end
    
    % printing methods
    % (psprint)
    methods (Static)
        % print, psprint
        function psprint(FileName,Resolution,Color,Level)
            % Print a postscript file into current directory
            % Package: @ds9
            % Description: Print a postscript file into current directory.
            % Input  : - PS file name.
            %          - Resolution. Default is 150.
            %          - Color. Default is 'cmyk'.
            %          - PS level. Default is 2.
            % Output : null
            % Reliable: 2
            
            Def.Resolution = 150;
            Def.Color      = 'cmyk';
            Def.Level      = 2;
            if (nargin==1)
                Resolution = Def.Resolution;
                Color      = Def.Color;
                Level      = Def.Level;
            elseif (nargin==2)
                Color      = Def.Color;
                Level      = Def.Level;
            elseif (nargin==3)
                Level      = Def.Level;
            elseif (nargin==4)
                % do nothing
            else
                error('Illegal number of input arguments: ds9.print(FileName,[Resolution,Color,Level])');
            end
            FileName = sprintf('%s%s%s',pwd,filesep,FileName);  % to print in current directory
            ds9.system('xpaset -p ds9 print destination file');
            ds9.system('xpaset -p ds9 print filename %s',FileName);
            ds9.system('xpaset -p ds9 print resolution %d',Resolution);
            ds9.system('xpaset -p ds9 print color %s',Color);
            ds9.system('xpaset -p ds9 print level %d',Level);
            ds9.system('xpaset -p ds9 print');
            
        end
        
    end
    
    % Region methods
    % (write_region, load_region, delete_region, save_region, plot, text)
    methods (Static)
        % construct a region file
        function FileName=write_region(Cat,varargin)
            % Write a regions file for a list of coordinates and properties
            % Package: @ds9
            % Description: Write a regions file for a list of coordinates
            %              and properties.
            % Input  : - A two column matrix of coordinates [X,Y] or
            %            [RA,Dec] (degrees), or an AstCat object.
            %            The AstCat object coordinates column names are
            %            defined by the 'ColName...' argumenents.
            %            For RA/Dec units must be degrees.
            %          * Arbitrary number of ...,key,val,... pairs.
            %            The following keywords are available:
            %            'FileName' - Output region file name.
            %                         Default is tempname.
            %            'Append'   - Append region file to an existing
            %                         region file. Default is false.
            %            'Coo'      - Coordinates type: 'image'|'fk5'.
            %                         Default is 'image' (i.e., pixels).
            %            'Units'    - If 'Coo' is 'fk5' than this specify
            %                         if the input coordinates are in 'deg'
            %                         or 'rad'. Default is 'deg'.
            %            'Marker'   - A string or a cell array of strings
            %                         of markers type.
            %                         Options are: 'circle'|'circ'|'o'
            %                                      'box'|'s'
            %                                      'ellipse'|'e'
            %                                      'line'|'l'
            %                                      'vector'|'v'
            %                                      'polygon'|'p'
            %                         Default is 'circle'.
            %             'Color'   - A string or a cell array of strings
            %                         of marker colors
            %                         ('red'|'blue'|'green'|'black'|
            %                          'white'|...).
            %                         Default is 'red'.
            %             'Width'   - A scalar or a vector of markers
            %                         width. Default is 1.
            %             'Size'    - Size array. Either one row, or a row 
            %                         per coordinate. Default is 10.
            %                         Note that the required number of
            %                         columns in size depands on the marker
            %                         type.
            %                         Following attributes are needed:
            %                           For 'circle':  [Radius]
            %                           For 'box':     [Width,Height]
            %                           For 'ellipse': [Major axis, Minor axis]
            %                           For 'vector':  [Length,PA]
            %                           For 'line':    [StartX,StartY,EndX,EndY]
            %                           For 'polygon': The coordinates are
            %                           used as the polygon verteces.
            %                         Default is 10.
            %             'Text'    - Text to plot with marker.
            %                         Default is ''.
            %             'Font'    - Font type. Default is 'helvetica'.
            %             'FontSize'- Font size. Default is 16.
            %             'FontStyle'-Font style. Default is 'normal'.
            %             'ColNameX'- Cell array of possible column names
            %                         of X coordinate in AstCat object.
            %                         If the first input in an AstCat
            %                         object then the first exitsing column
            %                         name will be used.
            %                         Default is
            %                         {'XWIN_IMAGE','X','X_IMAGE'}.
            %             'ColNameY'- Like 'ColNameX', but for the Y
            %                         coordinate.
            %                         Default is
            %                         {'YWIN_IMAGE','Y','Y_IMAGE'}.
            %             'ColNameRA'-Like 'ColNameX', but for the RA
            %                         coordinate.
            %                         Default is
            %                         {'ALPHAWIN_IMAGE','RA'}.
            %             'ColNameDec'-Like 'ColNameX', but for the Dec
            %                         coordinate.
            %                         Default is
            %                         {'DELTAWIN_IMAGE','DEC'}.
            % Output : - Region file name.
            % See also: ds9.plot
            % Example: FileName=write_region(Cat);
            %          write_region([X,Y],'Marker','s','Color','cyan')
            % Reliable: 2
            
            DefV.FileName        = tempname;  % use temp file name
            DefV.Append          = false;
            DefV.Coo             = 'image';   % 'image'|'fk5'
            DefV.Units           = 'deg';     % for 'image' this is always pix!
            DefV.Marker          = 'circle';  % 'circle'|'box'|...
            DefV.Color           = 'red';
            DefV.Width           = 1;
            DefV.Size            = 10;
            DefV.Text            = '';
            DefV.Font            = 'helvetica';  %'helvetica 16 normal'
            DefV.FontSize        = 16;
            DefV.FontStyle       = 'normal';
            DefV.ColNameX        = {'XWIN_IMAGE','X','X_IMAGE'};
            DefV.ColNameY        = {'YWIN_IMAGE','Y','Y_IMAGE'};
            DefV.ColNameRA       = {'ALPHAWIN_J2000','RA'};
            DefV.ColNameDec      = {'DELTAWIN_J2000','Dec'};
            
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            % check if region file exist
            if (exist(InPar.FileName,'file')==0)
               if (InPar.Append)
                   error('User requested to append region file, but file doesnt exist');
               end
            end
            
            
            IsXY = false;
            switch lower(InPar.Coo)
             case 'image'
                CooUnits    = '';
                IsXY        = true;
                InPar.Units = 'pix';
             case 'fk5'
                CooUnits = '"';
             otherwise
                error('Coo units is not supported');
            end
            
            % prepare catalog
            if (AstCat.isastcat(Cat))
                % read the coordinates from an AstCat object
                if (IsXY)
                    ColNameX   = select_exist_colnames(Cat,InPar.ColNameX(:));
                    ColNameY   = select_exist_colnames(Cat,InPar.ColNameY(:));
                else
                    ColNameX   = select_exist_colnames(Cat,InPar.ColNameRA(:));
                    ColNameY   = select_exist_colnames(Cat,InPar.ColNameDec(:));
                end
                
                X = col_get(Cat,ColNameX);
                Y = col_get(Cat,ColNameY);
            else
                % Assume Cat is a two column matrix [X,Y]
                X = Cat(:,1);
                Y = Cat(:,2);
            end
            
            % In case of spherical coordinates - convert to deg
            if (~IsXY)
                Factor = convert.angular(InPar.Units,'deg');
                X      = X.*Factor;
                Y      = Y.*Factor;
            end
                
            
            % Number of regions to plot
            Nreg    = numel(X);
            % Prep properties
            if (~iscell(InPar.Marker))
                InPar.Marker = {InPar.Marker};
            end
            Nmarker = numel(InPar.Marker);
            if (~iscell(InPar.Color))
                InPar.Color = {InPar.Color};
            end
            Ncolor = numel(InPar.Color);
            Nwidth = numel(InPar.Width);
            Nsize  = size(InPar.Size,1);
            if (~iscell(InPar.Text))
                InPar.Text = {InPar.Text};
            end
            Ntext = numel(InPar.Text);
            if (~iscell(InPar.Font))
                InPar.Font = {InPar.Font};
            end
            Nfont = numel(InPar.Font);
            Nfontsize = numel(InPar.FontSize);
            if (~iscell(InPar.FontStyle))
                InPar.FontStyle = {InPar.FontStyle};
            end
            Nfontstyle = numel(InPar.FontStyle);
            
            % Open file
            if (InPar.Append)
                % append header to an existing region file
                FID = fopen(InPar.FileName,'a');
            else
                FID = fopen(InPar.FileName,'w');
                fprintf(FID,'# Region file format: DS9 version 4.1\n');
                fprintf(FID,'# Written by Eran Ofek via ds9.write_region.m\n');
                fprintf(FID,'global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n');
                fprintf(FID,'%s\n',InPar.Coo);
            end

            % for each coordinate (refion)
            for Ireg=1:1:Nreg
                
                switch lower(InPar.Marker{min(Ireg,Nmarker)})
                    case {'circle','circ','o'}
                        fprintf(FID,'%s(%15.8f,%15.8f,%15.8f%s)',... 
                                    'circle', X(Ireg), Y(Ireg),...
                                    InPar.Size(min(Ireg,Nsize)),CooUnits);
                    case {'box','b','s'}
                        if (numel(InPar.Size)==1)
                            fprintf(FID,'%s(%15.8f,%15.8f,%15.8f%s,%15.8f%s,%9.5f)',... 
                                        'box',X(Ireg),Y(Ireg),...
                                        InPar.Size(min(Ireg,Nsize),1),CooUnits,...
                                        InPar.Size(min(Ireg,Nsize),1),CooUnits,...
                                        0);
                        else
                            fprintf(FID,'%s(%15.8f,%15.8f,%15.8f%s,%15.8f%s,%9.5f)',... 
                                        'box',X(Ireg),Y(Ireg),...
                                        InPar.Size(min(Ireg,Nsize),1),CooUnits,...
                                        InPar.Size(min(Ireg,Nsize),2),CooUnits,...
                                        InPar.Size(min(Ireg,Nsize),3));
                        end
                    case {'ellipse','e'}
                        fprintf(FID,'%s(%15.8f,%15.8f,%15.8f%s,%15.8f%s,%9.5f)',... 
                                    'ellipse',X(Ireg),Y(Ireg),...
                                    InPar.Size(min(Ireg,Nsize),1),CooUnits,...
                                    InPar.Size(min(Ireg,Nsize),2),CooUnits,...
                                    InPar.Size(min(Ireg,Nsize),3));
                    case {'vector','v'}
                        fprintf(FID,'# %s(%15.8f,%15.8f,%15.8f%s,%9.5f)',... 
                                    'vector',X(Ireg),Y(Ireg),...
                                    InPar.Size(min(Ireg,Nsize),1),CooUnits,...
                                    InPar.Size(min(Ireg,Nsize),2));
    
                    case {'line','l'}
                        fprintf(FID,'%s(%15.8f,%15.8f,%15.8f,%15.8f)',... 
                                    'line',...
                                    InPar.Size(min(Ireg,Nsize),1),...
                                    InPar.Size(min(Ireg,Nsize),2),...
                                    InPar.Size(min(Ireg,Nsize),3),...
                                    InPar.Size(min(Ireg,Nsize),4));
                                
                    case {'polygon','p'}
                        if (Ireg==1)
                            fprintf(FID,'%s(','polygon');
                        end
                        fprintf(FID,'%15.8f,%15.8f,',X(Ireg),Y(Ireg)); %Par.Size(Ireg,Isize));
                        if (Ireg==Nreg)
                            fprintf(FID,'%15.8f,%15.8f)',X(Ireg),Y(Ireg)); %Par.Size(Ireg,Nsize));
                        end
                    otherwise
                        error('Unknown Marker option');
                end
                
                % additional properties
                fprintf(FID,'# color=%s width=%d font="%s %d %s" text={%s}\n',...
                            InPar.Color{min(Ireg,Ncolor)},...
                            InPar.Width(min(Ireg,Nwidth)),...
                            InPar.Font{min(Ireg,Nfont)},...
                            InPar.FontSize(min(Ireg,Nfontsize)),...
                            InPar.FontStyle{min(Ireg,Nfontstyle)},...
                            InPar.Text{min(Ireg,Ntext)});
                        
            end
            fclose(FID); % close region file

            FileName = InPar.FileName;
            pause(0.2);
            
        end
        
        % load regions from file
        function load_region(FileName)
            % load regions file name into current ds9 frame
            % Package: @ds9
            % Description: load regions file name into current ds9 frame.
            % Input  : - FileName. If file name does not contains the
            %            full path then assume the file is in the current
            %            directory.
            % Output : null
            % Example: ds9.load_region('A.reg')
            % Reliable: 2
            
            % check if FileName contains directory name
            if (isempty(strfind(FileName,filesep)))
                % File Name without dir
                % add pwd to file name
                FileName = sprintf('%s%s%s',pwd,filesep,FileName);
            end
            
            ds9.system('xpaset -p ds9 regions load %s',FileName);
        end
        
        % delete regions from frame
        function delete_region
            % Delete all regions from ds9 frame
            % Package: @ds9
            % Description: Delete all regions from ds9 frame
            % Input  : null
            % Output : null
            % Example: ds9.delete_region
            % Reliable: 2
            
            ds9.system('xpaset -p ds9 regions delete all');
        end
        
        % save regions to file
        function save_region(FileName)
            % Save regions into a file
            % Package: @ds9
            % Description: Save regions into a file.
            % Input  : FileName. If file name does not contains the
            %            full path then save it in the current
            %            directory. Default is 'ds9.reg'.
            % Output : null
            % Example: ds9.save_region
            % Reliable: 2
            
            if (nargin==0)
                FileName = 'ds9.reg';
            end
            
            % check if FileName contains directory name
            if (isempty(strfind(FileName,filesep)))
                % File Name without dir
                % add pwd to file name
                FileName = sprintf('%s%s%s',pwd,filesep,FileName);
            end
            ds9.system('xpaset -p ds9 regions save %s',FileName);
        end
        
        function plotc(varargin)
            % Generate and plot a region file from a list of celestial coordinates
            % Package: @ds9
            % Input  : * see ds9.plot(...,'Coo','fk5')
            
            ds9.plot(varargin{:},'Coo','fk5');
            
        end
        
        % plot regions
        function plot(varargin)
            % Generate and plot a region file from a list of coordinates
            % Package: @ds9
            % Description: Generate and plot a region file from a list of
            %              coordinates.
            %              This is like the function ds9.write_region with
            %              the exception that it also load the region file
            %              into the current frame.
            %              By default this program deletes the region file
            %              name after plotting. Use ds9.write_region to
            %              create the file.
            % Input  : * Parameters to pass to ds9.write_region.m.
            %            See ds9.write_region.m for details.
            %            Alternatively, the second argument may be a string
            %            containing the marker and color
            %            (e.g., 'yo' - yellow circle; 'ks' - black square).
            %            Alternatively, the first argument is X coordinate,
            %            the second is Y coordinate and the rest as before.
            % Output : null
            % Example: ds9.plot([1000 1000]);
            %          ds9.plot([1020 1010],'Marker','s','Size',[10 10 45]);
            %          ds9.plot([900 900],'rs');
            %          ds9.plot(950,950);
            %          ds9.plot(950,950,'rs','size',[20 20 0]);
            %          ds9.plot(1000,1020,'ro','Text','Object X');
            %          ds9.plot(Sim);   % plot circles around AstCat object
            % Reliable: 2
            DeleteFile = true;
            
            Narg = numel(varargin);
            if (Narg>1)
                if (isnumeric(varargin{2}))
                    % treating ds9.plot(1,1,'o',...) like input
                    X = varargin{1};
                    Y = varargin{2};
                    varargin = varargin(2:end);
                    varargin{1} = [X,Y];
                    Narg = Narg - 1;
                end
                
                if ((Narg.*0.5)==floor(Narg.*0.5))
                    % treating ds9.plot([1 1],'rs',...) like input
                    
                    Size = 10;
                    if (~isempty(strfind(varargin{2},'o')))
                        Marker = 'o';
                    elseif (~isempty(strfind(varargin{2},'s')))
                        Marker = 's';
                        Size   = [10 10 0];
                    else
                        Marker = 'o';
                    end
                    if (~isempty(strfind(varargin{2},'r')))
                        Color = 'red';
                    elseif (~isempty(strfind(varargin{2},'b')))
                        Color = 'blue';
                    elseif (~isempty(strfind(varargin{2},'g')))
                        Color = 'green';
                    elseif (~isempty(strfind(varargin{2},'k')))
                        Color = 'black';
                    elseif (~isempty(strfind(varargin{2},'w')))
                        Color = 'white';
                    elseif (~isempty(strfind(varargin{2},'m')))
                        Color = 'magenta';
                    elseif (~isempty(strfind(varargin{2},'c')))
                        Color = 'cyan';  
                    elseif (~isempty(strfind(varargin{2},'y')))
                        Color = 'yellow'; 
                    else
                        Color = 'red';
                    end
                    varargin = varargin([1, 3:end]);
                    if any(strcmpi(varargin,'size'))
                        % Size already in varargin
                        varargin = [varargin, {'Marker',Marker,'Color',Color}];
                    else
                        % Size is not in varargin - use default
                        varargin = [varargin, {'Marker',Marker,'Color',Color,'Size',Size}];
                    end
                end
            end
                
            FileName = ds9.write_region(varargin{:});
            ds9.load_region(FileName);
            if (DeleteFile)
                delete(FileName);
                FileName = [];
            end
            
        end
                
        % plot line by x/y coordinates
        function line_xy(X,Y,varargin)
            % Plot a broken line
            % Package: @ds9
            % Description: Plot a broken line (curve) in ds9 image.
            % Input  : - Vector of X.
            %          - Vector of Y.
            %          * Additional parameters to pass to ds9.write_region
            %            function. The first argument may be a plot-like
            %            color indicator.
            % Output : null
            % Example: ds9.line_xy([356 400],[2000 2100],'Color','green')
            %          ds9.line_xy([356 400],[2000 2100],'r')
            %          ds9.line_xy([356 400],[2000 2100],'w','width',3)
            %          ds9.line_xy([360,400,500],[2000 2010 2100])
            % Reliable: 2

            
            DeleteFile = true;
            
            Narg  = numel(varargin);
            Color = 'red';
            if ~(Narg.*0.5==floor(Narg.*0.5))
                % assume 3rd parameter is a color indicator
                
                if (~isempty(strfind(varargin{1},'r')))
                    Color = 'red';
                elseif (~isempty(strfind(varargin{1},'b')))
                    Color = 'blue';
                elseif (~isempty(strfind(varargin{1},'g')))
                    Color = 'green';
                elseif (~isempty(strfind(varargin{1},'k')))
                    Color = 'black';
                elseif (~isempty(strfind(varargin{1},'w')))
                    Color = 'white';
                elseif (~isempty(strfind(varargin{1},'m')))
                    Color = 'magenta';
                elseif (~isempty(strfind(varargin{1},'c')))
                    Color = 'cyan';  
                elseif (~isempty(strfind(varargin{1},'y')))
                    Color = 'yellow'; 
                else
                    Color = 'red';
                end
                varargin = varargin(2:end);
            end
            
            X1 = X(1:end-1);
            Y1 = Y(1:end-1);
            X2 = X(2:end);
            Y2 = Y(2:end);
            
            FileName = ds9.write_region([X1(:), Y1(:)],'Color',Color,varargin{:},'Marker','line','Size',[X1(:), Y1(:), X2(:), Y2(:)]);
            ds9.load_region(FileName);
            if (DeleteFile)
                delete(FileName);
                FileName = [];
            end
            
        end
        
        function line_lt(X,Y,Length,Theta,varargin)
            % Plot multiple lines based on X,Y,length,theta
            % Package: @ds9
            % Description: Plot multiple lines based on X,Y,length,theta
            % Input  : - Vector of X starting points.
            %          - Vector of Y starting points.
            %          - Vector of line lengths.
            %          - Vector of line angles [deg] measured from X-axis.
            %          * Additional parameters to pass to ds9.write_region
            %            function. The first argument may be a plot-like
            %            color indicator.
            % Output : null
            % Example: ds9.line_lt(400,1900,100,10)
            %          ds9.line_lt(400,1900,100,10,'g','Width',3)
            %          ds9.line_lt([400,450],[1900,1910],[100,200],[10,45])
            
            Xend = X + Length.*cosd(Theta);
            Yend = Y + Length.*sind(Theta);
            
            N = numel(X);
            for I=1:1:N
                % for each line
                ds9.line_xy([X(I), Xend(I)],[Y(I), Yend(I)],varargin{:});
            end
        end
        
        % plot text
        function text(X,Y,Text,varargin)
            % plot text to ds9 current frame
            % Package: @ds9
            % Description: plot text to ds9 current frame in image
            %              coordinates position.
            % Input  : - Image coordinate X position (pixels), or J2000.0
            %            R.A. (deg or sexagesimal string). If string then
            %            will set 'Coo' to 'fk5' and assumes that also Dec
            %            is given in sexagesimal format.
            %          - Image coordinate Y position (pixels), or J2000.0
            %            Dec. (deg or sexagesimal string).
            %          - Text.
            %          * Arbitrary number of ...,key,val,... pairs.
            %            The following keywords are available:
            %             'Color'   - A string or a cell array of strings
            %                         of marker colors
            %                         ('red'|'blue'|'green'|'black'|
            %                          'white'|...).
            %                         Default is 'red'.
            %             'Width'   - A scalar or a vector of markers
            %                         width. Default is 1.
            %             'Font'    - Font type. Default is 'helvetica'.
            %             'FontSize'- Font size. Default is 16.
            %             'FontStyle'-Font style. Default is 'normal'.
            %             'Coo'     - Coordinate system: 'image'|'fk5'.
            %                         Default is 'image'.
            %                         
            % Output : null
            % see also: ds9.plot, ds9.write_region
            % Example: ds9.text(700,900,'Hello');
            %          ds9.text(700,980.1,'Hello','color','cyan');
            %          ds9.text(149.41455,69.758364,'Star','Coo','fk5'); % deg
            %          ds9.text('9:56:31.559','+69:49:27.60','Object','Color','blue')
            % Reliable: 2
            
            
                
            
            DefV.Color               = 'red';
            DefV.Width               = 1;
            DefV.Font                = 'helvetica';  %'helvetica 16 normal'
            DefV.FontSize            = 16;
            DefV.FontStyle           = 'normal';
            DefV.Coo                 = 'image';
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);
            
            if (ischar(X))
                InPar.Coo = 'fk5';
                IsSexagesimal = true;
            else
                IsSexagesimal = false;
            end
            
            if (IsSexagesimal)
                ds9.system('echo "%s; text %s %s # color=%s width=%d font={%s %d %s} text={%s}" | xpaset ds9 regions',...
                                InPar.Coo,X,Y,InPar.Color,InPar.Width,InPar.Font,InPar.FontSize,InPar.FontStyle,Text);
            else
                
                ds9.system('echo "%s; text %f %f # color=%s width=%d font={%s %d %s} text={%s}" | xpaset ds9 regions',...
                            InPar.Coo,X,Y,InPar.Color,InPar.Width,InPar.Font,InPar.FontSize,InPar.FontStyle,Text);
            end
        
            pause(0.2);
            
        end
        
        
    end
    
    % Tile methods
    methods (Static)
        % Set the tile display mode
        function tile(Par)
            % Set the tile the display mode of ds9
            % Package: @ds9
            % Description: Set the tile the display mode of ds9.
            % Input  : - Either a string to pass to the tile command,
            %            or this is a vector of the number of columns and
            %            rows in the tile display (i.e., layout command)
            %            Alternatively a third element may be provided with
            %            the grid gap (in pixels).
            %            Default is [3,2].
            % Output : null
            % Example: ds9.tile([3 3 10])
            %          ds9.tile('grid direction x')
            
            if (nargin==0)
                Par = [3 2];
            end
            
            if (ischar(Par))
                ds9.system('xpaset -p ds9 tile %s',Par);
            else
                if (numel(Par)==2)
                    ds9.system('xpaset -p ds9 tile grid layout %d %d',Par(1),Par(2));
                elseif (numel(Par)==3)
                    ds9.system('xpaset -p ds9 tile grid layout %d %d',Par(1),Par(2));
                    ds9.system('xpaset -p ds9 tile grid gap %d',Par(3));
                else
                    error('Numeric parameter may contain 2 or 3 values [Col, Rows, Gaps]');
                end
            end
            ds9.system('xpaset -p ds9 tile');
            
        end
        
        function single
            % Set to single image display mode
            % Package: @ds9
            % Input  : null
            % Output : null
            % Example: ds9.single
            % Reliable: 2
            
            ds9.system('xpaset -p ds9 single');
            
        end
        
        function blink(Interval)
            % Set to blink display mode
            % Package: @ds9
            % Input  : - Time interval. Default is 0.5s.
            % Output : null
            % Example: ds9.blink
            % Reliable: 2
        
            if (nargin==0)
                Interval = 0.5;
            end
            
            ds9.system('xpaset -p ds9 blink interval %f',Interval);
            ds9.system('xpaset -p ds9 blink');
            fprintf('To stop blink mode - execute: ds9.single\n');
            
        end
        
        % array
        
        % bin
        % block
        
        % blink
        
        % crop
        
        % crosshair
        
        % cursor
        
        % get data
        
        % Set grid
        
        % dssstsci
        function dss(RA,Dec,Size,Band,Frame,Save)
            % Get a DSS sky image from stsci
            % Package: @ds9
            % Description: Get a DSS sky image from stsci.
            % Input  : - J2000 RA [radians, [H M S], or sexagesimal string],
            %            or object name (object name if Dec is empty).
            %          - J2000 Dec [radians, [sign D M S], or sexagesimal
            %            string]. If empty then assume that the first
            %            argument is object name.
            %          - Image size in arcsec. Default is [900 900].
            %          - Band: 'red2'|'blue2'|'ir2'|'red1'|'blue1'|...
            %                  'quickv'|'gsc1'|'gsc2'.
            %            Default is 'red2'.
            %          - Frame number: 'new'|'current'. Default is 'new'.
            %          - Save. Defaukt is false.
            
            Def.Dec   = [];
            Def.Size  = [900 900];
            Def.Band  = 'red2';
            Def.Frame = 'new';
            Def.Save  = false;
            if (nargin==1)
                Dec     = Def.Dec;
                Size    = Def.Size;
                Band    = Def.Band;
                Frame   = Def.Frame;
                Save    = Def.Save;
            elseif (nargin==2)
                Size    = Def.Size;
                Band    = Def.Band;
                Frame   = Def.Frame;
                Save    = Def.Save;
            elseif (nargin==3)
                Band    = Def.Band;
                Frame   = Def.Frame;
                Save    = Def.Save;
            elseif (nargin==4)
                Frame   = Def.Frame;
                Save    = Def.Save;
            elseif (nargin==5)
                Save    = Def.Save;
            elseif (nargin==6)
                % do nothing
            else
                error('Illegal number of input arguments: dss(RA,[Dec,Size,Band,Frame,Save]');
            end
                
            if (isempty(Dec))
                % Assume RA contains object name
                ds9.system('xpaset -p ds9 dssstsci %s',RA);
            else
                RA  = celestial.coo.convertdms(RA,'gH','SH');
                Dec = celestial.coo.convertdms(Dec,'gD','SD');
                ds9.system('xpaset -p ds9 dssstsci coord %s %s sexagesimal',RA,Dec);
            end
            
            if (numel(Size)==1)
                Size = [Size Size];
            end
            ds9.system('xpaset -p ds9 dssstsci size %d %d arcsec');
            
            if (Save)
                ds9.system('xpaset -p ds9 dssstsci save yes');
            else
                ds9.system('xpaset -p ds9 dssstsci save no');
            end
            
            % frame
            ds9.system('xpaset -p ds9 dssstsci frame %s',Frame);
            
            % band
            Translation = {'red2','poss2ukstu_red';...
                           'ir2','poss2ukstu_ir';...
                           'blue2','poss2ukstu_blue';...
                           'red1','poss1_red';...
                           'blue1','poss1_blue';...
                           'quickv','quickv';...
                           'gsc1','phase2_gsc1';...
                           'gsc2','phase2_gsc2'};
            Iband = find(strcmp(Translation(:,1),Band));
            if (isempty(Iband))
                error('Illegal band name');
            end
            Band = Translation{Iband,2};
            ds9.system('xpaset -p ds9 dssstsci survey %s',Band);
            
            % get image
            ds9.system('xpaset -p ds9 dssstsci');
            
        end
        
        % dsssao
        % dsseso
        % dssstsci
        % nvss
        % vlas
        % vo
        % 2mass
        % skyview
        % url
        
        % catalog
        
        
        % header
        
        % iexem
        
    end % methods
        
    % lock and match methods
    % (lock_wcs, lock_xy, match_wcs, match_xy, match_scale, match_scalelimits)
    methods (Static)
            
        % lock by wcs coordinates
        function lock_wcs(Par)
            % Lock all images WCS to current frame
            % Package: @ds9
            % Description: Lock all images WCS to current frame.
            % Input  : - true|false. If not given than toggle.
            % Output : null
           
            if (nargin==0)
                % get lock status
                Ans = ds9.system('xpaget ds9 lock frame wcs');
                switch lower(Ans(1:3))
                    case 'wcs'
                        % toggle off
                        Par = false;
                    otherwise
                        Par = true;
                end
            end
            %ds9.system('xpaset -p ds9 lock frame wcs');
            if (Par)
                ds9.system('xpaset -p ds9 lock frame wcs');
            else
                ds9.system('xpaset -p ds9 lock frame none');
            end
                        
        end
        
        % lock by image coordinates
        function lock_xy(Par)
            % Lock all images x/y coordinayes to current frame
            % Package: @ds9
            % Description: Lock all images x/y coordinayes to current frame.
            % Input  : - true|false. If not given than toggle.
            % Output : null
           
            if (nargin==0)
                % get lock status
                Ans = ds9.system('xpaget ds9 lock frame image');
                switch lower(Ans(1:4))
                    case 'imag'
                        % toggle off
                        Par = false;
                    
                    otherwise
                        Par = true;
                        %error('Unknonw Answer from xpaget lock frame image');
                end
            end
            %ds9.system('xpaset -p ds9 lock frame image');
            if (Par)
                ds9.system('xpaset -p ds9 lock frame image');
            else
                ds9.system('xpaset -p ds9 lock frame none');
            end
                        
        end
        
        
        % match by wcs coordinates
        function match_wcs
            % Match the WCS coordinates of all frames to the current frame
            % Package: @ds9
            % Description: Match the WCS coordinates of all frames to the
            %              current frame.
            % Input  : null
            % Output : null
            
            ds9.system('xpaset -p ds9 match frame wcs');
        end
        
        % match by image coordinates
        function match_xy
            % Match the image coordinates of all frames to the current frame
            % Package: @ds9
            % Description: Match the image coordinates of all frames to the
            %              current frame.
            % Input  : null
            % Output : null
            
            ds9.system('xpaset -p ds9 match frame image');
        end
 
        % match by intensity scale
        function match_scale
            % Match the intensity scale of all frames to the current frame
            % Package: @ds9
            % Description: Match the intensity scale of all frames to the
            %              current frame.
            % Input  : null
            % Output : null
            
            ds9.system('xpaset -p ds9 match scale');
        end
            
        % match by intensity scale
        function match_scalelimits
            % Match the intensity scalelimits of all frames to the current frame
            % Package: @ds9
            % Description: Match the intensity scalelimits of all frames to
            %              the current frame.
            % Input  : null
            % Output : null
            
            ds9.system('xpaset -p ds9 match scalelimits');
        end
        
        % match by intensity scale
        function match_colorbar
            % Match the intensity colorbar of all frames to the current frame
            % Package: @ds9
            % Description: Match the intensity colorbar of all frames to
            %              the current frame.
            % Input  : null
            % Output : null
            
            ds9.system('xpaset -p ds9 match colorbar');
        end

        
    end % methods
        
    
    % xy2coo, coo2xy
    methods (Static)
        function [CooX,CooY]=coo2xy(RA,Dec,CooType)
            % Convert RA/Dec to X/Y (image) using ds9 tools
            % Package: @ds9
            % Input  : - RA [deg].
            %          - Dec [deg].
            %          - Output coordinate type: {'image'|'physical'}.
            %            Default is 'image.
            % Output : - X [pix].
            %          - Y [pix].
           
            if (nargin<3)
                CooType = 'image';
            end
            N    = numel(RA);
            CooX = zeros(N,1);
            CooY = zeros(N,1);
            for I=1:1:N
                %--- set coordinates of crosshair ---
                ds9.system(sprintf('xpaset -p ds9 crosshair %f %f wcs icrs',RA(I),Dec(I)));
                %--- get Coordinates of crosshair ---
                CooIm = ds9.system(sprintf('xpaget ds9 crosshair %s',CooType));
                %[CooX(I), CooY(I)] = strread(CooIm,'%f %f',1); %,'headerlines',4);
                Cell = textscan(CooIm,'%f %f'); %,'headerlines',4);
                CooX(I) = Cell{1};
                CooY(I) = Cell{2};
            end
            ds9.system(sprintf('xpaset -p ds9 mode pointer'));

        end
        
    end 
    
    % Interactive control and information (mouse)
    % (ginput, getpos, getcoo, getbox)
    methods (Static)
        
        function [CooX,CooY,Val,Key]=ginput(CooType,ExitMode,Mode)
            % Interactively get the coordinates (X/Y or WCS)
            % Package: @ds9
            % Description: Interactively get the coordinates (X/Y or WCS) and value
            %              of the pixel selected by the mouse (left click) or by clicking
            %              any character on the ds9 display.
            % Input  : - Coordinate type {'fk4'|'fk5'|'icrs'|'eq'|
            %                             'galactic'|'image'|'physical'},
            %            default is 'image'.
            %            'eq' is equivalent to 'icrs'.
            %          - Operation mode:
            %            If numeric than will return after the user clicked
            %            the specified number of times.
            %            if 'q' than will return if the user clicked 'q'.
            %            Default is 'q'.
            %          - Operation mode:
            %            'any'   - will return after any character or left click is
            %                      pressed (default).
            %            'key'   - will return after any character is pressed.
            %            'mouse' - will return after mouse left click is pressed.
            % Output : - X/RA or Galactic longitude. If celestial coordinates, the
            %            return the result in radians.
            %          - Y/Dec or Galactic latitude. If celestial coordinates, the
            %            return the result in radians.
            %          - Pixel value.
            %          - Cell array of string of clicked events.
            %            '<1>' for mouse left click.
            % Required: XPA - http://hea-www.harvard.edu/saord/xpa/
            % Tested : Matlab 7.0
            %     By : Eran O. Ofek                    Feb 2007
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [X,Y,V]=ds9.ginput(3,'fk5');   % return the WCS RA/Dec position
            %          [X,Y,V,Key] = ds9.ginput('icrs','q','key');
            %          [X,Y,V,Key] = ds9.ginput('icrs',10,'mouse');
            %          [X,Y,V,Key] = ds9.ginout('image',2,'any');
            % Reliable: 2
            
            
            RAD = 180./pi;

            Def.CooType  = 'image';
            Def.ExitMode = 'q';
            Def.Mode     = 'any';
            if (nargin==0)
               CooType  = Def.CooType;
               ExitMode = Def.ExitMode;
               Mode     = Def.Mode;
            elseif (nargin==1)
               ExitMode = Def.ExitMode;
               Mode     = Def.Mode;
            elseif (nargin==2)
               Mode     = Def.Mode;
            elseif (nargin==3)
               % do nothing
            else
               error('Illegal number of input arguments: ds9.ginput(CooType,Mode);');
            end

            switch lower(Mode)
             case 'mouse'
                Mode = '';
            end
            
            switch lower(CooType)
             case {'eq','equatorial'}
                CooType = 'icrs';
             otherwise
                % do nothing
            end

            switch lower(CooType)
             case {'fk4','fk5','icrs'}
                String = sprintf('wcs %s degrees',lower(CooType));
             case {'gal','galactic'}
                String = 'wcs galactic degrees';
             case {'image'}
                String = 'image';
             case {'physical'}
                String = 'physical';
             otherwise
                error('Unknown CooType option');
            end
            % remove "degrees" from coordinate string for crosshair command
            CrosshairString = strsplit(String, ' degrees');
            CrosshairString = CrosshairString{1};

%             CooX  = zeros(N,1);
%             CooY  = zeros(N,1);
%             Val   = zeros(N,1);
%             Key   = cell(N,1);
            Cont  = true;
            I     = 0;
            while Cont
               I = I + 1;
               %--- get Coordinates (from interactive mouse click) ---
               [Coo] = ds9.system('xpaget ds9 imexam %s coordinate %s',Mode,String);
               SpCoo = regexp(Coo,' ','split');
               Key{I} = SpCoo{1};
               CooX(I)   = str2double(SpCoo{2});
               CooY(I)   = str2double(SpCoo{3});

               %--- set crosshair position to Coordinates ---
               ds9.system('xpaset -p ds9 crosshair %d %d %s',CooX(I),CooY(I),CrosshairString);
               %--- get Coordinates of crosshair ---
               [CooIm] = ds9.system('xpaget ds9 crosshair image');

               %--- get Pixel value at crosshair position ---
               [ValStr] = ds9.system('xpaget ds9 data image %d %d 1 1 yes',CooX(I),CooY(I));
               ValStr = strtrim(ValStr);
               
               %--- Exit crosshair mode ---
               ds9.system('xpaset -p ds9 mode none');

               if (isempty(ValStr))
                    Val(I) = NaN;
               else
                    Val(I)          = sscanf(ValStr,'%f');
               end
               %CooVal       = sscanf(Coo,'%f %f');
               switch lower(CooType)
                case {'fk4','fk5','icrs','gal','galactic'}
                   CooX(I)            = CooX(I)./RAD;
                   CooY(I)            = CooY(I)./RAD;
                otherwise
                   % do nothing - image coordinates
               end

               if (isnumeric(ExitMode))
                   if (I>=ExitMode)
                       Cont = false;
                   end
               else
                   switch lower(ExitMode)
                       case 'q'
                           if (strcmpi(Key{I},'q'))
                               Cont = false;
                           end
                       otherwise
                           error('Unknown Mode option');
                   end
               end
            end
            
            if (isnumeric(ExitMode))
                CooX = CooX.';
                CooY = CooY.';
                Val  = Val.';
                Key  = Key.';
            else
                CooX = CooX(1:end-1).';
                CooY = CooY(1:end-1).';
                Val = Val(1:end-1).';
                Key = Key(1:end-1).';
            end
            
            ds9.system('xpaset -p ds9 mode pointer');

        end
        
        function [CooX,CooY,Value,Key]=getpos(varargin)
            % Get X,Y position and pixel value
            % Package: @ds9
            % Description: Get X,Y position and pixel value in clicked
            %              position.
            % Input  : - Operation mode:
            %            If numeric than will return after the user clicked
            %            the specified number of times.
            %            if 'q' than will return if the user clicked 'q'.
            %            Default is 'q'.
            %          - Operation mode:
            %            'any'   - will return after any character or left click is
            %                      pressed (default).
            %            'key'   - will return after any character is pressed.
            %            'mouse' - will return after mouse left click is pressed.
            % Output : - X
            %          - Y
            %          - Pixel value.
            %          - Cell array of string of clicked events.
            %            '<1>' for mouse left click.
            % Required: XPA - http://hea-www.harvard.edu/saord/xpa/
            % Tested : Matlab 7.0
            %     By : Eran O. Ofek                    Feb 2007
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [X,Y,V]=ds9.getpos(3);
            %          [X,Y,V,Key] = ds9.getpos(1,'key');
            %          [X,Y,V,Key] = ds9.getpos(2,'mouse');
            %          [X,Y,V,Key] = ds9.getpos(1,'any');
            % Reliable: 2
            
            [CooX,CooY,Value,Key]=ds9.ginput('image',varargin{:});
            
        end
            
        function [CooX,CooY,Value,Key]=getcoo(varargin)
            % Interactively get the coordinates (WCS)
            % Package: @ds9
            % Description: Interactively get the coordinates (WCS) and value
            %              of the pixel selected by the mouse (left click) or by clicking
            %              any character on the ds9 display.
            % Input  : - Operation mode:
            %            If numeric than will return after the user clicked
            %            the specified number of times.
            %            if 'q' than will return if the user clicked 'q'.
            %            Default is 'q'.
            %          - Operation mode:
            %            'any'   - will return after any character or left click is
            %                      pressed (default).
            %            'key'   - will return after any character is pressed.
            %            'mouse' - will return after mouse left click is pressed.
            % Output : - RA or Galactic longitude [rad].
            %          - Y/Dec or Galactic latitude [rad].
            %          - Pixel value.
            %          - Cell array of string of clicked events.
            %            '<1>' for mouse left click.
            % Required: XPA - http://hea-www.harvard.edu/saord/xpa/
            % Tested : Matlab 7.0
            %     By : Eran O. Ofek                    Feb 2007
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [X,Y,V]=ds9.getcoo(3);   % return the WCS RA/Dec position
            %          [X,Y,V,Key] = ds9.getcoo('q','key');
            %          [X,Y,V,Key] = ds9.getcoo(2,'mouse');
            %          [X,Y,V,Key] = ds9.getcoo('q','any');
            % Reliable: 2
            
            
            [CooX,CooY,Value,Key]=ds9.ginput('icrs',varargin{:});
            
            
        end
        
        function [MatVal,MatX,MatY]=getbox(Coo,Method,CooType)
            % Get the pixel values in a specified box region
            % Package: @ds9
            %------------------------------------------------------------------------------
            % getbox function                                                          ds9
            % Description: Get from the ds9 display the pixel values in a specified
            %              box region.
            % Input  : - Box ccordinates (in 'image' coordinates):
            %            [Xmin, Xmax, Ymin,Ymax]  (section)
            %            or
            %            [X_left_corner,  Y_left_corner,  X_width,
            %            Y_height] (corner)
            %            or
            %            [X_center,       Y_center,       X_semiwidth,
            %            Y_semiheight] (center).
            %            Alternatively, if empty than the user have to
            %            define two corners of the box using the mouse left
            %            click. Default is empty.
            %          - Method of box position {'section'|'corner'|'center'},
            %            default is 'section'.
            %            Note that in case of Method='center' and CooType='image',
            %            box size will be 2*semiwidth+1.
            %          - Coordinate type {'image'|'physical'|'fk4'|'fk5'|'icrs'},
            %            default is 'image'.
            %            THIS OPTION IS NOT SUPPORTED!
            % Output : - Matrix of pixel values in box.
            %          - Matrix of X position in 'image' coordinates, corresponding
            %            to the matrix of pixel values.
            %          - Matrix of Y position in 'image' coordinates, corresponding
            %            to the matrix of pixel values.
            % Reference: http://hea-www.harvard.edu/RD/ds9/ref/xpa.html
            % Tested : Matlab 7.0
            %     By : Eran O. Ofek                    Feb 2007
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [MatV,MatX,MatY]=ds9.getbox([101 110 101 110]);
            % Reliable: 2
            %------------------------------------------------------------------------------
            Def.Coo      = [];
            Def.Method   = 'section';
            Def.CooType  = 'image';

            if (nargin==0)
                Coo      = Def.Coo;
                Method   = Def.Method;
                CooType  = Def.CooType;
            elseif (nargin==1)
                Method   = Def.Method;
                CooType  = Def.CooType;   
            elseif (nargin==2)
                CooType  = Def.CooType;
            elseif (nargin==3)
                % do nothing
            else
               error('Illegal number of input arguments: ds9.getbox(Coo,Method,CooType)');
            end
            
            if (isempty(Coo))
                % get box coordinates from mouse clicks
                fprintf('Select two corners of the box using the mouse left clicks\n');
                [X,Y]  = ds9.ginput(CooType,2,'mouse');
                Method = 'section';
                Coo    = [min(X), max(X), min(Y), max(Y)];
            end            
            
            Coo = round(Coo);
            
            % convert Coordinate to 'corner' meethod
            switch lower(Method)
                 case 'section'
                     CornerCoo = [Coo(1), Coo(3), Coo(2)-Coo(1)+1, Coo(4)-Coo(3)+1];
                 case 'corner'
                     CornerCoo = Coo;
                 case 'center'
                     CornerCoo = [Coo(1) - Coo(3), Coo(2) - Coo(4), 2.*Coo(3)+1, 2.*Coo(4)+1];
                otherwise
                     error('Unknown Method option');
            end

            SizeMat = [CornerCoo(4), CornerCoo(3)];   % [SizeY, SizeX]
            MatVal  = zeros(SizeMat(1),SizeMat(2));

            switch CooType
                case 'image'

                    [Res] = ds9.system('xpaget ds9 data %s %f %f %f %f no',CooType,round(CornerCoo));

                    ResMat = sscanf(Res,'%d,%d = %f\n',[3 SizeMat(1).*SizeMat(2)]);
                    ResMat = ResMat.';
                    % relative position in matrix (relative to corner)
                    RelX = ResMat(:,1) - min(ResMat(:,1)) + 1;
                    RelY = ResMat(:,2) - min(ResMat(:,2)) + 1;

                    Ind         = sub2ind(SizeMat,RelY,RelX);
                    MatVal(Ind) = ResMat(:,3);
                    %MatVal      = MatVal; %.';

                case 'fk5'
                    error('NOT SUPPORTED');
                    %[~,Res] = ds9_system(sprintf('xpaget ds9 data %s %f %f %f %f no',CooType,round(CornerCoo)));
                    %ResMat = sscanf(Res,'%f,%f = %f\n',[3 SizeMat(1).*SizeMat(2)]);
                    %ResMat = ResMat.';


                otherwise
                    error('CooType option is currently unsupported');
            end

            if (nargout>1)
               [MatX,MatY] = meshgrid((1:1:SizeMat(2)),(1:1:SizeMat(1)));
               MatX        = MatX + CornerCoo(1)-1;
               MatY        = MatY + CornerCoo(2)-1;
            end

        end
        
    end % methods
    
    % Interactive examination
    % (imexam,
    methods (Static)
        
        % plot and return a line profile along the x-axis
        function [Res]=imexam(Image,varargin)
            % ds9 image examination utility
            % Package: @ds9
            % Description: Interactive image examination in ds9.
            % Input  : - Optional image. This may be a SIM object
            %            containing a catalog.
            %            If not given then read the image from the current
            %            ds9 frame.
            %          * Arbitrary number of ...,keyword,value,... pairs.
            %            The following keywords are available:
            %            'Plot' - Display plot. Default is true.
            %            'PlotPar - Cell array of additional arguments to
            %                     pass to the plot command.
            %                     Default is {}.
            %            'AperRad' - Vector of aperture radii in which to
            %                     calculate aperture photometry.
            %                     Default is [2 4 8 12 16].
            %            'Annulus' - Aperture photometry background annulus
            %                     [inner outer] radius.
            %                     Default is [16 22].
            %            'ZP'    - Photometry zero point. Default is 22.
            %            'SemiLen' - Cuts line semi-length [pix].
            %                      Default is 15.
            %            'MeanFun' - Radial plots annular mean function
            %                      handle. Default is @mean.
            %            'Radius' - Radial plots radius [pix].
            %                      Default is 7.
            %            'Sigma'  - PSF sigma guess for first moment
            %                      estimation [pix]. Default is 1.5.
            %            'MaxIter' - Maximum number of centering
            %                      iterations. Default is 3.
            %            'BackPar' - Cell array of additional arguments to
            %                      pass to the SIM/background function.
            %                      Default is {'Block',[128 128]}.
            %            'ExtractionFun' - Source extraction function
            %                      handle. Default is @mextractor.
            %            'ExtractionFunPar' - Cell array of additional
            %                      arguments to pass to the source
            %                      extraction function.
            %            'SearchRad' - Nearest object catalog search
            %                      radius [pix]. Default is 50.
            %            'Field' - SIM field from which to extract cut
            %                      values. Default is 'Im'.
            % Output : - A structure array (element per click) with the
            %            retrieved numerical information.
            
            DefV.Plot             = true;
            DefV.PlotPar          = {};
            DefV.AperRad          = [2 4 8 12 16];
            DefV.Annulus          = [16 22];
            DefV.ZP               = 22;
            DefV.SemiLen          = 15;
            DefV.MeanFun          = @mean;
            DefV.Radius           = 7;
            DefV.Sigma            = 1.5;
            DefV.MaxIter          = 3;
            DefV.BackPar          = {'Block',[128 128]};
            DefV.ExtractionFun    = @mextractor;
            DefV.ExtractionFunPar = {};
            DefV.SearchRad        = 50;
            DefV.Field            = SIM.ImageField;
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);

            if (nargin==0)
                Image = [];
            end
            
            if (isempty(Image))
                % Read entire image from ds9
                Image = ds9.read2sim;
            else
                if (SIM.issim(Image))
                    % image already in SIM format - do nothing
                else
                    Image = images2sim(Image);
                end
            end
            
            Cont = true;
            Ind  = 0;
            while Cont
                Ind = Ind + 1;
                fprintf('\n');
                fprintf('Click h for help, q to abort\n');
                [CooX,CooY,~,Key]=ds9.getpos(1);
                Res(Ind).Type = Key;
                Res(Ind).CooX = CooX;
                Res(Ind).CooY = CooY;
                
                switch Key{1}
                    case {'Right','Left','Up','Down'}
                        % move cursor
                        IndPos = strcmp(Key{1},{'Right','Left','Up','Down'});
                        Pos = [1 0; -1 0; 0 -1; 0 1]; % note that up and down are inverted (ds9)
                        % Move cursor relative position
                        ds9.system('xpaset -p ds9 cursor %d %d',Pos(IndPos,:));
                        
                    case {'h','?','H'}
                        % display help
                        fprintf('\n --- ds9.imexam Menu:\n')
                        fprintf('q   - Abort\n');
                        fprintf('h   - This help\n');
                        fprintf('Left mouse click - recenter on mouse click\n');
                        fprintf('Arrows - move cursor\n');
                        fprintf('v   - Plot vector cut between 2 points\n');
                        fprintf('x,j - Plot vector cut along x axis\n');
                        fprintf('y,i - Plot vector cut along y axis\n');
                        fprintf('s   - Surface plot around point\n');
                        fprintf('c   - Countour plot around point\n');
                        fprintf('m   - First and second moments\n');
                        fprintf('r   - Radial profile with centering\n');
                        fprintf('R   - Radial profile without centering\n');
                        fprintf('a   - Aperture photometry with centering\n');
                        fprintf('A   - Aperture photometry without centering\n');
                        fprintf('p   - PSF photometry with centering\n');
                        fprintf('P   - PSF photometry without centering\n');
                        fprintf('b   - Estimate local background and noise\n');
                        fprintf('S   - Return the nearest source found using mextractor\n');
                        fprintf('e   - Edit/set parameters (and show help for editing)\n');
                    case 'q'
                        % Abort
                        Cont = false;
                    case '<1>'
                        % center on mouse click
                        fprintf('Recenter image on X=%9.3f  Y=%9.3f\n',CooX,CooY);
                        ds9.pan('to',CooX,CooY,'image');
                        
                    case 'v'
                        % plot vector cut between two points
                        fprintf('Click on second point for vector cut plot\n');
                        [CooX2,CooY2,~,~]=ds9.getpos(1);
                        cla;  % clear axis
                        Res(Ind).Res  = vector_prof(Image,[CooX,CooY],[CooX2,CooY2],...
                                                    'Plot',InPar.Plot,'PlotPar',InPar.PlotPar,'Field',InPar.Field);
                    case {'x','j'}
                        % plot x-axis vector around a single point
                        cla;  % clear axis
                        Res(Ind).Res = vector_prof(Image,[CooX-InPar.SemiLen, CooY],...
                                                         [CooX+InPar.SemiLen, CooY],...
                                                   'Plot',InPar.Plot,'PlotPar',InPar.PlotPar,'Field',InPar.Field,...
                                                   'AxisType','x');
                    case {'y','i'}
                        % plot y-axis vector around a single point
                        cla;  % clear axis
                        Res(Ind).Res = vector_prof(Image,[CooX, CooY-InPar.SemiLen],...
                                                         [CooX, CooY+InPar.SemiLen],...
                                                   'Plot',InPar.Plot,'PlotPar',InPar.PlotPar,'Field',InPar.Field,...
                                                   'AxisType','y');
                    case 's'
                        % plot surface
                        cla;  % clear axis
                        Res(Ind).Res = local_surface(Image,[CooX, CooY],InPar.SemiLen,...
                                                     'Plot',InPar.Plot,'PlotPar',InPar.PlotPar,'Field',InPar.Field,...
                                                     'PlotFun',@surface);
                                                 
                    case 'c'
                        % plot contour
                        cla; % clear axis
                        Res(Ind).Res = local_surface(Image,[CooX, CooY],InPar.SemiLen,...
                                                     'Plot',InPar.Plot,'PlotPar',InPar.PlotPar,'Field',InPar.Field,...
                                                     'PlotFun',@contour);
                    case 'm'
                        % first and second moments
                        Res(Ind).Res = moments(Image,[CooX,CooY],...
                                               'Field',InPar.Field,'Radius',InPar.Radius,...
                                               'Sigma',InPar.Sigma,'MaxIter',InPar.MaxIter);
                        Table = astcat_array2table(Res(Ind).Res);
                        fprintf('Clicked position: X=%9.3f   Y=%9.3f\n',Res(Ind).CooX,Res(Ind).CooY);
                        fprintf('1st and 2nd moments (with %d iterations centering):\n',InPar.MaxIter);
                        disp(Table.Cat)
                    case 'r'
                        % radial profile with centering
                        ResM = moments(Image,[CooX,CooY],...
                                               'Field',InPar.Field,'Radius',InPar.Radius,...
                                               'Sigma',InPar.Sigma,'MaxIter',InPar.MaxIter);
                        X = col_get(ResM,'XWIN_IMAGE');
                        Y = col_get(ResM,'YWIN_IMAGE');
                        fprintf('Clicked position: X=%9.3f   Y=%9.3f\n',Res(Ind).CooX,Res(Ind).CooY);
                        fprintf('Centered position after %d iterations: X=%9.3f   Y=%9.3f\n',InPar.MaxIter,X,Y);
                        cla;  % clear axis
                        Res(Ind).Res = rad_prof(Image,[X, Y],...
                                                       InPar.Radius,...
                                                 'Plot',InPar.Plot,'PlotPar',InPar.PlotPar,...
                                                 'MeanFun',InPar.MeanFun,...
                                                 'Field',InPar.Field);
                        Res(Ind).Res.X = X;
                        Res(Ind).Res.Y = Y;
                        
                    case 'R'
                        % radial profile without centering
                        cla;  % clear axis
                        Res(Ind).Res = rad_prof(Image,[CooX, CooY],...
                                                       InPar.Radius,...
                                                 'Plot',InPar.Plot,'PlotPar',InPar.PlotPar,...
                                                 'MeanFun',InPar.MeanFun,...
                                                 'Field',InPar.Field);
                    case 'a'
                        % Aperture photometry with centering
                        
                        CatField = AstCat.CatField;
                        
                        ResM = moments(Image,[CooX,CooY],...
                                               'Field',InPar.Field,'Radius',InPar.Radius,...
                                               'Sigma',InPar.Sigma,'MaxIter',InPar.MaxIter);
                        X = col_get(ResM,'XWIN_IMAGE');
                        Y = col_get(ResM,'YWIN_IMAGE');
                        
                        fprintf('Clicked position: X=%9.3f   Y=%9.3f\n',Res(Ind).CooX,Res(Ind).CooY);
                        fprintf('Centered position after %d iterations: X=%9.3f   Y=%9.3f\n',InPar.MaxIter,X,Y);
                        
                        [Res(Ind).Res,ColAper] = aper_phot(Image,[X, Y],...
                                                 'AperRad',InPar.AperRad,...
                                                 'Annulus',InPar.Annulus);
                                             
                        Naper = numel(InPar.AperRad);
                        fprintf('Aperture photometry\n');
                        for Iaper=1:1:Naper
                            fprintf('  AperPhot Luptitude (rad=%6.2f): %6.3f +/- %6.3f\n',InPar.AperRad(Iaper),...
                                                 convert.flux2mag(Res(Ind).Res.(CatField)(1,ColAper.Aper(Iaper)),InPar.ZP),...
                                                 Res(Ind).Res.(CatField)(1,ColAper.AperErr(Iaper))./...
                                                 Res(Ind).Res.(CatField)(1,ColAper.Aper(Iaper)));
                            fprintf('  AperPhot flux (rad=%6.2f): %6.3f +/- %6.3f\n',InPar.AperRad(Iaper),...
                                                                             Res(Ind).Res.(CatField)(1,ColAper.Aper(Iaper)),...
                                                                             Res(Ind).Res.(CatField)(1,ColAper.AperErr(Iaper)));
                        end
                    case 'A'
                        % Aperture photometry without centering
                        CatField = AstCat.CatField;
                        
                        fprintf('Clicked position: X=%9.3f   Y=%9.3f\n',Res(Ind).CooX,Res(Ind).CooY);

                        
                        [Res(Ind).Res,ColAper] = aper_phot(Image,[CooX, CooY],...
                                                 'AperRad',InPar.AperRad,...
                                                 'Annulus',InPar.Annulus);
                                             
                        Naper = numel(InPar.AperRad);
                        fprintf('Aperture photometry\n');
                        for Iaper=1:1:Naper
                            fprintf('  AperPhot Luptitude (rad=%6.2f): %6.3f +/- %6.3f\n',InPar.AperRad(Iaper),...
                                                 convert.flux2mag(Res(Ind).Res.(CatField)(1,ColAper.Aper(Iaper)),InPar.ZP),...
                                                 Res(Ind).Res.(CatField)(1,ColAper.AperErr(Iaper))./...
                                                 Res(Ind).Res.(CatField)(1,ColAper.Aper(Iaper)));
                            fprintf('  AperPhot flux (rad=%6.2f): %6.3f +/- %6.3f\n',InPar.AperRad(Iaper),...
                                                                             Res(Ind).Res.(CatField)(1,ColAper.Aper(Iaper)),...
                                                                             Res(Ind).Res.(CatField)(1,ColAper.AperErr(Iaper)));
                        end
                                                                             

                    case 'p'
                        % PSF photometry with centering
                    case 'P'
                        % PSF photometry without centering
                        
                        CatField = AstCat.CatField;
                        
                        [Res(Ind).Res,ColPSF] = psf_phot(Image,[CooX, CooY]);
                                                 
                                             
                        fprintf('PSF photometry\n');
                        fprintf('  PSF phot: %6.3f +/- %6.3f\n',Res(Ind).Res.(CatField)(1,ColPSF.Flux),...
                                                                Res(Ind).Res.(CatField)(1,ColPSF.FluxErr));
                    case 'b'
                        % Estimate local background and noise
                        BackField = SIM.BackField;
                        ErrField  = SIM.ErrField;
                        
                        if (~isfield_populated(Image,BackField) || ~isfield_populated(Image,ErrField))
                            % populate the background field
                            Image = background(Image,InPar.BackPar{:});
                        end
                        
                        % get background
                        if (numel(Image.(BackField))==1)
                            Res(Ind).Res.Back = Image.(BackField);
                        else
                            Res(Ind).Res.Back = Image.(BackField)(round(CooY),round(CooX));
                        end
                        % get noise
                        if (numel(Image.(ErrField))==1)
                            Res(Ind).Res.Err = Image.(ErrField);
                        else
                            Res(Ind).Res.Err = Image.(ErrField)(round(CooY),round(CooX));
                        end
                        
                        fprintf('Clicked position: X=%9.3f   Y=%9.3f\n',Res(Ind).CooX,Res(Ind).CooY);
                        fprintf('  Background : %f\n',Res(Ind).Res.Back);
                        fprintf('  Noise      : %f\n',Res(Ind).Res.Err);
                        
                    case 'S'
                        % Return the nearest source found using mextractor
                        % check if catalog exist
                        CatField = AstCat.CatField;
                        if (~isfield_populated(Image,CatField))
                            Image = InPar.ExtractionFun(Image,InPar.ExtractionFunPar{:});
                        end
                        [AstC,Dist] = near_coo(Image,CooX,CooY,InPar.SearchRad,'RadiusUnits','pix');
                        D = col_get(Dist,'Dist');
                        [~,MinI] = min(D);
                        NearCat = row_select(Image,MinI);
                        % show catalog of nearest source
                        fprintf('  Nearest source\n');
                        show(NearCat);
                        
                    case 'e'
                        % Edit/set parameters (and show help for editing)
                        fprintf('Edit the internal ds9.imexam parameters\n');
                        fprintf('Possible parameters:\n');
                        fprintf('  Plot    - plot figure [default is true]\n');
                        fprintf('  AperRad - Aperture phot radii [default is [2 4 8 12 16]]\n');
                        fprintf('  Annulus - Annulus inner/outer radii [default is [16 22]]\n');
                        fprintf('  ZP      - Photometric ZP [default is 22]\n');
                        fprintf('  PlotPar - Plot parameters [default is {}]\n');
                        fprintf('  SemiLen - Line/surface semi length [default is 15]\n');
                        fprintf('  Radius  - Moments and radial radius [default is 7]\n');
                        fprintf('  MeanFun - Function for radial profile mean [default is @mean]\n');
                        fprintf('  Sigma   - Moments sigma of Gaussian weighting [default is 1.5]\n');
                        fprintf('  MaxIter - Number of moments iterations [default is 3]\n');
                        fprintf('  BackPar - Additional pars to pass to SIM/background [default is {}]\n');
                        fprintf('  ExtractionFun- Source extraction prog [default is @mextractor]\n');
                        fprintf('  ExtractionFunPar- Additional pars to pass to extraction prog [default is {}]\n');
                        fprintf('  SearchRad- Sources search radius [default is 50]\n');
                        fprintf('  Field   - SIM field on which to operate [default is Im]\n');
                        
                        Ans = input('Insert parameter and its value (e.g., "SearchRad=100") : ','s');
                        if (~isempty(Ans))
                            eval(sprintf('InPar.%s',Ans));
                        end
                                                                       
                    otherwise
                        fprintf('  Unknown ds9.imexam click option - try again\n');
                end
                        
                        
                        
                        
            
            
            end
        end
            
    end  % end methods
    
    % Interactive ploting
    % (imark, iline, ipoly)
    methods (Static)
        
        % Interactive marker plot
        function [X,Y,Val]=imark(varargin)
            % Interactive plot symbols
            % Package: @ds9
            % Description: Interactive plot symbols in right click
            %              coordinates. Use 'q' click to abort.
            %              Click 'c' to switch to circle plotting, and
            %              's' for box plotting.
            % Input  : * Any of the additional marker type parameters of
            %            the ds9.plot function. I.e., the arguments
            %            following the X,Y coordinates.
            % Output : - Vector of clicked X positions.
            %          - Vector of clicked Y positions.
            %          - Vector of values in clicked positions.
            % Example: ds9.imark;
            %          ds9.imark('Marker','s','Size',[10 10 45]);
            %          ds9.imark('rs');
            %          ds9.imark('rs','size',[20 20 0]);
            %          ds9.imark('ro','Text','Object X');
            % Reliable: 2
        
            Cont = true;
            while Cont
                
                [X,Y,Val,Key]=ds9.getpos(1,'any');
                Marker = [];
                
                switch Key{1}
                    case 'q'
                        % abort
                        Cont = false;
                    otherwise
                        switch lower(Key{1})
                            case '<1>'
                                % left click - mark with previous marker
                                Cont = true;
                            case 'c'
                                % plot circle
                                Marker = 'o';
                            case 's'
                                Marker = 's';
                            otherwise
                                % abort
                                Cont = true;
                        end
                end
                   
                % plot
                if (Cont)
                    if (isempty(Marker))
                        ds9.plot(X,Y,varargin{:});
                    else
                        % use user-defined marker
                        ds9.plot(X,Y,varargin{:},'Marker',Marker);
                    end
                end
                
            end
        end
        
        function [X,Y,ValLine,LineX,LineY]=iline(varargin)
            % Interactively plot a line
            % Package: @ds9
            % Description: Interactively plot a line in ds9 between two
            %              points defined by mouse left clicks.
            %              Also return the interpolated values in the image
            %              along the line.
            % Input  : * Any of the additional marker type parameters of
            %            the ds9.plot function. I.e., the arguments
            %            following the X,Y coordinates.
            % Output : - X coordinates of start and end points.
            %          - Y coordinates of start and end points.
            %          - Interpolated image values along the line.
            %          - X coordinates along the line.
            %          - Y coordinates along the line.
            % Example: [X,Y,ValLine,LineX,LineY]=ds9.iline;
            % Reliabel: 2
            
            InterpMethod = 'linear';
            
            fprintf('Click on two points using mouse left click or any keyboard key\n');
            [X,Y,~,~]=ds9.getpos(2,'any');
            ds9.plot(1,1,varargin{:},'Size',[X(1), Y(1), X(2), Y(2)],'Marker','line');
            if (nargout>2)
                % get values along the line
                
                [MatVal,MatX,MatY]=ds9.getbox([min(X), max(X), min(Y), max(Y)],'section');
                
                Dist       = Util.Geom.plane_dist(X(1),Y(1),X(2),Y(2));
                RoundDist  = round(Dist);
                LineX      = (X(1):(X(2)-X(1))./RoundDist:X(2)).';
                LineY      = (Y(1):(Y(2)-Y(1))./RoundDist:Y(2)).';
                ValLine = interp2(MatX,MatY,MatVal,LineX,LineY,InterpMethod);
            end
                
        end
        
        function [X,Y]=ipoly(varargin)
            % Interactively plot a polygon
            % Package: @ds9
            % Description: Interactively plot a polygon. Verteces are
            %              defined by mouse or keybord clicks.
            %              'q' to finish and abort.
            % Input  : * Any of the additional marker type parameters of
            %            the ds9.plot function. I.e., the arguments
            %            following the X,Y coordinates.
            % Output : - X coordinates of verteces.
            %          - Y coordinates of verteces.
            % Example: [x,y]=ds9.ipoly;
            % Reliable: 2
            
            [X,Y]=ds9.getpos;
            X(end+1) = X(1);
            Y(end+1) = Y(1);
            N = numel(X);
            for I=1:1:N-1
                ds9.plot(1,1,varargin{:},'Size',[X(I), Y(I), X(I+1), Y(I+1)],'Marker','line');
            end
            
        end
                    
    end % methods
    
    % Interactive coordinates to image in browser
    methods (Static)
        function [RA,Dec,Link]=sdssnavi(Browser)
            % Open SDSS navigator for clicked position
            % Package: @ds9
            % Description: Click on a position in an image displayed in ds9 and this
            %              program will open the SDSS navigator web page for the
            %              coordinates.
            % Input  : - Broweser type:
            %            'browser' - existing system browser (default).
            %            'new'     - Matlab internal browser.
            % Output : - J2000.0 RA of position [radians].
            %          - J2000.0 Dec of position [radians].
            %          - Link to SDSS navigator.
            % Required: XPA - http://hea-www.harvard.edu/saord/xpa/
            % Tested : Matlab 7.11
            %     By : Eran O. Ofek                    May 2011
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [RA,Dec,Link]=ds9.sdssnavi;
            % Reliable: 2
            %------------------------------------------------------------------------------
            Def.Browser = 'browser';
            if (nargin==0)
               Browser = Def.Browser;
            elseif (nargin==1)
               % do nothing
            else
               error('Illegal number of input arguments');
            end

            % get RA/Dec
            fprintf('Click on coordinates to get SDSS navigator in browser\n');
            [RA,Dec]=ds9.getcoo(1);

            % get Link to SDSS navigator
            [Link] = VO.SDSS.navigator_link(RA,Dec);
            web(Link{1},sprintf('-%s',Browser));
            Link = Link{1};
        end

           function [RA,Dec,Link]=nedlink(Browser)
            % Open NED link for clicked position
            % Package: @ds9
            % Description: Click on a position in an image displayed in ds9 and this
            %              program will open the NED coordinate search web page for the
            %              coordinates.
            % Input  : - Broweser type:
            %            'browser' - existing system browser (default).
            %            'new'     - Matlab internal browser.
            % Output : - J2000.0 RA of position [radians].
            %          - J2000.0 Dec of position [radians].
            %          - Link to NED search page.
            % Required: XPA - http://hea-www.harvard.edu/saord/xpa/
            % Tested : Matlab 7.11
            %     By : Eran O. Ofek                    May 2011
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: [RA,Dec,Link]=ds9.sdssnavi;
            % Reliable: 2
            %------------------------------------------------------------------------------
            SearchRad = 3;  % arcmin
            Def.Browser = 'browser';
            if (nargin==0)
               Browser = Def.Browser;
            elseif (nargin==1)
               % do nothing
            else
               error('Illegal number of input arguments');
            end

            % get RA/Dec
            fprintf('Click on coordinates to get SDSS navigator in browser\n');
            [RA,Dec]=ds9.getcoo(1);

            % get Link to SDSS navigator
            [Link] = VO.NED.ned_link(RA,Dec,SearchRad);
            web(Link{1},sprintf('-%s',Browser));
            Link = Link{1};
        end

    end % methods
    
    % Interactive coordinate to local/external catalog
    methods (Static)
        function [Cat,RA,Dec]=sdsscat(SearchRadius)
            % Get SDSS catalog near clicked position
            % Package: @ds9
            % Description: Get SDSS catalog near clicked position
            % Input  : - Search radius [arcsec]. Default is 10.
            % Output : - Catalog containing SDSS sources near clicked
            %            position.
            %          - Clicked J2000.0 RA coordinate [rad].
            %          - Clicked J2000.0 Dec coordinate [rad].
            % Example: [RA,Dec,Cat]=ds9.sdsscat;
            % Reliable: 2
            
            RAD = 180./pi;
            
            Def.SearchRadius = 10;
            if (nargin==0)
                SearchRadius = Def.SearchRadius;
            end
            
            % get RA/Dec
            fprintf('Click on position to get SDSS catalog near coordinates\n');
            [RA,Dec]=ds9.getcoo(1);
            
            [Cat]=get_sdss(RA,Dec,SearchRadius./(RAD.*3600));
            
        end
        
        
        
        
    end % methods
    
    % Interact with SIM catalogs
    methods (Static)
        

        function [Cat,MinDist,Units,XRA,YDec,PixVal]=nearestcat(Sim,varargin)
            % Get the nearest source in a SIM/AstCat object
            % Package: @ds9
            % Description: Get the nearest source in a SIM/AstCat object
            %              to the clicked position.
            % Input  : - A SIM or AstCat object.
            %          * Arbitrary number of pairs of ...,key,val,...
            %            arguments. The following keywords are available:
            %            'ColX'  - Cell array of possible keywords
            %                      containing the X coordinate column name
            %                      on which to search.
            %                      The first existing column name will be
            %                      used.
            %                      Default is
            %                      {'XWIN_IMAGE','X','ALPHAWIN_J2000','RA'}
            %            'ColY'  - Like 'ColX' but for the Y coordinate.
            %                      Default is 
            %                      {'YWIN_IMAGE','Y','DELTAWIN_J2000','Dec'}
            %            'ColUnits' - Units corresponding to 'ColX'.
            %                      Default is {'pix','pix','deg','ra'}.
            %            'SrcFun'  - Source extraction function to use.
            %                      Default is @mextractor
            %            'SrcPar'  - Cell array of additional parameters
            %                      to pass to the source extraction
            %                      function. Default is {}.
            %            'SearchRad' - Search radius [rad] or [pix].
            %                      Units defined by availability of column
            %                      in the catalog.
            %                      If empty then return only the nearest
            %                      source.
            %                      Default is empty.
            %            'Verbose' - Default is true.
            % Output : - AstCat object containing the nearest source to
            %            clicked position.
            %          - Distance to nearest source [rad] or [pix].
            %          - Units of RA/Dec (and distance).
            %          - RA or X coordinate of clicked position.
            %          - Dec or Y coordinate of clicked position.
            %          - Value at clicked position.
            % See also: ds9.nearcat
            % Example: [Cat,MinDist,RA,Dec]=ds9.nearestcat(Sim);
            % Reliable: 2
            
            if (~AstCat.isastcat(Sim))
                error('First input argument must be an AstCat object');
            end
            
            RAD = 180./pi;
            CatField     = AstCat.CatField;
            ColCellField = AstCat.ColCellField;
            
            DefV.ColX                = {'XWIN_IMAGE','X','ALPHAWIN_J2000','RA'};
            DefV.ColY                = {'YWIN_IMAGE','Y','DELTAWIN_J2000','Dec'};
            DefV.ColUnits            = {'pix','pix','deg','rad'};
            DefV.SrcFun              = @mextractor;
            DefV.SrcPar              = {};
            DefV.SearchRad           = [];
            DefV.Verbose             = true;
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);
            
            if (~isfield_populated(Sim,ColCellField))
                % Cat is not populated
                % run mextractor
                fprintf('Catalog is not populated - extract sources\n');
                Sim = InPar.SrcFun(Sim,InPar.SrcPar{:});
            end
            [~,ColIndX,UseIndX] = select_exist_colnames(Sim,InPar.ColX(:));
            [~,ColIndY]         = select_exist_colnames(Sim,InPar.ColY(:));
            ColUnits            = InPar.ColUnits{UseIndX};
            
            X = col_get(Sim,ColIndX);
            Y = col_get(Sim,ColIndY);
            switch lower(ColUnits)
                case 'deg'
                    X = X./RAD;
                    Y = Y./RAD;
                    ColUnits = 'rad';
                otherwise
                    % do nothing - already in pix or rad
            end
            
            Cont = true;
            Ind  = 0;
            while Cont
                Ind = Ind + 1;
                switch lower(ColUnits)
                    case {'deg','rad'}
                        % get RA/Dec - spherical coordinates
                        fprintf('Click on position to get SIM catalog near coordinates - q to abort\n');
                        [RA,Dec,Val,Key]=ds9.getcoo(1);

                        D = sphere_dist(X,Y,RA,Dec);
                    case 'pix'
                        % get X/Y - planner coordinates
                        fprintf('Click on position to get SIM catalog near coordinates - q to abort\n');
                        [RA,Dec,Val,Key]=ds9.getpos(1);

                        D = Util.Geom.plane_dist(X,Y,RA,Dec);
                    otherwise
                        error('Unknown ColUnits option');
                end
                
                switch lower(Key{1})
                    case 'q'
                        % abort
                        Cont = false;
                    otherwise
                        [~,MinInd]   = min(D);
                        MinDist(Ind) = D(MinInd);
                        if (isempty(InPar.SearchRad))
                            % found sources
                            IndF     = MinInd;
                        else
                            IndF     = D<InPar.SearchRad;
                        end
                        Cat(Ind)     = row_select(Sim,IndF);
                        Units{Ind}   = ColUnits;
                        XRA(Ind)     = RA;
                        YDec(Ind)    = Dec;
                        PixVal(Ind)  = Val;

                        if (InPar.Verbose),
                            fprintf('\n');
                            fprintf('Clicked position:  %f %f [%s]\n',RA,Dec,Units{Ind})
                            fprintf('Distance: %f [%s]\n',MinDist(Ind),Units{Ind});
                            fprintf('Clicked position value: %f\n',PixVal(Ind));
                            fprintf('-----------------------\n');
                            disp(array2table(Cat(Ind).(CatField),'VariableNames',Cat(Ind).(ColCellField)))
                        end
                end
            end
        end
        
        
        function [Cat,MinDist,Units,XRA,YDec,PixVal]=nearcat(Sim,SearchRad,varargin)
            % Get the nearest source in a SIM/AstCat object
            % Package: @ds9
            % Description: Get the nearest source in a SIM/AstCat object
            %              to the clicked position.
            % Input  : - A SIM or AstCat object.
            %          - Search radius. Default is 100 (in units of queried
            %            column).
            %          * Arbitrary number of pairs of ...,key,val,...
            %            arguments. The following keywords are available:
            %            'ColX'  - Cell array of possible keywords
            %                      containing the X coordinate column name
            %                      on which to search.
            %                      The first existing column name will be
            %                      used.
            %                      Default is
            %                      {'XWIN_IMAGE','X','ALPHAWIN_J2000','RA'}
            %            'ColY'  - Like 'ColX' but for the Y coordinate.
            %                      Default is 
            %                      {'YWIN_IMAGE','Y','DELTAWIN_J2000','Dec'}
            %            'ColUnits' - Units corresponding to 'ColX'.
            %                      Default is {'pix','pix','deg','ra'}.
            %            'SrcFun'  - Source extraction function to use.
            %                      Default is @mextractor
            %            'SrcPar'  - Cell array of additional parameters
            %                      to pass to the source extraction
            %                      function. Default is {}.
            %            'Verbose' - Default is true.
            % Output : - AstCat object containing the nearest source to
            %            clicked position.
            %          - Distance to nearest source [rad] or [pix].
            %          - Units of RA/Dec (and distance).
            %          - RA or X coordinate of clicked position.
            %          - Dec or Y coordinate of clicked position.
            %          - Value at clicked position.
            % See also: ds9.nearcat
            % Example: [Cat,MinDist,RA,Dec]=ds9.nearcat(Sim,100);
            % Reliable: 2
          
            
            RAD = 180./pi;
            CatField     = AstCat.CatField;
            ColCellField = AstCat.ColCellField;
            
            if (nargin==1)
                SearchRad = 100;
            end
            
            DefV.ColX                = {'XWIN_IMAGE','X','ALPHAWIN_J2000','RA'};
            DefV.ColY                = {'YWIN_IMAGE','Y','DELTAWIN_J2000','Dec'};
            DefV.ColUnits            = {'pix','pix','deg','rad'};
            DefV.SrcFun              = @mextractor;
            DefV.SrcPar              = {};
            DefV.Verbose             = true;
            InPar = InArg.populate_keyval(DefV,varargin,mfilename);
            
            [Cat,MinDist,Units,XRA,YDec,PixVal] = ds9.nearestcat(Sim,varargin{:},'SearchRad',SearchRad);

        end
        
        
        
    end % methods
    
    % Interact with SIM images
    methods (Static)
        
        function Res=simval(Sim)
            % Interactively get values from SIM images
            % Package: @ds9
            % Description: Interactively get values from SIM images (image,
            %              background, error, weight and mask) at clicked
            %              positions. Click 'q' to abort.
            % Input  : - A single element SIM object.
            % Output : - A structure array with the results: clicked
            %            positions and pixel values.
            % Example: Res=simval(Sim);
            % Reliable: 
            
            Verbose = true;
            
            ImageField  = SIM.ImageField;
            BackField   = SIM.BackField;
            ErrField    = SIM.ErrField;
            WeightField = SIM.WeightField;
            MaskField   = MASK.MaskField;
            
            Cont = true;
            Ind  = 0;
            while Cont
                % get info from mouse/key click
                fprintf('Click on pixel position - q to abort\n');
                [X,Y,Val,Key]=ds9.getpos(1);
                
                if (strcmpi(Key{1},'q'))
                    % abort
                    Cont = false;
                else
                    % Continue
                    Ind = Ind + 1;
                    
                    Res(Ind).X   = X;
                    Res(Ind).Y   = Y;
                    Res(Ind).Val = Val;
                    Res(Ind).Key = Key{1};
                    FX = floor(X);
                    FY = floor(Y);
                    
                    % get SIM data
                    % Image
                    if (isfield_populated(Sim,ImageField))
                        FI = sub2ind(size(Sim.(ImageField)),FY,FX);
                        Res(Ind).ImageVal = Sim.(ImageField)(FI);
                    else
                        Res(Ind).ImageVal = [];
                    end
                    % background
                    if (isfield_populated(Sim,BackField))
                        if (numel(Sim.(BackField))==1)
                            FI = 1;
                        else
                            FI = sub2ind(size(Sim.(BackField)),FY,FX);
                        end
                        Res(Ind).BackVal = Sim.(BackField)(FI);
                    else
                        Res(Ind).BackVal = [];
                    end
                    % error
                    if (isfield_populated(Sim,ErrField))
                        if (numel(Sim.(ErrField))==1)
                            FI = 1;
                        else
                            FI = sub2ind(size(Sim.(ErrField)),FY,FX);
                        end
                        Res(Ind).ErrVal = Sim.(ErrField)(FI);
                    else
                        Res(Ind).ErrVal = [];
                    end 
                    % weight
                    if (isfield_populated(Sim,WeightField))
                        if (numel(Sim.(WeightField))==1)
                            FI = 1;
                        else
                            FI = sub2ind(size(Sim.(WeightField)),FY,FX);
                        end
                        Res(Ind).WeightVal = Sim.(WeightField)(FI);
                    else
                        Res(Ind).WeightVal = [];
                    end 
                    % mask
                    if (isfield_populated(Sim,MaskField))
                        FI = sub2ind(size(Sim.(MaskField)),FY,FX);
                        Res(Ind).MaskVal = Sim.(MaskField)(FI);
                        % get bit mask name
                        Res(Ind).BitNames = mask2bitname(Sim);
                    else
                        Res(Ind).MaskVal  = [];
                        Res(Ind).BitNames = {};
                    end
                    
                    % show report
                    if (Verbose)
                        fprintf('\n');
                        fprintf('Clicked pixel: %f %f\n',X,Y);
                        fprintf('ds9 image value: %f\n',Val);
                        if (~isempty(Res(Ind).ImageVal))
                            fprintf('SIM image value: %f\n',Res(Ind).ImageVal);
                        end
                        if (~isempty(Res(Ind).BackVal))
                            fprintf('SIM background value: %f\n',Res(Ind).BackVal);
                        end
                        if (~isempty(Res(Ind).ErrVal))
                            fprintf('SIM error value: %f\n',Res(Ind).ErrVal);
                        end
                        if (~isempty(Res(Ind).WeightVal))
                            fprintf('SIM weight value: %f\n',Res(Ind).WeightVal);
                        end
                        if (~isempty(Res(Ind).MaskVal))
                            fprintf('SIM mask value: %f\n',Res(Ind).MaskVal);
                            fprintf('   List of bit mask names\n',Res(Ind).MaskVal);
                            for Im=1:1:numel(Res(Ind).BitNames)
                                fprintf('   %s\n',Res(Ind).BitNames{Im});
                            end
                        end
                        fprintf('\n');
                    end
                    
                end
                
            end
        end
        
    end
    
    
end % end class
            
