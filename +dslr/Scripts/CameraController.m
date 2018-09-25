%Controller for tethered DSLR cameras using <a href=http://digicamcontrol.com/>digiCamControl</a> (Windows app) 
% C = CameraController    -create class and auto-detect digiCamControl
% C = CameraController(ip)  -address of pc running digiCamControl webserver
% C = CameraController(fold)  -folder containing digiCamControl app
% 
%Instruction:
%1.Install and run digiCamControl, BETA v2.0.69 or greater, from:
%  https://sourceforge.net/projects/digicamcontrol/files/latest/download
%2.Enable webserver: File>Settings>Webserver>Enable>RESTART APP.
%3.Connect one or more cameras using USB cable (or WiFi if supported).
%4.For full control set camera to (M) and lens to (MF).
%5.Ensure camera is working through digiCamControl.
%6.Try the examples bellow and read this help.
% 
%Remarks:
%-This class can be used to control supported cameras, stream liveview,
% capture photos & video, download captured files, change settings such as
% ISO, exposure, focus, aperture(fnumber), white balance, compression, etc. 
%-digiCamControl is a multi purpose, free, open source, but Windows only
% application that can control a host of <a href=http://digicamcontrol.com/cameras>supported cameras</a>. 
%-This class communicates with camera(s) via digiCamControl's included 
% <a href=http://digicamcontrol.com/doc/userguide/web>HTTP webserver</a>(recommended) or <a href=http://digicamcontrol.com/doc/userguide/remoteutil>CMD Utility</a>. 
%-The webserver is much faster and allows camera(s) to be controlled from
% any Windows/Linux computer on the network or via the <a href=http://digicamcontrol.com/doc/userguide/settings#webserver>internet</a>. 
%-Visit http://digiCamControl.com for <a href=http://digicamcontrol.com/doc>documentation</a>, <a href=http://digicamcontrol.com/phpbb/>forums</a> and to <a href=http://digicamcontrol.com/donate>donate</a>.
%-Method in this class are Capitalised and have additional descriptions.
%-When this class is created it does a one of retrieval of allowed camera
% options. Redefine this class when swapping cameras.
% 
%Limitations:
%-This class cannot download old photos, user has to use digiCamControl app
% manually: <a href=http://digicamcontrol.com/doc/userguide/interface/downph>digiCamControl>Download photos</a>. 
%-This class can only stream liveview (low-rez, noisy, ~15Hz) from 
% <a href=http://digicamcontrol.com/cameras>supported cameras</a>. However digiCamControl does support "Open Broadcaster
% Software" (OBS) and "XSplit", see <a href=http://digicamcontrol.com/doc/usecases/live>Streaming</a> and <a href=http://digicamcontrol.com/phpbb/search.php?keywords=%5BOBS+%7C+XSplit+%7C+streaming%5D&terms=any&author=&sc=1&sf=all&sr=posts&sk=t&sd=d&st=0&ch=300&t=0&submit=Search>Search Forums</a> for info.
%-This class does not know when capture+download finish, if its > ~3sec.
%-No alphanumeric characters found in some Nikon camera properties are
% being removed. These properties can be queried but can not be set. 
% eg "-", "." in "center-weighted_area" "active_d-lighting" "long_exp._nr"
%-digiCamControl issues: http://digicamcontrol.com/phpbb/viewforum.php?f=4
% 
%Camera Settings:
%-Some settings will not have affect if camera is not in Manual mode (M).
%-To control focus ensure lens is set to Manual Focus (MF):
%-Focus step size & speed can be modified in: <a href=http://digicamcontrol.com/doc/userguide/settings#live-view>File>Settings>Live view</a>
%-Note: Lenses use servo motors which have no discrete physical 'steps'. To
% achieve a specific focus reproducibly try to go to the lens's physical
% limit, in either direction, and apply a set change from there.
% 
%Image Capture:
%-To reduce capture latency from 0.3-0.6 sec to ~0.05s ensure webserver is
% enabled, File>Settings>Webserver>Enable>RESTART APP
%-To measure delay and variance try imaging the computer's own clock by
% calling the "Clock" method provided with this class, C.Clock, however
% I do not know how to measure monitor display latency and variance.
%-Cmd('CaptureAll') will trigger all connected cameras but there will be a
% lag of 0.005-0.020 sec between consecutive cameras.
%-To record video turn on live preview using Cmd('LiveViewWnd_Show') and
% user Cmd('StartRecord') and Cmd('StopRecord').
%-The Capture method blocks code except, but only if acquisition + download
% take more then ~3 sec, then digiCamControl returns without error.
% 
%Image Download:
%-Download is affected by Transfer mode (in app) and session settings.
%-Transfer mode is set via the main app to: PC & Camera | Camera only, if
% set to Camera only some session settings will be ignored. Set it to PC
% & Camera and use session setting "deletefileaftertransfer" if needed.
%-session setting "filenametemplate" works only if "useoriginalfilename" is
% disabled, and it is applied to downloaded files only, not camera files.
% It supports many useful [tags], eg: [Date yyyy-MM-dd], [Time hh-mm-ss],
% [Date yyyy-MM-dd-hh-mm-ss], [Exif.Photo.ExposureTime],
% [Exif.Photo.FNumber], [Exif.Photo.ISOSpeedRatings], etc 
% (for a full list go to: Session>Edit Current Session>File Name Template)
%-"filenametemplate" can be set when calling the Capture method and applies
% to all connected cameras. To distinguish cameras use [Camera Name] or
% [Camera Counter 4 digit] tags in the template.
%-"folder" does not support [tags], instead use "\" in "filenametemplate".
% 
%Ex: download settings, see also "Transfer" in bottom left of GUI
% C = CameraController;
% C.session.folder = 'C:\DSLR';
% C.session.filenametemplate = '[Camera Name]\[Date yyyy-MM-dd-hh-mm-ss]';
% C.session.useoriginalfilename = 0; %ignores "filenametemplate"
% C.session.downloadthumbonly = 0; %not working (v2.0.72.9)
% C.session.downloadonlyjpg = 0; %only used if "PC+CAM"
% C.session.deletefileaftertransfer = 1; %only has affect if Transfer="Cam+PC" and affectively converts it to "PC only"
% C.session.asksavepath = 0; %dialog popup for after capture
% C.session.allowoverwrite = 0; %overwrite if file exists
% C.session.lowercaseextension = 1; %use "*.jpg" instead of "*.JPG"
% 
%Ex: camera settings
% C = CameraController;
% C.camera.isonumber = 100;
% C.camera.fnumber = 4;
% C.camera.shutterspeed = 1/200;
% C.camera.compressionsetting = 'Large Fine JPEG';
% C.camera.drive_mode = 'Single-Frame Shooting';
% 
%Ex: simple capture
% C = CameraController; %initialise
% C.Capture %capture (filename set by "session.filenametemplate") 
% C.Capture('MyPhoto') %capture (set custom filename)
% C.Capture('[Time hh-mm-ss]') %capture (use time tag as filename)
% file = C.lastfile %get last downloaded filenames
% 
%Ex: timed capture
% C = CameraController;
% time = ceil(now*24*60*6)/24/60/60; %upcoming whole second
% file = [datestr(time,'yyyy-mm-dd_HHMMSS.FFF') '_' C.property.devicename]; %timestamp & camera name
% C.Capture(file,time); %capture
% datestr(time)
% 
%Ex: two cameras
% C = CameraController;
% C.session.filenametemplate = '[Camera Name]'; %set filename pattern
% C.Cameras(1), C.property.devicename = 'Cam1'; %change camera name
% C.Cameras(2), C.property.devicename = 'Cam2';
% C.Cmd('CaptureAll')
% 
%Ex: focus stacking
% C = CameraController;
% C.Cmd('LiveViewWnd_Show'), pause(1)  %turn on live preview
% for k = 0:2                          %take 3 photos
%     C.Focus(-2,'small',1)            %two small step towards near focus
%     C.Capture(num2str(k,'Focus%g')); %capture and number the photos
% end
% C.Cmd('LiveViewWnd_Hide')            %turn off live preview to save battery
% 
%Ex: stream live view
%To remove rectangle: Live View>Display>Show focus rectangle
%To reduce lag enable: Live View>Display>No processing
% C = CameraController;
% C.Cmd('LiveViewWnd_Show'); %start live view
% C.Cmd('All_Minimize'); %minimise digiCamControl
% pause(3) %wait for live view
% clf, h = imshow(C.LiveView); %prepare figure
% uicontrol('str','Capture','call','C.Capture') %capture button
% while ishandle(h) %loop until closed
%     set(h,'cdata',C.LiveView) %update live view
%     drawnow %update display
% end
% C.Cmd('LiveViewWnd_Hide'); %stop live view
% 
%Ex: debugging
% C = CameraController;
% C.Clock %show a clock with milliseconds
% C.dbg = 2; %display commands and replies
% C.Capture %capture photo
% 
%Serge 2017
% Questions/bugs/fixes: <a href="http://mailto:s3rg3y@hotmail.com">s3rg3y@hotmail.com</a>
 
%Change Log:
%v1.3.1 (2018-02-24)
%-bunch of minor stuff
%v1.3 (2017-07-22)
%-Support a remote http server
%-Better error handling
%-Allow commas in filenames
%-Minor changes and better help
 
%Update ReadMe.txt:
% fprintf(fopen('ReadMe.txt','w'),'%s',regexprep(help('CameraController'),{'<a .*?>|</a>' '\n ' ' *' '\n\n'},{'' '\n' ' ' '\n \n'}));fclose all
 
%Example webserver commands:
% http://localhost:5513                                                     %primitive http GUI 
% http://localhost:5513/?SLC=CaptureNoAf&param1=Test\[Time%20hh-mm-ss]      %capture and set filename
% http://localhost:5513/?CMD=Capture                                        %capture and display controls webpage with currently selected (previous) photo
% http://localhost:5513/?CMD=CaptureAll                                     %capture with all connected cameras
% http://localhost:5513/?SLC=capture&camera=255076227371                    %capture with specified camera, NOT WORKING, Cam1 fires regardless of number
% http://localhost:5513/?SLC=capture&param1=filename&param2=                %param2 is ???
% http://localhost:5513/preview.jpg                                         %preview (~500k) currently selected photo
% http://localhost:5513/?CMD=LiveViewWnd_Show                               %start and display live preview
% http://localhost:5513/liveview.jpg                                        %live preview current frame
% http://localhost:5513/?CMD=LiveViewWnd_Hide                               %stop live preview
% http://localhost:5513/image/IMG_1200.jpg                                  %download image from hdd (must already have been downloaded from camera)
% http://localhost:5513/thumb/large/IMG_1145.jpg                            %thumb large from hdd (must already have been downloaded from camera)
% http://localhost:5513/thumb/small/IMG_1145.jpg                            %thumb small from hdd (must already have been downloaded from camera)
% http://localhost:5513/session.json                                        %current session data
% http://localhost:5513/?SLC=Get&param1=lastcaptured
% http://localhost:5513/?SLC=List&param1=camera
% http://localhost:5513/?SLC=List&param1=camera.fnumber
% http://localhost:5513/?SLC=Set&param1=session.folder&param2=c:\pictures
% http://localhost:5513/?SLC=Set&param1=session.filenametemplate&param2=capture1
 
%Example CameraControlCmd.exe commands: work when digiCamControl is off, but very SLOW!
% system('"C:\Program Files (x86)\digiCamControl\CameraControlCmd.exe" /filename E:\test\test.jpg /capture')
% system('"C:\Program Files (x86)\digiCamControl\CameraControlCmd.exe" /captureallnoaf')
%See also: http://digicamcontrol.com/doc/userguide/cmd
%Applies to all cameras:
% /help                      - this screen
% /capture                   - capture photo
% /capturenoaf               - capture photo without autofocus
% /captureall                - capture photo with all connected devices
% /captureallnoaf            - capture photo without autofocus with all devices
% /format                    - format camera card(s)
% /session session_name      - use session [session_name]
% /preset preset_name        - use preset [preset_name]
% /folder path               - set the photo save folder
% /filenametemplate template - set the photo save file name template
% /filename fileName         - set the photo save file name
% /counter number            - set the photo initial counter
% /wait [mseconds]           - wait for a keypress or milliseconds
% /nop                       - force past usage with no parameters 
% /verbose                   - lots of status messages 
%Applies to main camera:
% /export filename.txt       - export current connected camera properties 
% /iso isonumber             - set the iso number ex. 100 200 400 
% /aperture aperture         - set the aperture number ex. 9,5 8,0 
% /shutter shutter speed     - set the shutter speed ex. "1/50" "1/250" 1s 3s 
% /ec compensation           - set the exposure comp. -1,5 +2 
% /compression compression   - set the compression Ex: JPEG_(NORMAL) RAW_+_JPEG_(FINE) 
%Nikon only:
% /comment comment           - set in camera comment string 
% /copyright copyright       - set in camera copyright string 
% /artist artist             - set in camera artist string 
 
classdef CameraController < handle
    %% Properties
    properties (Dependent = true) %accessed by other methods
        camera   %GET/SET camera  settings: fnumber, isonumber, shutterspeed, compressionsetting, drive_mode, ...
        session  %GET/SET session settings: folder, filenametemplate, counter, downloadonlyjpg, downloadthumbonly, deletefileaftertransfer, useoriginalfilename, ...
        property %GET/SET device  settings: serialnumber, devicename, nodownload, counter, counterinc, captureinsdram, ...
    end
    properties (SetAccess = private) %read only
        connection  %protocol used to communicate with digiCamControl (set by CheckConnection)
        options     %list of valid camera options
        cmds        %list of valid commands
        lastfile    %last downloaded filename
        lasterr     %last error message
    end
    properties (Hidden = true) %hidden
        dcc %ip/hostname of PC running digiCamControl webserver OR folder on this computer with digiCamControl CMD utility
        dbg %debug level: 1-print requests, 2-print replies (interfears with using TAB for outo-compleate)
    end
    
    %% Constructor
    methods (Hidden = true)
        function C = CameraController(dcc,dbg)
            if nargin<1 || isempty(dcc), C.dcc = ''; else, C.dcc = dcc; end
            if nargin<2 || isempty(dbg), C.dbg = 1;  else, C.dbg = dbg; end
            if C.CheckConnection(C.dcc) %prints error msgs, if any
                disp(['connection type: ' C.connection])
                disp(['digiCamControl:  ' C.dcc])
                [name,serials] = C.CheckCamera;
                if ~isempty(name)
                    t = sprintf('%s, ',serials{:});
                    disp(['camera serials:  ' t(1:end-2)])
                    disp(['current camera:  ' name])
                end
            end
        end
    end
    
    %% Methods
    methods
        function [status,err] = CheckConnection(C,dcc)
            %Check if digiCamControl responds to HTTP or CMD commands.
            %If successful set C.connection to 'HTTP' or 'CMD' else ''.
            % C.CheckConnection    -try to connect to 'localhost' via http
            %                       else dcc install folder via cmd
            % C.CheckConnection(ip)   -use custom webserver address
            % C.CheckConnection(fold)   -use custom cmd utility folder
            % [status,err] = C.CheckConnection
            C.connection = ''; err = ''; %init
            if nargin<2 || isempty(dcc) %try both HTTP and CMD
                C.dcc = 'localhost';
                [status,err1] = C.TestHTTP(C.dcc);
                if status
                    C.connection = 'HTTP'; %success
                else
                    [C.dcc,err2] = C.FindDCC;
                    if ~isempty(err2)
                        err = [err1 ' + ' err2];
                    else
                        [status,err3] = C.TestCMD(C.dcc);
                        if status
                            C.connection = 'CMD'; %success (slow connection)
                            C.Error('Turn on HTTP webserver to reduce latency and to stream LiveView')
                        else
                            err = [err1 ' + ' err3];
                        end
                    end
                end
            else %try one
                switch upper(dcc)
                    case 'HTTP', C.dcc = 'localhost'; %default HTTP
                    case 'CMD',  [C.dcc,err] = C.FindDCC;
                    otherwise,   C.dcc = dcc;
                end
                if ~isempty(C.dcc)
                    if isdir(C.dcc) %#ok<ISDIR> %assume cmd
                        [status,err] = C.TestCMD(C.dcc);
                        if status
                            C.connection = 'CMD';
                        end
                    else
                        [status,err] = C.TestHTTP(C.dcc);
                        if status
                            C.connection = 'HTTP';
                        end
                    end
                end
            end
            C.Error(err,nargout<2)
        end
        function [name,serials,err] = CheckCamera(C)
            %Detect camera name (if any) and update its allowed setting
            name = ''; %init
            [serials,err] = C.Cameras;
            if ~isempty(serials) %is a camera connected
                [name,err] = C.Get('property.devicename');
                C.options = C.Options; %cache camera options
                C.cmds = C.List('Cmds'); %cache commands
            end
            C.Error(err,nargout<3)
        end
        
        function [out,err] = Capture(C,file,time,mode,lag,block)
            %Capture photo, now or at set timed, with one or all cameras
            % Capture         -capture photo now
            % Capture(file)      -filename, if downloaded (no extension)
            % Capture(file,time)    -start time (datenum) or -delay (sec)
            % Capture(file,time,mode)   -{'noaf'} 'af' 'all'
            % Capture(file,time,mode,lag)  -capture lag (sec) {0.05 or 0.4}
            % Capture(file,time,mode,lag,blck) -wait until finished {true}
            % [out,err] = Capture(..)          -return error messages
            %file: filename, if downloaded to pc (no extension)
            %      eg ['[Date yyyy-MM-dd]\' datestr(now,'HH-MM-SS.FFF')]
            %time: if empty - capture immediately
            %      if positive or string - absolute time of when to capture
            %         (ie datenum or datestr, eg '2018-05-25 23:38:44.50')
            %      if negative - a delay for capture in seconds (eg -10)
            %mode: 'noaf' - capture without autofocus (default)
            %      'af'   - autofocus and if successful then capture
            %      'all'  - capture with all connected cameras
            %lag:  start capture this many seconds ahead of specified time
            %      if time is absolute, default: 0.05 (HTTP) or 0.4 (CMD)
            %block:if true wait for capture to be completed (default),
            %      but HTTP may time out and return while capture and
            %      download is still occurring.
            if nargin<2 || isempty(file), file = '';   end %custom file name
            if nargin<3 || isempty(time), time = [];   end %start capture now or at this absolute time
            if nargin<4 || isempty(mode), mode = 'no'; end %capture mode, 'NoAF' 'AF' 'All'
            if nargin<5 || isempty(lag),  lag  = [];   end %for timed capture start capture this many seconds ahead of specified time
            if nargin<6 || isempty(block),block = true;end %wait for capture to be completed before returning
            if ~isempty(file) && (strcmpi(mode(1:2),'al') || any(file==' ') && isequal(C.connection,'CMD'))
                %If user wants a custom filename then set it now if:
                %1) CaptureAll is being used
                %2) filename has a space & cmd utility is being used (SLOW! CMD commands take around 0.4s to run!)
                C.Run(['Set session.filenametemplate ' file]) %set the file
                file = '';
            end
            if ~isempty(time) %is this a timed capture
                if isscalar(time) && time<=0 %delay of minus this many seconds
                    pause(abs(time))
                else %time must be datestr, datevec or datenum
                    if isempty(lag)
                        switch C.connection
                            case 'HTTP', lag = 0.05;
                            case 'CMD',  lag = 0.40;
                        end
                    end
                    time = datenum(time)-lag/24/60/60; %time = datenum(time)-(fudg+focus*0.03-~wait*0.04)/24/60/60;
                    while time>now %loop allows clock adjustments during wait
                        pause(min((time-now)*24*60*60,0.5))
                    end
                end
            end
            switch lower(mode(1:2))
                case 'no', [out,err] = C.Run('CaptureNoAf',file,[],block);
                case 'af', [out,err] = C.Run('Capture'    ,file,[],block);
                case 'al', [out,err] = C.Cmd('CaptureAll'); %note CaptureAll does not allow extra argument
            end
            if ~nargout
                clear out
            end
        end
        
        function [I,err] = LiveView(C)
            err = ''; I = []; %init
            if strcmp(C.connection,'HTTP')
                try
                    I = imread(['http://' C.dcc ':5513/liveview.jpg'],'jpg');
                catch e
                    if strcmp(e.identifier,'MATLAB:imagesci:imread:readURL')
                        err = 'HTTP connection timed out';
                    elseif strcmp(e.identifier,'MATLAB:imagesci:jpeg_depth:unhandledLibraryError')
                        err = 'Live view not active';
                        %alternative is to turn on live view now and wait
                        %upto ~4 sec for image, however when live view is
                        %turned off the old live view image is cached and
                        %there is no way to tell that live view is actualy
                        %off. So onus is on user to ensure live view is
                        %turned on and is running.
                        % warning('off','MATLAB:imagesci:jpeg_depth:libraryMessage') %supress repeated warnings
                        % C.Cmd('LiveViewWnd_Show'); %start live view
                        % C.Cmd('All_Minimize');     %minimise live view window
                        % for k = 1:40
                        %     pause(0.1)
                        %     try
                        %         I = imread(['http://' C.dcc ':5513/liveview.jpg'],'jpg');
                        %     end
                        %     if ~isempty(I)
                        %         break
                        %     end
                        % end
                        % if isempty(I)
                        %     err = 'Live view not active';
                        % end
                    else
                        err = e.message;
                    end
                end
            elseif strcmp(C.connection,'CMD')
                err = 'LiveView only works with HTTP webserver';
            else
                err = 'Check connection';
            end
            C.Error(err,nargout<2)
        end
        
        function [out,err] = Cameras(C,val)
            %List connected cameras or select a camera
            % SN = Cameras        -list of connected cameras serials
            % SN = Cameras(index)  -set current camera using list index
            % SN = Cameras(serial)  -set current camera using serial number
            % [SN,err] = Cameras(.)  -return error string
            %Use property.serialnumber to get current camera's serial
            [serials,err] = C.List('cameras'); %all cameras serial numbers
            if isempty(serials) || strcmpi(serials,'OK')
                err = 'No camera(s) detected';
                out = '';
            elseif nargin>1 %select a specific camera
                if isnumeric(val) %index selection mode
                    try
                        val = serials{val};
                    catch %#ok<CTCH>
                        out = '';
                        err = 'Invalid camera index';
                        C.Error(err,nargout<2)
                        return
                    end
                end
                old = C.Get('property.serialnumber'); %current camera serial number
                [out,err] = C.Set('camera',val,old,serials); %change camera if different
                if ~strcmp(out,old) %has selected camera changed
                    C.options = C.Options; %cache new camera options
                    C.cmds = C.List('Cmds'); %cache single-line-commands
                end
            else
                out = serials;
            end
            C.Error(err,nargout<2)
        end
        
        function [out,err] = Sessions(C,val)
            %List available session or set current session
            % names = Sessions    -list available session names
            % name = Sessions(name)  -set current session
            % [name,err] = Sessions(.)  -return error string
            %Note: use session.name to get current session
            [list,err] = C.List('sessions'); %all session names
            if nargin>1
                if isnumeric(val) %index selection mode
                    try
                        val = list{val};
                    catch %#ok<CTCH>
                        out = '';
                        err = 'Invalid camera index';
                        C.Error(err,nargout<2)
                        return
                    end
                end
                old = C.Get('session'); %current session name
                [out,err] = C.Set('session',val,old,list);
            else
                out = list;
            end
        end
        
        function [status,err] = Cmd(C,cmd)
            %Run a single line command, or list available commands.
            % cmds = Cmd         -list of single-line-commands (cellstr)
            % [cmd,err] = Cmd(cmd) -run command and return errors (string) 
            status = 0; %init
            [I,mch] = C.Match(cmd,C.cmds);
            if sum(I) == 1
                [~,err] = C.Run('Do',mch); %Do commands do not return anything
                if isempty(err)
                    status = 1;
                end
            else
                err = sprintf('Invalid command: options are:%s',sprintf('\n''%s''',C.cmds{:}));
            end
            C.Error(err,nargout<2)
        end
        
        function [status,err] = Focus(C,Num,Mode,Wait)
            %Adjust camera focus, or auto-focus
            % Focus([])   -auto focus, lens must be set to AF
            % Focus(Num)     -number of steps, +ve=far field -ve=near field
            % Focus(Num,Mode)   -type of step {'small'} 'med' 'large'
            % Focus(Num,Mode,Wait)  -time delay per step (sec)
            %Starts live view, lens can be in MF|AF, camera can be in M|A..
            %Step size can be set in: File>Settings>Live view
            C.Cmd('LiveViewWnd_Show') %can skip if LiveView is on
            if nargin<2 || isempty(Num)
                [status,err] = C.Cmd('LiveView_Focus'); %auto focus, user must wait for focus to finish manually
            elseif Num ~= floor(Num)
                err = 'Focus step must be an integer';
            elseif Num==0
                %do nothing
            else
                if nargin<3 || isempty(Mode) %default step mode
                    Mode = 'small';
                end
                if nargin<4 || isempty(Wait) %default delay
                    switch lower(Mode(1))
                        case {1 's'}, Wait = 1; %adjust these defaults
                        case {2 'm'}, Wait = 5;
                        case {3 'l'}, Wait = 15;
                    end
                end
                if Num>0, cmd = 'P'; %move towards far focus
                else,     cmd = 'M'; %move towards near focus
                end
                switch lower(Mode(1))
                    case {1 's'} %do nothing
                    case {2 'm'}, cmd = [cmd cmd];
                    case {3 'l'}, cmd = [cmd cmd cmd];
                end
                for k = 1:abs(Num)
                    [status,err] = C.Cmd(['LiveView_Focus_' cmd]); %send command
                    pause(Wait) %wait for focus adjustment to hopefully finish
                end
            end
            C.Error(err,nargin<2)
        end
 
        function Clock(~,run_in_this_session)
            %Display a clock with miliseconds for time tests
            % Clock   -start in new MatLab session so as not to block code
            % Clock(1)  -start using this session, blocking code execution
            if nargin<2 || ~run_in_this_session
                fprintf('Starting another MatLab session to display clock...\n')
                %UGLY! can this be done using threads?
                %!matlab -nodesktop -nosplash -minimize -r "try C=CameraController;C.Clock(1);catch,exit,end" & 
                !matlab -nosplash -minimize -r "try C=CameraController;com.mathworks.mde.desk.MLDesktop.getInstance.getMainFrame.hide;C.Clock(1);catch,exit,end" &
            else
                clf(figure(1)), axis off
                set(gcf,'color','k','name','Clock','numb','off','menu','n','tool','n')
                h0 = text(0.5,0.95,'time'    ,'fontsize',60,'hor','cen','color','w');
                h1 = text(0.5,0.7 ,'seconds' ,'fontsize',60,'hor','cen','color','w');
                h2 = text(0.5,0.5 ,'tenths'  ,'fontsize',60,'hor','cen','color','w');
                h3 = text(0.5,0.2 ,'hundreds','fontsize',60,'hor','cen','color','w');
                h4 = text(0.5,-0.05,'Delay:' ,'fontsize',30,'hor','cen','color','w');
                old = now; %time of previous frame
                num = 10; %measure time elapsed between n frames to display average
                lag = nan(1,num);
                idx = 1; %circular counter
                %addlistener(h4,'String','PreSet',@(~,~)H.UpdateClock(h0,h1,h2,h3,h4,old,num,OLD,idx)); %this isn't working execution
                while true
                    pause(0.0001)
                    new = now;
                    p1 = floor(mod(new*24*60*60    ,10)); %sec
                    p2 = floor(mod(new*24*60*60*10 ,10)); %1/10 sec
                    p3 = floor(mod(new*24*60*60*100,10)); %1/100 sec
                    set(h0,'str',datestr(new,'HH:MM:SS.FFF'))
                    set(h1,'pos',[p1/10 0.6 0],'str',num2str(p1))
                    set(h2,'pos',[p2/10 0.4 0],'str',num2str(p2))
                    set(h3,'pos',[p3/10 0.2 0],'str',num2str(p3))
                    lag(idx) = new-old;   %update LAG amounts
                    idx = mod(idx,num)+1; %update circular counter
                    old = new;            %update oprevious frame time
                    set(h4,'str',num2str(nanmean(lag)*24*60*60,'Delay: >%.3f sec'))
                    drawnow
                end
            end
        end
        
    end
    
%     methods
%         function a=UpdateClock(H,h0,h1,h2,h3,h4,o,n,D,c)
% %             try %when figure closes, fail silently
%                 while true
%                     pause(0.0001)
%                     t = now;
%                     p1 = floor(mod(t*24*60*60    ,10)); %sec
%                     p2 = floor(mod(t*24*60*60*10 ,10)); %1/10 sec
%                     p3 = floor(mod(t*24*60*60*100,10)); %1/100 sec
%                     set(h1,'pos',[p1/10 0.7 0],'str',num2str(p1))
%                     set(h2,'pos',[p2/10 0.5 0],'str',num2str(p2))
%                     set(h3,'pos',[p3/10 0.3 0],'str',num2str(p3))
%                     d = t - o;
%                     D(c) = d;
%                     c = mod(c,n)+1;
%                     o = t;
%                     t = datestr(t,'HH:MM:SS.FFF');
%                     set(h0,'str',t)
%                     set(h4,'str',num2str(nanmean(D)*24*60*60,'%.3fs'))
%                     drawnow
%                 end
% % H.UpdateClock(h0,h1,h2,h3,h4,o,n,D,c)
% %             end
%         end
%     end
    
    %% Hidden methods
    methods (Hidden = true)
        function [out,err] = Run(C,cmd,prp,val,block)
            %Send HTTP or CMD command to digiCamControll
            % Run(cmd)       -command (string)
            % Run(cmd,prp)       -argument or property {''} (string)
            % Run(cmd,prp,val)       -property value {''}
            % Run(cmd,prp,val,block) -wait for command to finish {true}
            % [out,err] = Run(..)    -output reply and error strings
            if nargin<3 || isempty(prp),  prp  = '';  end
            if nargin<4 || isempty(val),  val  = '';  end
            if nargin<5 || isempty(block),block= true;end
            out = ''; err = ''; %init
            switch C.connection
                case 'HTTP' %webserver
                    cmd = ['http://' C.dcc ':5513/?SLC=' cmd]; 
                    if ~isempty(prp)
                        cmd = [cmd '&param1=' regexprep(prp,{' ' '=' ';'},{'%20' '%3D' '%3B'})]; %allow spaces and equal signs [www.w3schools.com/tags/ref_urlencode.asp]
                    end
                    if ~isempty(val)
                        val = regexprep(val,{' ' '='},{'%20' '%3D'}); %replace spaces and equal signs
                        val = regexprep(val,'[<>"|?*]',''); %windows filenames forbid [<>:"/\|?*], [/] is needed for fractions, [:\] is used for subfolders
                        if ~isempty(val)
                            cmd = [cmd '&param2=' val];
                        end
                    end
                    if C.dbg >= 2 %display HTTP command
                        %disp(['>> ' cmd]) %plain text
                        %disp(['>> urlread(' cmd ')']) %matlab command
                        %disp(['<a href="' cmd '">' cmd '</a> ']) %clickable link that opens browser 
                        disp(['<a href="matlab:urlread(''' cmd ''')">' cmd '</a> ']) %clickable link that runs in matlab
                    end
                case 'CMD'
                    cmd = ['"' fullfile(C.dcc,'CameraControlRemoteCmd.exe') '" /c ' cmd ]; %spaces are not allowed in filename, this also prohibits the use of most tags
                    if ~isempty(prp)
                        cmd = [cmd ' "' strrep(prp,' ','_') '"'];
                    end
                    if ~isempty(val)
                        val = regexprep(val,'[<>"|?*;]',''); %windows filename forbids [<>:"/\\|?*], dcc cmd forbids [;], [/] is needed for fractions, [:\] is used for subfolders
                        if ~isempty(val)
                            cmd = [cmd ' ' val];
                        end
                    end
                    if C.dbg >= 2 %display CMD command
                        %disp(['>> ' cmd]) %plain text
                        disp(['>> system(''' cmd ''')']) %matlab command
                        % disp(['<a href="' cmd '">' cmd '</a> ']) %clickable link that opens browser 
                        % disp(['<a href="matlab:system(''"' cmd ''')">' cmd '</a> ']) %clickable link that runs in matlab (broken) 
                    end
            end
            
            if strcmp(C.connection,'HTTP')
                if block %wait while command executes
                    [out,status] = urlread(cmd); %send httm request and read reply
                        if C.dbg >= 3
                            disp(out) %display replies
                        end
                        if ~isempty(out) && out(end)==10 %remove trailing linefeeds
                            out(end) = [];
                        end
                        if isempty(out)
                            %should this issue an error???
                        end
                        out = strrep(out,'Cannot perform runtime binding on a null reference','No camera detected'); %improve some error msgs
                        out = strrep(out,'Unknow ','Unknown ');
                        if ~status || strcmp(out,'Unknown parameter') || strncmpi(out,'Wrong value',11) || strcmp(out,'No camera detected')
                            err = out;
                            out = '';
                        end
                else %dont wait
                    [~,~] = system('start /B curl http://localhost:5513/?SLC=CaptureNoAf'); %ignore output
                end
            else %CMD
                if block
                    [failed,out] = system(cmd); %run cmd command and read the reply
                    if C.dbg >= 3
                        disp(['ans =' 10 out]) %display replies
                    end
                    [out,err] = C.CleanCMD(out);
                    if C.dbg >= 4
                        disp(['ans =' 10 out]) %display replies
                    end
                    if failed || strncmpi(out,'error',5)
                        err = out;
                        out = '';
                    end
                else
                    try
                        out = java.lang.Runtime.getRuntime().exec(cmd);
                    catch e
                        err = sprintf(e.message);
                    end
                end
            end
            C.Error(err,nargout<2)
        end
        function [out,err] = Get(C,prop)
            [out,err] = C.Run('Get',prop);
            out = regexprep(out,{'True' 'False'},{'true' 'false'}); %HACK
            C.Error(err,nargout<2)
        end
        function [new,err] = Set(C,prp,val,old,opt)
            % [new,err] = Set(C,prp,val,old,opt)
            if nargin<4, old = ''; end %current value, to avoid superfluos set commands
            if nargin<5, opt = {}; end %default valid options
            new = ''; err = ''; %init
            if ~ischar(val)
                val = mat2str(val); %convert numbers to strings
            end
            if ~isempty(opt) %are valid options known
                [~,val] = C.Match(val,opt); %find match
                if isempty(val) %no match found
                    err = sprintf('Invalid assignment: options are:%s',sprintf(' ''%s''',opt{:}));
                    C.Error(err,nargout<2)
                    return
                end
            end
            if strcmpi(old,val) %is old value known and is it different
                new = val; %do nothing
            else
                [~,err] = C.Run('Set',prp,val); %send command
                new = C.Get(prp); %varify value after set (can skip this)
                if ~isequal(new,val) && ~isequal(str2num(lower(new)),str2num(lower(val))) %#ok<ST2NM> %verify success, allows: 'True'='true'=true=1
                    err = sprintf('Set command failed: %s',err);
                end
            end
            C.Error(err,nargout<2)
        end
        function [out,err] = List(C,cmd)
            [out,err] = C.Run('List',cmd);
            if ~isempty(out) && isempty(err)
                out = regexprep(out,{'True' 'False'},{'true' 'false'}); %HACK
                out = regexp(strtrim(out),'\n','split')';
            end
            C.Error(err,nargout<2)
        end
        function [s,err] = Options(C)
            %Get a list of valid camera options as struct of cellstr
            params = fieldnames(C.camera); %list of parameters
            for k = 1:numel(params)
                [s.(params{k}),err] = C.List(['camera.' params{k}]); %list options for each parameter
                s.(params{k}) = s.(params{k})';
            end
        end
    end
    
    %% get/set methods
    %Set methods do not know which sub-field(s) were set, to avoid setting
    %all fields they GET current values and SET only those that changed.
    methods
        function s = get.camera(C) %Get current camera settings as struct, empty if no camera
            s = C.List('camera'); %camera settings as cellstr, eg 'camera.fnumber=4.0'
            if ~isempty(s)
                s = regexp(s,'camera\.(.*?)=(.*)','tokens','once'); %split fields and values, eg {{'fnumber' '4.0'};...}
                s = cat(1,s{:})'; %form a cellstr table 2-by-n
                s(1,:) = regexprep(s(1,:),'[^\w]',''); %remove "." "-" from field names (Nikon), set methods will not work for affected fields
                s = struct(s{:}); %make a struct
                if isfield(s,'exposurestatus')
                    s = rmfield(s,'exposurestatus'); %"exposurestatus" is read only and does not appear to change with a Canon
                end
            end
        end
        function s = get.session(C)
            s = C.List('session');
            if ~isempty(s)
                s = regexp(s,'session\.(.*?)=(.*)','tokens','once');
                s = cat(1,s{:})';
                %s(2,strcmpi(s(2,:),'True' )) = {true};
                %s(2,strcmpi(s(2,:),'False')) = {false};
                s(2,:) = regexprep(s(2,:),{'^False$' '^True$'},{'false' 'true'}); %change case for consistancy
                s = struct(s{:});
            end
        end
        function s = get.property(C)
            s = C.List('property');
            if ~isempty(s)
                s = regexp(s,'property\.(.*?)=(.*)','tokens','once');
                s = cat(1,s{:})';
                %s(2,strcmpi(s(2,:),'True' )) = {true};
                %s(2,strcmpi(s(2,:),'False')) = {false};
                s(2,:) = regexprep(s(2,:),{'^False$' '^True$'},{'false' 'true'}); %change case for consistency
                s = struct(s{:});
            end
        end
        function c = get.lastfile(C)
            c = C.Get('lastcaptured');
            if any(strcmp(c,{'-' '?'}))
                c = '';
            end
        end
        %Due to the way this class is structured attempts to set a property
        %will first trigger a get method (not needed) and then set method.
        function set.camera(C,new)
            %dcc can list valid option for each camera.parameter, so we
            %can do a bit of common sense check to guess which value the
            %user wanted or tell the user valid options if we can't figure
            %it out. But there are a few strange transient parameter that
            %make this tricky, eg camera.focusmode & camera.bracketing
            %which are sometimes not reported by "Get camera".
            old = C.camera; %get current settings
            if ~isempty(old)
                for f = fieldnames(new)' %step through parameters
                    if ~isfield(C.options,f{1}) %have options for this parameter been cached
                        [list,err] = C.List(['camera.' f{1}]); %is this a valid parameter
                        if isempty(err)
                            C.options.(f{1}) = list'; %update valid parameter options
                        else
                            C.Error(err) %display error
                            continue
                        end
                    end
                    if ~isfield(old,f{1})
                        old.(f{1}) = ''; %init
                    end
                    C.Set(['camera.' f{1}],new.(f{1}),old.(f{1}),C.options.(f{1})); %Set each property if it changed and if value is valid
                end
            end
        end
        function set.session(C,new)
            old = C.session;
            for f = fieldnames(new)'
                C.Set(['session.',f{1}],new.(f{1}),old.(f{1}),'');
            end
        end
        function set.property(C,new)
            old = C.property;
            for f = fieldnames(new)'
                C.Set(['property.',f{1}],new.(f{1}),old.(f{1}),'');
            end
        end
    end
    
    %% Helper functions
    methods (Static = true, Hidden = true)
        function [status,err] = TestHTTP(ip)
            %Test HTTP webserver communication
            % [status,err] = TestHTTP(ip)
            status = 0; err = ''; %init
            try
                t = java.net.URL([],['http://' ip ':5513'],sun.net.www.protocol.http.Handler).openConnection; %does this work on linux?
                %t.setConnectTimeout(0.5); t.setReadTimeout(0.5); %timeout not working, defaults to ~2.5 seconds
                t.getInputStream;
                status = 1; %success
            catch e
                if     strfind(e.message,'connect timed out'), err = 'HTTP connection timed out';
                elseif strfind(e.message,'Connection refused'),err = 'HTTP connection refused';
                elseif strfind(e.message,'UnknownHost'),       err = 'HTTP unknown address';
                elseif strfind(e.message,'Permission denied'), err = 'HTTP permission denied';
                else,                                          err = ['HTTP ' e.message];
                end
            end
            %errors not displayed
        end
        function [fold,err] = FindDCC
            %Check default install location for digiCamControl app
            % [fold,err] = FindDCC
            err = ''; %init
            if     isdir(fullfile(char(java.lang.System.getenv('ProgramFiles(x86)')),'digiCamControl')) %#ok<ISDIR>
                fold   = fullfile(char(java.lang.System.getenv('ProgramFiles(x86)')),'digiCamControl');
            elseif isdir(fullfile(char(java.lang.System.getenv('ProgramFiles'     )),'digiCamControl')) %#ok<ISDIR>
                fold   = fullfile(char(java.lang.System.getenv('ProgramFiles'     )),'digiCamControl');
            else
                fold = '';
                err = 'digiCamControl folder not found';
            end
            %errors not displayed
        end
        function [status,err] = TestCMD(dccfolder)
            %Test CMD utility communication
            % [status,err] = TestCMD(dccfolder)
            status = 0; err = ''; %init
            [~,t] = system('tasklist /FI "imagename eq CameraControl.exe"'); %can skip this test, but will take a long time to fail if app is not running
            if ~contains(t,'CameraControl.exe')
                err = 'digiCamControl is not running';
            else
                exe = fullfile(dccfolder,'CameraControlRemoteCmd.exe');
                if ~exist(exe,'file')
                    err = 'CameraControlRemoteCmd.exe not found';
                else
                    [fail,msg] = system(['"' exe '" /c Get "property.devicename"']);
                    if ~fail
                        status = 1; %success
                    else
                        err = ['CMD ' msg];
                    end
                end
            end
            %errors not displayed
        end
        function [str,err] = CleanCMD(str)
            %clean up cmd responce (faster then using /clean argument)
            str = regexp(str,'(?<=:;response:).*?(?=[;]*\n)','match','once');
            if strncmpi(str,'error',5) %was an error generated
                err = regexprep(str,'error;message:',''); %clean error message
                str = '';
            else
                err = '';
                if ~isempty(str) && strcmp(str(1),'[') %is string a list
                    str = regexprep(str,{'","' '^[' ']$'},{'"\n"' '' ''}); %replace commas with linefeed, remove start end braces
                end
                str = regexprep(str,{'\\\\' '"'},{'\\' ''}); %remove any escape character and quotes
            end
            %errors not displayed
        end
        function [I,mch] = Match(str,opt)
            %Flexible comparison of cellstr to string, return best match
            % [I,mch] = Match(str,opt)
            I = strcmpi(str,opt); %compare whole string
            if ~any(I) %compare start of string
                num = str2num(str); %#ok<ST2NM> eg '1/200'=0.005
                if ~isempty(num)
                    I = cellfun(@(x)isequal(str2num(x),num),opt); %#ok<ST2NM> may want to do a rounding tolerance
                end
                if ~any(I) %try numeric comparison
                    I = strncmpi(opt,str,numel(str));
                    if ~any(I) %try fragment of string comparison
                        I = contains(lower(opt),lower(str)); 
                    end
                end
            end
            if sum(I)==1
                mch = opt{I};
            else
                mch = '';
            end
        end
    end
    methods (Access = private)
        function Error(C,err,display)
            if ~isempty(err) && (nargin<3 || display)
                C.lasterr = err; %save last error
                fprintf(2,'%s\n',err); %display error using red text, but do not abbort execution
            end
        end
    end
    
    %% Hide handle methods by overloading them
    methods (Hidden = true)
        function l = addlistener(varargin), l = addlistener@handle(varargin{:}); end
        function     notify     (varargin),     notify@handle     (varargin{:}); end
        function     delete     (varargin),     delete@handle     (varargin{:}); end
        function h = findobj    (varargin), h = findobj@handle    (varargin{:}); end
        function p = findprop   (varargin), p = findprop@handle   (varargin{:}); end
        function b = eq(varargin),          b = eq@handle(varargin{:});          end
        function b = ne(varargin),          b = ne@handle(varargin{:});          end
        function b = lt(varargin),          b = lt@handle(varargin{:});          end
        function b = le(varargin),          b = le@handle(varargin{:});          end
        function b = gt(varargin),          b = gt@handle(varargin{:});          end
        function b = ge(varargin),          b = ge@handle(varargin{:});          end
        %isvalid is a sealed method and cannot be overloaded
    end
end