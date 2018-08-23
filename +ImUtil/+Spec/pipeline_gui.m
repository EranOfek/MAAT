function varargout = pipeline_gui(varargin)
% Spectrocopy pipeline GUI
%Input: AllData, AllPeaks
%Output: AllData, AllPeaks


% Palomar P200 spectrocopy pipeline GUI
% Author: Adam Rubin, Weizmann Institute of Science
% Last revision: 2014-03-02



% PIPELINE_GUI MATLAB code for pipeline_gui.fig
%      PIPELINE_GUI, by itself, creates a new PIPELINE_GUI or raises the existing
%      singleton*.
%
%      H = PIPELINE_GUI returns the handle to a new PIPELINE_GUI or the handle to
%      the existing singleton*.
%
%      PIPELINE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PIPELINE_GUI.M with the given input arguments.
%
%      PIPELINE_GUI('Property','Value',...) creates a new PIPELINE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pipeline_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pipeline_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pipeline_gui

% Last Modified by GUIDE v2.5 23-Mar-2014 12:30:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pipeline_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @pipeline_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before pipeline_gui is made visible.
function pipeline_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pipeline_gui (see VARARGIN)

% put data into handles
handles.AllData=varargin{1};
handles.AllPeaks=varargin{2};

% Setup first run parameters
handles.Status_Params.Current_Image=1;
handles.Status_Params.Current_Source=1;
handles.Status_Params.Source_Plot_Flag='BS';
handles.Status_Params.Selected_Background=0;
handles.Status_Params.First_Run=1;
handles.fig=figure;	

handles.Status_Params.annotation_handle=annotation('textbox',[0.15,0.8,.1,.1]);
set(handles.Status_Params.annotation_handle,'FontSize',15,'FontWeight','bold');

set(handles.BS,'value',1);
% set(handles.Source_List,'value',handles.Status_Params.Current_Source);

guidata(hObject,handles);
update_display(hObject, eventdata, handles);


uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = pipeline_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Write AllData and AllPeaks to the output vars
varargout{1} = handles.AllData;
varargout{2} = handles.AllPeaks;
figure1_CloseRequestFcn(hObject, eventdata, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
delete(hObject);
delete(handles.fig);
end


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
switch eventdata.Character
    case 'n'
        if handles.Status_Params.Current_Source < length(get(handles.Source_List,'String'))
            Next_Source_Callback(hObject, eventdata, handles);
        end
    case 'p'
        if handles.Status_Params.Current_Source > 1
            Previous_Source_Callback(hObject, eventdata, handles);
        end
    case 'r'
        Remove_Source_Callback(hObject, eventdata, handles);
    case 'N'
        Next_Image_Callback(hObject, eventdata, handles);
    case 'P'
	Previous_Image_Callback(hObject, eventdata, handles);
    case 'a'
        Add_Source_Callback(hObject, eventdata, handles);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Callback functions%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --- Executes on selection change in Source_List.
function Source_List_Callback(hObject, eventdata, handles)
% hObject    handle to Source_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Status_Params.Current_Source=get(hObject,'value');
guidata(hObject,handles);
update_display(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function Source_List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Source_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Add_Source.
function Add_Source_Callback(hObject, eventdata, handles)
% hObject    handle to Add_Source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(0, 'currentfigure', handles.fig);

[X,~,button]=ginput(1);

Current_Image=handles.Status_Params.Current_Image;

if isempty(handles.Peaks)
    Current_Source=1;
else
    Current_Source=length(handles.Peaks)+1;
end


[~,minindex]=min(abs(handles.Data.X-X));
% [~,minindex]=max(handles.Data.Ybs(handles.Data.X>=(X-10) & handles.Data.X<=(X+10)));

handles.Status_Params.Current_Source=Current_Source;
handles.Peaks(Current_Source).X=handles.Data.X(minindex);
handles.Peaks(Current_Source).OptimAperRad=[10 10];
handles.Peaks(Current_Source).Y=handles.Data.Y(handles.Peaks(Current_Source).X);
handles.Peaks(Current_Source).Ybs=handles.Data.Ybs(handles.Peaks(Current_Source).X);
X=handles.Peaks(Current_Source).X;
OptimAperRad=handles.Peaks(Current_Source).OptimAperRad;
handles.Peaks(Current_Source).Back=[[X-2*OptimAperRad(1),X-OptimAperRad(1)];[X+OptimAperRad(2),X+2*OptimAperRad(2)]];

[~,order]=sort([handles.Peaks.X]);
handles.Peaks=handles.Peaks(order);
handles.Status_Params.Current_Source=order(end);


handles.AllPeaks{Current_Image}=handles.Peaks;

guidata(hObject,handles);
update_display(hObject, eventdata, handles);


% --- Executes on button press in Remove_Source.
function Remove_Source_Callback(hObject, eventdata, handles)
% hObject    handle to Remove_Source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    Current_Image=handles.Status_Params.Current_Image;
    Current_Source=handles.Status_Params.Current_Source;

    if length(handles.Peaks)==1
        handles.Peaks(1)=[];
        handles.AllPeaks{Current_Image}=handles.Peaks;
        guidata(hObject,handles);
        update_display(hObject, eventdata, handles);
        return
    end



    handles.Peaks(Current_Source)=[];
    Source_List_Strings=get(handles.Source_List,'String');
    Source_List_Strings(end)=[];
    set(handles.Source_List,'String',Source_List_Strings);

    if Current_Source > length(handles.Peaks)
        handles.Status_Params.Current_Source=length(handles.Peaks);
    end
    
    
    handles.AllPeaks{Current_Image}=handles.Peaks;

    
    
    guidata(hObject,handles);
    update_display(hObject, eventdata, handles);

    
% --- Executes on button press in Next_Image.
function Next_Image_Callback(hObject, eventdata, handles)
% hObject    handle to Next_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
	NImages=length(handles.AllData);
	handles.Status_Params.Current_Image=handles.Status_Params.Current_Image+1;

	while isempty(handles.AllData{handles.Status_Params.Current_Image})
	handles.Status_Params.Current_Image=handles.Status_Params.Current_Image+1;
		if handles.Status_Params.Current_Image == NImages+1
		Previous_Image_Callback(hObject, eventdata, handles)
		return
		end
	end
		


    	handles.Status_Params.Current_Source=1;
        handles.Status_Params.Selected_Background=0;
        set(0, 'currentfigure', handles.fig);
        set(gca,'XLimMode','auto','YLimMode','auto');
	handles.Status_params.First_Run=1;
	guidata(hObject,handles);
	update_display(hObject, eventdata, handles);

    
% --- Executes on button press in Previous_Image.
function Previous_Image_Callback(hObject, eventdata, handles)
% hObject    handle to Previous_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
	
	handles.Status_Params.Current_Image=handles.Status_Params.Current_Image-1;

	while isempty(handles.AllData{handles.Status_Params.Current_Image})
	handles.Status_Params.Current_Image=handles.Status_Params.Current_Image-1;
		if handles.Status_Params.Current_Image == 0
		Next_Image_Callback(hObject, eventdata, handles)
		return
		end
	end



    	handles.Status_Params.Current_Source=1;
        handles.Status_Params.Selected_Background=0;
                handles.Status_Params.Selected_Background=0;
        set(0, 'currentfigure', handles.fig);
        set(gca,'XLimMode','auto','YLimMode','auto');
	handles.Status_params.First_Run=1;
	guidata(hObject,handles);
	update_display(hObject, eventdata, handles);


% --- Executes on button press in Next_Source.
function Next_Source_Callback(hObject, eventdata, handles)
% hObject    handle to Next_Source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    handles.Status_Params.Current_Source=handles.Status_Params.Current_Source+1;
    handles.Status_Params.Selected_Background=0;
    guidata(hObject,handles);
    update_display(hObject,eventdata,handles);


% --- Executes on button press in Previous_Source.
function Previous_Source_Callback(hObject, eventdata, handles)
% hObject    handle to Previous_Source (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.Status_Params.Current_Source=handles.Status_Params.Current_Source-1;
    handles.Status_Params.Selected_Background=0;
    guidata(hObject,handles);
    update_display(hObject,eventdata,handles);

% --- Executes on button press in BS.
function BS_Callback(hObject, eventdata, handles)
% hObject    handle to BS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BS
% 
    switch get(hObject,'value')
        case 0
            handles.Status_Params.Source_Plot_Flag='NoBS';
        case 1
            handles.Status_Params.Source_Plot_Flag='BS';
    end
    
guidata(hObject,handles);
update_display(hObject, eventdata, handles);

% --- Executes on button press in Replace_Background.
function Replace_Background_Callback(hObject, eventdata, handles)
% hObject    handle to Replace_Background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Current_Source=handles.Status_Params.Current_Source;
Selected_Background=handles.Status_Params.Selected_Background;
Current_Image=handles.Status_Params.Current_Image;
Back=handles.Peaks(Current_Source).Back;

if Selected_Background ~=0
    selection=getrect(handles.fig);
    Back(Selected_Background,:)=[selection(1),selection(1)+selection(3)];
    handles.Peaks(Current_Source).Back=Back;
    handles.AllPeaks{Current_Image}=handles.Peaks;
end

handles.Status_Params.Selected_Background = 0;

guidata(hObject,handles);
update_display(hObject, eventdata, handles);

% --- Executes on button press in Select_Background.
function Select_Background_Callback(hObject, eventdata, handles)
% hObject    handle to Select_Background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	guidata(hObject,handles);
		
	set(0, 'currentfigure', handles.fig);
	Selection=ginput(1);
		
	handles.Status_Params.Selected_Background_X=Selection(1);
	
	handles=Which_Background(hObject, eventdata, handles);

	guidata(hObject,handles);
    
    update_display(hObject, eventdata, handles);
    
    plot_invchildren(get(handles.fig,'CurrentAxes'));





% --- Executes on button press in Replace_Aperture.
function Replace_Aperture_Callback(hObject, eventdata, handles)
% hObject    handle to Replace_Aperture (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

	set(0, 'currentfigure', handles.fig);

%	v = getrect(handles.fig);
%	Current_Source=handles.Status_Params.Current_Source;
%	Current_Image=handles.Status_Params.Current_Image;
%	X=handles.Peaks(Current_Source).X;
%	OptimAperRad=[(X - v(1)), (v(1) + v(3) - X)];
%	handles.Peaks(Current_Source).OptimAperRad=OptimAperRad;
%	handles.AllPeaks{Current_Image}=handles.Peaks;

	v = ginput(2);
	v = sort(v);
	Current_Source=handles.Status_Params.Current_Source;
	Current_Image=handles.Status_Params.Current_Image;
	X=handles.Peaks(Current_Source).X;
	OptimAperRad=[(X - v(1)), (v(2) - X)];
	handles.Peaks(Current_Source).OptimAperRad=OptimAperRad;
	handles.AllPeaks{Current_Image}=handles.Peaks;

   	guidata(hObject,handles);
	update_display(hObject, eventdata, handles);



% --- Executes on button press in Close.
function Close_Callback(hObject, eventdata, handles)
% hObject    handle to Close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%AllData=handles.AllData;
%AllPeaks=handles.AllPeaks;
%
%if exist('PipelineOUT.mat','file')
%	disp('PipelineOUT.mat already exits!');
%else
%	save('PipelineOUT.mat','AllData','AllPeaks');
	close(handles.figure1);
%end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Assist functions %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% Plot sources
function plot_all_sources(hObject, eventdata, handles)
    Current_Source = handles.Status_Params.Current_Source;
    for ii=1:length(handles.Peaks)
        if ii == Current_Source
            continue
        end
        hold on
        X=handles.Peaks(ii).X;
        OptimAperRad=handles.Peaks(ii).OptimAperRad;
	if length(OptimAperRad)==1
		OptimAperRad=[OptimAperRad OptimAperRad];
	end
        switch handles.Status_Params.Source_Plot_Flag
            case 'BS'
                Y=min(handles.Data.Ybs(handles.Data.X < X+OptimAperRad(2) & handles.Data.X >= X-OptimAperRad(1)));
            case 'NoBS'
                Y=min(handles.Data.Y(handles.Data.X < X+OptimAperRad(2) & handles.Data.X >= X-OptimAperRad(1)));
        end
        hplot=plot([X-OptimAperRad(1),X+OptimAperRad(2)],[Y,Y],'r');
        set(hplot,'linewidth',2);
        hold off
    end
    
% Plot data
function plot_data(hObject, eventdata, handles)
%	XLim=[min(handles.Data.X),max(handles.Data.X)];
    switch handles.Status_Params.Source_Plot_Flag
        case 'BS'
            plot(handles.Data.X,handles.Data.Ybs,'b');
%	YLim=[min(handles.Data.Ybs),max(handles.Data.Ybs)];
        case 'NoBS'
            plot(handles.Data.X,handles.Data.Y,'b');
%	YLim=[min(handles.Data.Y),max(handles.Data.Y)];
    end



% Plot source
function plot_source(hObject, eventdata, handles)
    hold on
    Source_Num=handles.Status_Params.Current_Source;
    X=handles.Peaks(Source_Num).X;
    OptimAperRad=handles.Peaks(Source_Num).OptimAperRad;
	if  length(OptimAperRad)==1
		OptimAperRad=[OptimAperRad OptimAperRad];
	end
    switch handles.Status_Params.Source_Plot_Flag
        case 'BS'
            Y=min(handles.Data.Ybs(handles.Data.X < X+OptimAperRad(2) & handles.Data.X >= X-OptimAperRad(1)));
        case 'NoBS'
            Y=min(handles.Data.Y(handles.Data.X < X+OptimAperRad(2) & handles.Data.X >= X-OptimAperRad(1)));
    end
    hplot=plot([X-OptimAperRad(1),X+OptimAperRad(2)],[Y,Y],'g');
    set(hplot,'linewidth',3,'color',[0,0.8,0]);
    hold off

%Plot source background region
function plot_source_background(hObject, eventdata, handles)
    hold on
    Source_Num=handles.Status_Params.Current_Source;
%     X=handles.Peaks(Source_Num).X;
%     DX=handles.Peaks(Source_Num).DX;
    switch handles.Status_Params.Source_Plot_Flag
        case 'BS'
            Y=handles.Data.Ybs;
        case 'NoBS'
            Y=handles.Data.Y;
    end
    

    X=handles.Data.X;
    
    for ii=1:2
        X1=handles.Peaks(Source_Num).Back(ii,1);
        X2=handles.Peaks(Source_Num).Back(ii,2);
        indtemp=X>=X1 & X<=X2;
        ind=[find(indtemp,1,'first'),find(indtemp,1,'last')];
        switch ii
            case handles.Status_Params.Selected_Background
                hplot=plot(X(indtemp),Y(indtemp),'m');
                set(hplot,'linewidth',10,'MarkerSize',10);
            otherwise
                hplot=plot(X(ind),Y(ind),'xblack');
                set(hplot,'linewidth',2,'MarkerSize',10);
        end
    end
    hold off
    
function write_textbox(hObject, eventdata, handles)
    Current_Source=handles.Status_Params.Current_Source;
    Peak=handles.Peaks(Current_Source);
    handles.Text_Box_String={...
        ['X = ' num2str(Peak.X)],...
        ['Ap = ' num2str(round(Peak.OptimAperRad))],...
        ['Y = ' num2str(Peak.Y)],...
        ['Ybs = ' num2str(Peak.Ybs)],...
	['Background ranges:'],...
	[num2str(round(Peak.Back(1,:)))],...
	[num2str(round(Peak.Back(2,:)))]};
    set(handles.Data_Display,'String',handles.Text_Box_String);
    guidata(hObject,handles);
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Update display%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Updates the display according  to the paramaeters in Status_Params

function update_display(hObject, eventdata, handles)

    ReleaseFocus(handles.figure1);
    
    %Get zoom settings of figure
    set(0, 'currentfigure', handles.fig);
    XLim=get(gca,'XLim');
    YLim=get(gca,'YLim');
   
    %set current image
    handles.Data=handles.AllData{handles.Status_Params.Current_Image};
    handles.Peaks=handles.AllPeaks{handles.Status_Params.Current_Image};



	NSources = length(handles.Peaks);
	NImages = length(handles.AllData);
	if isempty(handles.Data)
	Object_Name = '';
	else
	Object_Name = handles.Data.Object{1};
	end
	Current_Source = handles.Status_Params.Current_Source;
	Current_Image = handles.Status_Params.Current_Image;
	Source_String = {['Image ' num2str(Current_Image) '/' num2str(NImages)],...
			['Object name: ' Object_Name],...
			['Source ' num2str(Current_Source) '/' num2str(NSources)]...
			};

	set(handles.Status_Params.annotation_handle,'String',Source_String,'FontSize',20,'FontWeight','bold');



	if isempty(handles.Data)
		handles.Status_Params.Current_Image=handles.Status_Params.Current_Image+1;
        	guidata(hObject,handles);
		update_display(hObject, eventdata, handles);
		return
	end

    if isempty(handles.Peaks)
%         handles.Status_Params.Current_Source=0;
        plot_data(hObject, eventdata, handles);
        set(handles.Source_List,'Enable','off');
        set(handles.Next_Source,'Enable','off');
        set(handles.Previous_Source,'Enable','off');
        set(handles.Remove_Source,'Enable','off');
        set(handles.Select_Background,'Enable','off');
        set(handles.Replace_Background,'Enable','off');
        set(handles.Replace_Aperture,'Enable','off');
        switch handles.Status_Params.Current_Image 
        case 1
            set(handles.Previous_Image,'Enable','off');
            set(handles.Next_Image,'Enable','on');
        case length(handles.AllPeaks) 
            set(handles.Next_Image,'Enable','off');
            set(handles.Previous_Image,'Enable','on'); 
        otherwise
            set(handles.Next_Image,'Enable','on');
            set(handles.Previous_Image,'Enable','on');
        end
        guidata(hObject,handles);
        return
    end
    
%     handles.Status_Params.Current_Source=1;
    
    set(handles.Next_Source,'Enable','on');
    set(handles.Select_Background,'Enable','on');
    set(handles.Replace_Background,'Enable','on');
    set(handles.Remove_Source,'Enable','on');
    % Populate source list
    for ii=1:length(handles.Peaks)
        Source_Strings{ii}=['Source ' num2str(ii)];
    end

    set(handles.Source_List,'String',Source_Strings,'value',handles.Status_Params.Current_Source,'Enable','on');

%     set(handles.Source_List,'value',handles.Status_Params.Current_Source);

    set(0, 'currentfigure', handles.fig);

%	if length(handles.Peaks(handles.Status_Params.Current_Source).OptimAperRad) == 1
%		handles.Peaks(handles.Status_Params.Current_Source).OptimAperRad=[handles.Peaks(handles.Status_Params.Current_Source).OptimAperRad(1) handles.Peaks(handles.Status_Params.Current_Source).OptimAperRad(1)];
%		handles.AllPeaks{handles.Status_Params.Current_Image}=handles.Peaks;
%		guidata(hObject, handles);
%	end

    plot_data(hObject, eventdata, handles);
    plot_all_sources(hObject, eventdata, handles);
    plot_source(hObject, eventdata, handles);
    plot_source_background(hObject, eventdata, handles);
    write_textbox(hObject, eventdata, handles);
    

    
    
    switch handles.Status_Params.Current_Source 
        case 1
            set(handles.Previous_Source,'Enable','off');
            set(handles.Next_Source,'Enable','on');
        case length(get(handles.Source_List,'String')) 
            set(handles.Next_Source,'Enable','off');
            set(handles.Previous_Source,'Enable','on'); 
        otherwise
            set(handles.Next_Source,'Enable','on');
            set(handles.Previous_Source,'Enable','on');
    end
    
    if length(get(handles.Source_List,'String'))==1
            set(handles.Next_Source,'Enable','off');
            set(handles.Previous_Source,'Enable','off'); 
    end
    
    
    
        switch handles.Status_Params.Current_Image 
        case 1
            set(handles.Previous_Image,'Enable','off');
            set(handles.Next_Image,'Enable','on');
        case length(handles.AllPeaks) 
            set(handles.Next_Image,'Enable','off');
            set(handles.Previous_Image,'Enable','on'); 
        otherwise
            set(handles.Next_Image,'Enable','on');
            set(handles.Previous_Image,'Enable','on');
        end
    
        
    switch handles.Status_Params.Selected_Background
        case 0
            set(handles.Replace_Background,'Enable','off');
        otherwise
            set(handles.Replace_Background,'Enable','on');
    end
    
    if length(handles.AllPeaks)==1
            set(handles.Next_Image,'Enable','off');
            set(handles.Previous_Image,'Enable','off'); 
    end
    
    set(0, 'currentfigure', handles.fig);
    if handles.Status_Params.First_Run==0
    set(gca,'XLim',XLim);
    set(gca,'YLim',YLim);
    end
    

    handles.Status_Params.First_Run=0;
    guidata(hObject,handles);
    
    figure(handles.figure1);
    
function [handles] = Which_Background(hObject, eventdata, handles)


	Current_Source=handles.Status_Params.Current_Source;
	X_Selected=handles.Status_Params.Selected_Background_X;
	Back=handles.Peaks(Current_Source).Back;
	handles.Status_Params.Selected_Background=0;
	
	Overlap = X_Selected >= Back(:,1) & X_Selected <= Back(:,2);

	if ~any(Overlap)
		disp('WARNING: No background selected.');
		handles.Status_Params.Selected_Background=0;
	else	
		handles.Status_Params.Selected_Background=find(Overlap);
    end
    
function ReleaseFocus(fig)
        set(findobj(fig, 'Type', 'uicontrol'), 'Enable', 'off');
        drawnow;
        set(findobj(fig, 'Type', 'uicontrol'), 'Enable', 'on');


% --- Executes on selection change in Survey_List.
function Survey_List_Callback(hObject, eventdata, handles)
% hObject    handle to Survey_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Survey_List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Survey_List

Survey_String_List = get(handles.Survey_List,'String');
Survey_List_Value = get(handles.Survey_List,'value');

switch Survey_String_List{Survey_List_Value}
    case '2mass'
        set(handles.Filter_List,'String',{'j','h','k'});
    case 'dssstsci'
        set(handles.Filter_List,'String',{'poss2ukstu_red','poss2ukstu_ir','poss2ukstu_blue',  'poss1_blue,poss1_red','all','quickv','phase2_gsc2','phase2_gsc1'});
    case 'dsseso'
        set(handles.Filter_List,'String',  {'DSS1','DSS2-red','DSS2-blue','DSS2-infrared'});
    case 'skyview'
        set(handles.Filter_List,'String',{'sdssi','sdssr','sdssg','sdssu','sdssg'});
end


% --- Executes during object creation, after setting all properties.
function Survey_List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Survey_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

Survey_Name  = {'2mass','dssstsci','dsseso','skyview'} ;
set(hObject,'String',Survey_Name);
guidata(hObject, handles);



% --- Executes on selection change in Filter_List.
function Filter_List_Callback(hObject, eventdata, handles)
% hObject    handle to Filter_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Filter_List contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Filter_List


% --- Executes during object creation, after setting all properties.
function Filter_List_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Filter_List (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

Survey_Filter   =  {'j','h','k'};
set(hObject,'String',Survey_Filter);
guidata(hObject, handles);




% --- Executes on button press in Display_Slit.
function Display_Slit_Callback(hObject, eventdata, handles)
% hObject    handle to Display_Slit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

RAD = 180./pi;
Current_Image=handles.Status_Params.Current_Image;

RA=handles.AllData{Current_Image}.RA{1};
Dec=handles.AllData{Current_Image}.Dec{1};
PA=handles.AllData{Current_Image}.PA{1};


Survey_Name_List=get(handles.Survey_List,'String');
Survey_Name_Value=get(handles.Survey_List,'value');
Filter_Name_List=get(handles.Filter_List,'String');
Filter_Name_Value=get(handles.Filter_List,'value');

SurveyName=Survey_Name_List{Survey_Name_Value};
SurveyFilter=Filter_Name_List{Filter_Name_Value};


ds9_slit(RA./RAD,Dec./RAD,'PA',PA,'StartDS9',true,'Survey',SurveyName,'Filter',SurveyFilter)
