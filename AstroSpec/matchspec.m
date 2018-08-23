function varargout = matchspec(varargin)
%--------------------------------------------------------------------------
% matchspec function                                             AstroSpec
% Description: A GUI utility to inspect and match spectrum with templates.
% Input  : - Spectra to inspect [Wavelength, Intensity, Error]
%            Error is optional.
%          - Method to plot spectra: {'stairs' | 'plot' | 'errorxy'},
%            default is 'stairs'.
%          - Template spectra [Wavelength, Intensity, Error].
%          - Template color, default is 'r'.
% Output : null
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jul 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Spec=get_spectra('Gal_S0'); matchspec(Spec);
% Reliable: 2
%--------------------------------------------------------------------------

% MATCHSPEC M-file for matchspec.fig
%      MATCHSPEC, by itself, creates a new MATCHSPEC or raises the existing
%      singleton*.
%
%      H = MATCHSPEC returns the handle to a new MATCHSPEC or the handle to
%      the existing singleton*.
%
%      MATCHSPEC('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MATCHSPEC.M with the given input arguments.
%
%      MATCHSPEC('Property','Value',...) creates a new MATCHSPEC or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before matchspec_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to matchspec_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help matchspec

% Last Modified by GUIDE v2.5 17-Jul-2005 19:25:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @matchspec_OpeningFcn, ...
                   'gui_OutputFcn',  @matchspec_OutputFcn, ...
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

%--------------------------------------------------------------------
ColW  = 1;
ColI  = 2;
ColE  = 3;

%---------------------------------------------
%--- Plot Intial Screen and save variables ---
%---------------------------------------------
Narg = length(varargin);
if (Narg>0 && size(varargin{1},1)>1),

    Spec = varargin{1};
    
   if (Narg==1),
      PlotMethod    = 'stairs';
      Template      = [];
      TemplateColor = 'r';
   elseif (Narg==2),
      PlotMethod    = varargin{2};
      Template      = [];
      TemplateColor = 'r';
   elseif (Narg==3),
      PlotMethod    = varargin{2};
      Template      = varargin{3};
      TemplateColor = 'r';
   elseif (Narg==4),
      PlotMethod    = varargin{2};
      Template      = varargin{3};
      TemplateColor = varargin{4};
      % do nothing
   else
      error('Illegal number of input arguments');
   end

   %-----------------
   %--- Plot Spec ---
   %-----------------
   Hfig = figure;
   hold on;
   box on;
   set(gca,'FontSize',12);
   H = xlabel('Wavelength');
   set(H,'FontSize',14);
   H = ylabel('Flux');
   set(H,'FontSize',14);


   switch PlotMethod
    case 'plot'
       Hspec = plot(Spec(:,ColW),Spec(:,ColI));
    case 'stairs'
       Hspec = stairs(Spec(:,ColW),Spec(:,ColI));
    case 'errorxy'
       if (size(Spec,2)>2),
          errorxy(Spec);
       else
          disp('Spectrum doesnot contain error - cannot use errorxy');
       end

    otherwise
       disp('Unknown PlotMethod Option');
       return;
   end

   %---------------------
   %--- Plot Template ---
   %---------------------
   if (isempty(Template)),
      % do nothing
      Htemp = [];
   else
      Htemp = stairs(Template(:,ColW),Template(:,ColI),TemplateColor);

      %set(findobj(gcbf,'Tag','pushbutton_Color'),'ForegroundColor',TemplateColor);

   end


   %------------------------------
   %--- Store data in UserData ---
   %------------------------------
   UserData.Spec          = Spec;
   UserData.Template{1}   = Template;
   if (isempty(Htemp)),
      UserData.Ntemp      = 0;
   else
      UserData.Ntemp      = 1;
   end
   UserData.Hspec         = Hspec;
  
   UserData.Htemp{1}      = Htemp;
   UserData.CurTemp       = 1;
   UserData.ScaleTemp{1}  = [0 1];

   set(matchspec,'UserData',UserData);


   %------------------------
   %--- Set Hold On Axis ---
   %------------------------
   %HoldOnAxis = get(findobj(gcbf,'Tag','checkbox_HoldOnAxis'),'Value');
   %if (isempty(HoldOnAxis)),
   %   % do nothing
   %else
   %   switch HoldOnAxis
   %    case 0
   %       % do not hold on axis
   %       set(gca,'XLimMode','Auto');
   %       set(gca,'YLimMode','Auto');
   %    case 1
   %       % hold on axis
   %       set(gca,'XLimMode','Manual');
   %       set(gca,'YLimMode','Manual');
   %    other wise
   %       error('Unknown HoldOnAxis Option');
   %   end
   %end
end



% --- Executes just before matchspec is made visible.
function matchspec_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to matchspec (see VARARGIN)

% Choose default command line output for matchspec
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes matchspec wait for user response (see UIRESUME)
% uiwait(handles.figure1);




% --- Outputs from this function are returned to the command line.
function varargout = matchspec_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_Redshift_Callback(hObject, eventdata, handles)
% hObject    handle to slider_Redshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


%--- Read Min Redshift ---
MinRedshift = str2num(get(findobj(gcbf,'Tag','edit_MinRedshift'),'String'));

%--- Read Max Redshift ---
MaxRedshift = str2num(get(findobj(gcbf,'Tag','edit_MaxRedshift'),'String'));

%--- Get redshift scoler position ---
RedshiftScrolerPos = get(hObject,'Value');
Redshift           = (MaxRedshift - MinRedshift).*RedshiftScrolerPos + MinRedshift;

%--- set redshift in redshift text window ---
set(findobj(gcbf,'Tag','edit_Redshift'),'String',sprintf('%8.5f',Redshift));

redisplay_spec;






% --- Executes during object creation, after setting all properties.
function slider_Redshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_Redshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_MinRedshift_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MinRedshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MinRedshift as text
%        str2double(get(hObject,'String')) returns contents of edit_MinRedshift as a double


% --- Executes during object creation, after setting all properties.
function edit_MinRedshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MinRedshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_MaxRedshift_Callback(hObject, eventdata, handles)
% hObject    handle to edit_MaxRedshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_MaxRedshift as text
%        str2double(get(hObject,'String')) returns contents of edit_MaxRedshift as a double


% --- Executes during object creation, after setting all properties.
function edit_MaxRedshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_MaxRedshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Redshift_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Redshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Redshift as text
%        str2double(get(hObject,'String')) returns contents of edit_Redshift as a double

redisplay_spec;





% --- Executes during object creation, after setting all properties.
function edit_Redshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Redshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_Scaling.
function popupmenu_Scaling_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Scaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_Scaling contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Scaling


redisplay_spec('Method');




% --- Executes during object creation, after setting all properties.
function popupmenu_Scaling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_Scaling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String', {'None','Mean','Median','StD','Range','Min','Max'});


% --- Executes on selection change in popupmenu_ZeroShift.
function popupmenu_ZeroShift_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_ZeroShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_ZeroShift contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_ZeroShift

redisplay_spec('Method');


    
    
    

% --- Executes during object creation, after setting all properties.
function popupmenu_ZeroShift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_ZeroShift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String', {'None','Mean','Median','StD','Min','Max','Fit'});


function edit_ShiftVal_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ShiftVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ShiftVal as text
%        str2double(get(hObject,'String')) returns contents of edit_ShiftVal as a double

redisplay_spec;


% --- Executes during object creation, after setting all properties.
function edit_ShiftVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ShiftVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ScaleVal_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ScaleVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ScaleVal as text
%        str2double(get(hObject,'String')) returns contents of edit_ScaleVal as a double

redisplay_spec;


% --- Executes during object creation, after setting all properties.
function edit_ScaleVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ScaleVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_MarkLines.
function popupmenu_MarkLines_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_MarkLines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_MarkLines contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_MarkLines


% --- Executes during object creation, after setting all properties.
function popupmenu_MarkLines_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_MarkLines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%set(hObject,'String', {'QSO Emission','Galaxy','SN-Ia','SN-II','Stellar Late','Stellar Early'});


% --- Executes on selection change in popupmenu_LoadTemplate.
function popupmenu_LoadTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_LoadTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_LoadTemplate contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_LoadTemplate

ColW = 1;
ColI = 2;

TemplateOption = get(findobj(gcbf,'Tag','popupmenu_LoadTemplate'),'Value');
TemplateName   = index2templatename(TemplateOption);

Template       = get_spectra(TemplateName,[0 0],[0 0],0,0);

UserData       = get(matchspec,'UserData');

UserData.Ntemp = UserData.Ntemp + 1;
UserData.CurTemp = UserData.Ntemp;
UserData.Template{UserData.Ntemp} = Template;

   
%---------------------
%--- Plot Template ---
%---------------------
TemplateColor = get(findobj(gcbf,'Tag','pushbutton_Color'),'ForegroundColor');


Hfig = get(get(UserData.Hspec,'Parent'),'Parent');
figure(Hfig);
Htemp         = stairs(Template(:,ColW),Template(:,ColI));

set(Htemp,'Color',TemplateColor);

UserData.Htemp{UserData.Ntemp} = Htemp;
set(matchspec,'UserData',UserData);




% --- Executes during object creation, after setting all properties.
function popupmenu_LoadTemplate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_LoadTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String', {'QSO_LBQS','QSO_FBQS','QSO_FBQS_RL','QSO_FBQS_RQ','QSO_SDSS','QSO_HBal','QSO_LBal',...
    'Gal_E','Gal_S0','Gal_Sa','Gal_Sb','Gal_Sc','Gal_StarBurst1','Gal_Bulge',...
    'ukk0v','ukk2v','ukk3v','ukk4v','ukk5v','ukk7v','ukk0iii','ukk1iii','ukk2iii','ukk3iii','ukk4iii','ukk5iii'});

% 'O5V','O7B0V','B34V','B6V','A13V','A57V','A8V','F67V',...
%    'F89V','G12V','G68V','K4V','K5V','M2V','M4V','O8I','B1I','A03I','A79I','F7I','G01I','K35I',...

% --- Executes on button press in pushbutton_CCF.
function pushbutton_CCF_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_CCF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%--- Read CCF parameters ---
CCF_Step    = str2num(get(findobj(gcbf,'Tag','edit_CCF_Step'),'String'));

%--- Read Min Redshift ---
MinRedshift = str2num(get(findobj(gcbf,'Tag','edit_MinRedshift'),'String'));

%--- Read Max Redshift ---
MaxRedshift = str2num(get(findobj(gcbf,'Tag','edit_MaxRedshift'),'String'));

%--- Read ShiftMethod ---
ShiftOption = get(findobj(gcbf,'Tag','popupmenu_ZeroShift'),'Value');
ShiftMethod = index2shiftmethod(ShiftOption);

%--- Read ScaleMethod ---
ScaleOption = get(findobj(gcbf,'Tag','popupmenu_Scaling'),'Value');
ScaleMethod = index2scalemethod(ScaleOption);

SigmaClip            = 2.5;

%--- Read Poly ---
Poly        = str2num(get(findobj(gcbf,'Tag','edit_PolyDivide'),'String'));


%--- Get UserData ---
UserData = get(matchspec,'UserData');  
CurTemp  = UserData.CurTemp;
Template = UserData.Template{CurTemp};
Htemp    = UserData.Htemp{CurTemp};
Spec     = UserData.Spec;
Hspec    = UserData.Hspec;


if (isempty(Htemp) | isempty(Hspec)),
   disp(sprintf('Need both spectra and template in order to cross correlate'))
else
   RedshiftVec          = [MinRedshift:CCF_Step:MaxRedshift].';
   [CC,Npt,RedshiftVec] = ccf_spectra(Spec,Template,RedshiftVec,ShiftMethod,ScaleMethod,SigmaClip,Poly);
   
   %--- assignin variables ---
   disp(sprintf('Cross correlation results are available in session variables: '))
   disp(sprintf('     matchspec_CC'))
   disp(sprintf('     matchspec_Npt'))
   disp(sprintf('     matchspec_RedshiftVec'))
   disp(sprintf('Check out ccf_spectra.m for more details'))

   assignin('base','matchspec_CC',CC);
   assignin('base','matchspec_Npt',Npt);
   assignin('base','matchspec_RedshiftVec',RedshiftVec);
   
   
   figure;
   plot(RedshiftVec, CC);
   set(gca,'FontSize',14);
   H = xlabel('Redshift');
   set(H,'FontSize',18);
   H = ylabel('Cross-Correlation');
   set(H,'FontSize',18);
end








function edit_SourceE_BV_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SourceE_BV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SourceE_BV as text
%        str2double(get(hObject,'String')) returns contents of edit_SourceE_BV as a double

redisplay_spec;



% --- Executes during object creation, after setting all properties.
function edit_SourceE_BV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SourceE_BV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SourceR_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SourceR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SourceR as text
%        str2double(get(hObject,'String')) returns contents of edit_SourceR as a double

redisplay_spec;




% --- Executes during object creation, after setting all properties.
function edit_SourceR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SourceR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_BlackBodyTemp_Callback(hObject, eventdata, handles)
% hObject    handle to edit_BlackBodyTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_BlackBodyTemp as text
%        str2double(get(hObject,'String')) returns contents of edit_BlackBodyTemp as a double


% --- Executes during object creation, after setting all properties.
function edit_BlackBodyTemp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_BlackBodyTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ObserverE_BV_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ObserverE_BV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ObserverE_BV as text
%        str2double(get(hObject,'String')) returns contents of edit_ObserverE_BV as a double

redisplay_spec;




% --- Executes during object creation, after setting all properties.
function edit_ObserverE_BV_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ObserverE_BV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ObserverR_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ObserverR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ObserverR as text
%        str2double(get(hObject,'String')) returns contents of edit_ObserverR as a double

redisplay_spec;




% --- Executes during object creation, after setting all properties.
function edit_ObserverR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ObserverR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_BlackBodyTemp.
function pushbutton_BlackBodyTemp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_BlackBodyTemp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ColW = 1;
ColI = 2;
WaveVec = [1000:10:20000].';

Temp = str2num(get(findobj(gcbf,'Tag','edit_BlackBodyTemp'),'String'));
Template  = [WaveVec, black_body(Temp,WaveVec)];

UserData       = get(matchspec,'UserData');

UserData.Ntemp = UserData.Ntemp + 1;
UserData.CurTemp = UserData.Ntemp;
UserData.Template{UserData.Ntemp} = Template;

   
%---------------------
%--- Plot Template ---
%---------------------
TemplateColor = get(findobj(gcbf,'Tag','pushbutton_Color'),'ForegroundColor');


Hfig = get(get(UserData.Hspec,'Parent'),'Parent');
figure(Hfig);
Htemp         = stairs(Template(:,ColW),Template(:,ColI));

set(Htemp,'Color',TemplateColor);

UserData.Htemp{UserData.Ntemp} = Htemp;
set(matchspec,'UserData',UserData);







% --- Executes on button press in popupmenu_LoadSNTemplate.
function popupmenu_LoadSNTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_LoadSNTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ColW = 1;
ColI = 2;

TemplateOption = get(findobj(gcbf,'Tag','popupmenu_LoadSNTemplate'),'Value');
TemplateName   = index2sn_templatename(TemplateOption);

[Template,SN_AgeRange] = get_spectra(TemplateName,[0 0],[0 0],0,0);
set(findobj(gcbf,'Tag','SN_MinAge'),'String',num2str(SN_AgeRange(1)));
set(findobj(gcbf,'Tag','SN_MaxAge'),'String',num2str(SN_AgeRange(2)));

UserData       = get(matchspec,'UserData');

UserData.Ntemp = UserData.Ntemp + 1;
UserData.CurTemp = UserData.Ntemp;
UserData.Template{UserData.Ntemp} = Template;

   
%------------------------
%--- Plot SN Template ---
%------------------------
TemplateColor = get(findobj(gcbf,'Tag','pushbutton_Color'),'ForegroundColor');


Hfig = get(get(UserData.Hspec,'Parent'),'Parent');
figure(Hfig);
Htemp         = stairs(Template(:,ColW),Template(:,ColI));

set(Htemp,'Color',TemplateColor);

UserData.Htemp{UserData.Ntemp} = Htemp;
set(matchspec,'UserData',UserData);







% --- Executes during object creation, after setting all properties.
function popupmenu_LoadSNTemplate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_LoadSNTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


%set(hObject,'String',{'SN1984l','SN1987l','SN1990b','SN1990u','SN1991ar','SN1991bg','SN1991t','SN1992h','SN1993j','SN1994ak','SN1994d','SN1994i','SN1994y','SN1995d','SN1996cb','SN1998bp','SN1998bw','SN1998de','SN1998dt','SN1998es','SN1998s','SN1999ee','SN1999em','SN2000cx','SN2001x','SN2002ap','SN1996L','SN1997cy','SN1998S','SN1999ex'});
set(hObject,'String',{''}); %SN1984l','SN1987l','SN1990b','SN1990u','SN1991ar','SN1991bg','SN1991t','SN1992h','SN1993j','SN1994ak','SN1994d','SN1994i','SN1994y','SN1995d','SN1996cb','SN1998bp','SN1998bw','SN1998de','SN1998dt','SN1998es','SN1998s','SN1999ee','SN1999em','SN2000cx','SN2001x','SN2002ap','SN1996L','SN1997cy','SN1998S','SN1999ex'});




% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in checkbox_HoldOnAxis.
function checkbox_HoldOnAxis_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_HoldOnAxis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_HoldOnAxis

HoldOnAxis = get(findobj(gcbf,'Tag','checkbox_HoldOnAxis'),'Value');
switch HoldOnAxis
 case 0
    % do not hold on axis
    set(gca,'XLimMode','Auto');
    set(gca,'YLimMode','Auto');
 case 1
    % hold on axis
    set(gca,'XLimMode','Manual');
    set(gca,'YLimMode','Manual');
 otherwise
    error('Unknown HoldOnAxis Option');
end



% --- Executes on button press in checkbox_HoldOnTemplate.
function checkbox_HoldOnTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_HoldOnTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_HoldOnTemplate



function edit_SNAge_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SNAge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SNAge as text
%        str2double(get(hObject,'String')) returns contents of edit_SNAge as a double


%--- Change SN Age accordingly ---
SN_Age = str2num(get(findobj(gcbf,'Tag','edit_SNAge'),'String'));

%--- Get SN Age ---

ColW = 1;
ColI = 2;

TemplateOption = get(findobj(gcbf,'Tag','popupmenu_LoadSNTemplate'),'Value');
TemplateName   = index2sn_templatename(TemplateOption);

[Template,SN_AgeRange] = get_spectra(TemplateName,[0 0],[0 0],0,0,SN_Age);
%set(findobj(gcbf,'Tag','SN_MinAge'),'String',num2str(SN_AgeRange(1)));
%set(findobj(gcbf,'Tag','SN_MaxAge'),'String',num2str(SN_AgeRange(2)));

UserData       = get(matchspec,'UserData');

%UserData.Ntemp = UserData.Ntemp + 1;
UserData.CurTemp = UserData.Ntemp;
UserData.Template{UserData.Ntemp} = Template;

set(matchspec,'UserData',UserData);

%--- Redisplay spectra ---
redisplay_spec;





% --- Executes during object creation, after setting all properties.
function edit_SNAge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SNAge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_GetFromCoo.
function pushbutton_GetFromCoo_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_GetFromCoo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_PolyDivide_Callback(hObject, eventdata, handles)
% hObject    handle to edit_PolyDivide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_PolyDivide as text
%        str2double(get(hObject,'String')) returns contents of edit_PolyDivide as a double

redisplay_spec;





% --- Executes during object creation, after setting all properties.
function edit_PolyDivide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_PolyDivide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes on button press in pushbutton_Color.
function pushbutton_Color_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


TemplateColor = uisetcolor;
%--- Set current template color ---
set(findobj(gcbf,'Tag','pushbutton_Color'),'ForegroundColor',TemplateColor);

UserData = get(matchspec,'UserData');
CurTemp  = UserData.CurTemp;
Htemp    = UserData.Htemp{CurTemp};
set(Htemp,'Color',TemplateColor);






% --- Executes on selection change in popupmenu_LineType.
function popupmenu_LineType_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_LineType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_LineType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_LineType

%--- Read LineType ---
LineTypeOption = get(findobj(gcbf,'Tag','popupmenu_LineType'),'Value');
LineType       = index2linetype(LineTypeOption);

%--- Set current template LineType ---
UserData = get(matchspec,'UserData');
CurTemp  = UserData.CurTemp;
Htemp    = UserData.Htemp{CurTemp};
set(Htemp,'LineStyle',LineType);






% --- Executes during object creation, after setting all properties.
function popupmenu_LineType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_LineType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String', {'-','--','-.',':'});




% --- Executes on button press in radiobutton_RedshiftTemplate.
function radiobutton_RedshiftTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_RedshiftTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_RedshiftTemplate






%-------------------------------------------
% Convert index to LineType
% 1 - '-'
% 2 - '--'
% 3 - '-.'
% 4 - ':'
%-------------------------------------------
function LineType = index2linetype(Index);

switch Index
 case 1
    LineType = '-';
 case 2
    LineType = '--';
 case 3
    LineType = '-.';
 case 4
    LineType = ':';
 otherwise
    error('Unknown Index (LineType) Option');
end
       

%-------------------------------------------
% Convert index to ShiftMethod
% 1 - 'none'
% 2 - 'mean'
% 3 - 'median'
% 4 - 'std'
% 5 - 'min'
% 6 - 'max'
% 7 - 'fit'
%-------------------------------------------
function ShiftMethod = index2shiftmethod(Index);

switch Index
 case 1
    ShiftMethod = 'none';
 case 2
    ShiftMethod = 'mean';
 case 3
    ShiftMethod = 'median';
 case 4
    ShiftMethod = 'std';
 case 5
    ShiftMethod = 'min';
 case 6
    ShiftMethod = 'max';
 case 7
    ShiftMethod = 'fit';
 otherwise
    error('Unknown Index (ShiftMethod) Option');
end
       



%-------------------------------------------
% Convert index to ScaleMethod
% 1 - 'none'
% 2 - 'mean'
% 3 - 'median'
% 4 - 'range'
% 5 - 'std'
% 6 - 'min'
% 7 - 'max'
% 8 - 'fit'
%-------------------------------------------
function ScaleMethod = index2scalemethod(Index);

switch Index
 case 1
    ScaleMethod = 'none';
 case 2
    ScaleMethod = 'mean';
 case 3
    ScaleMethod = 'median';
 case 4
    ScaleMethod = 'std';
 case 5
    ScaleMethod = 'range';
 case 6
    ScaleMethod = 'min';
 case 7
    ScaleMethod = 'max';
 otherwise
    error('Unknown Index (ScaleMethod) Option');
end



%-------------------------------------------
% Convert index to TemplateName
%-------------------------------------------
function TemplateName = index2templatename(Index);

switch Index
 case 1
    TemplateName = 'QSO_LBQS';
 case 2
    TemplateName = 'QSO_FBQS';
 case 3
    TemplateName = 'QSO_FBQS_RL';
 case 4
    TemplateName = 'QSO_FBQS_RQ';
 case 5
    TemplateName = 'QSO_SDSS';
 case 6
    TemplateName = 'QSO_HBal';
 case 7
    TemplateName = 'QSO_LBal';
 case 8
    TemplateName = 'Gal_E';
 case 9
    TemplateName = 'Gal_S0';
 case 10
    TemplateName = 'Gal_Sa';
 case 11
    TemplateName = 'Gal_Sb';
 case 12
    TemplateName = 'Gal_Sc';
 case 13
    TemplateName = 'Gal_StarBurst1';
 case 14
    TemplateName = 'Gal_Bulge';
 case 15
    TemplateName = 'ukk0v';
 case 16
    TemplateName = 'ukk2v';
 case 17
    TemplateName = 'ukk3v';
 case 18
    TemplateName = 'ukk4v';
 case 19
    TemplateName = 'ukk5v';
 case 20
    TemplateName = 'ukk7v';
 case 21
    TemplateName = 'ukk0iii';   
 case 22
    TemplateName = 'ukk1iii';   
 case 23
    TemplateName = 'ukk2iii';   
 case 24
    TemplateName = 'ukk3iii'; 
 case 25
    TemplateName = 'ukk4iii';   
 case 26
    TemplateName = 'ukk5iii';    
    
 otherwise
    error('Unknonw index of TemplateName Option');
end



%-------------------------------------------
% Convert index to SN TemplateName
%-------------------------------------------
function TemplateName = index2sn_templatename(Index);

switch Index
 case 1
    TemplateName = 'SN1984l';
 case 2
    TemplateName = 'SN1987l';
 case 3
    TemplateName = 'SN1990b';
 case 4
    TemplateName = 'SN1990u';
 case 5
    TemplateName = 'SN1991ar';
 case 6
    TemplateName = 'SN1991bg';
 case 7
    TemplateName = 'SN1991t';
 case 8
    TemplateName = 'SN1992h';
 case 9
    TemplateName = 'SN1993j';
 case 10
    TemplateName = 'SN1994ak';
 case 11
    TemplateName = 'SN1994d';
 case 12
    TemplateName = 'SN1994i';
 case 13
    TemplateName = 'SN1994y';
 case 14
    TemplateName = 'SN1995d';
 case 15
    TemplateName = 'SN1996cb';
 case 16
    TemplateName = 'SN1998bp';
 case 17
    TemplateName = 'SN1998bw';
 case 18;
    TemplateName = 'SN1998de';
 case 19
    TemplateName = 'SN1998dt';
 case 20
    TemplateName = 'SN1998es';
 case 21
    TemplateName = 'SN1998s';
 case 22
    TemplateName = 'SN1999ee';
 case 23
    TemplateName = 'SN1999em';
 case 24
    TemplateName = 'SN2000cx';
 case 25
    TemplateName = 'SN2001x';
 case 26
    TemplateName = 'SN2002ap';
 case 27
    TemplateName = 'SN1996L';
 case 28
    TemplateName = 'SN1997cy';
 case 29
    TemplateName = 'SN1998S';
 case 30
    TemplateName = 'SN1999ex';
 otherweise
    error('Unknown SN Template index');
end



%-------------------------------------------
% Template Scaling
% Input : - Template
%         - Spec
%         - Shift Method or value
%         - Scale Method or value
%         - Polynomial divsion scaling (override shift+stretch)
% Outpu : - New scaled Template
%-------------------------------------------
function [ScaledTemplate,ScalePar] = scale_template(Template,Spec,ShiftMethod,ScaleMethod,Poly);

SigmaClip     = 2;
ColW = 1;
ColI = 2;
InterpMethod = 'linear';

if (Poly==0),
   % no polynomial
   if (isstr(ShiftMethod)==0),
      % assume both Shift and Scale are numbers!
      Shift = ShiftMethod;
      Scale = ScaleMethod;
   else
      [Shift,Scale] = call_find_shift_scale_spec(Spec,Template,ShiftMethod,ScaleMethod,SigmaClip);
   end
   ScaledTemplate = [Template(:,ColW), (Template(:,ColI)-Shift)./Scale];
   ScalePar       = [Shift, Scale];
else
   % polynomial scaling
   % interpolate Template to Spec

   [Factor] = find_shift_scale_spec(Spec,Template,[],[],[],Poly);
   ScaledTemplate = Template;
   ScaledTemplate(:,ColI) = ScaledTemplate(:,ColI).*Factor(:,2);
   ScalePar       = Factor;

   %NewTempI = interp1(Template(:,ColW),Template(:,ColI),Spec(:,ColW),InterpMethod);
   %Factor   = Spec(:,ColI)./NewTempI;
   %% poly fit Factor
   %In0            = find(Factor>0 & isnan(Factor)==0);
   %Par            = polyfit(Spec(In0,ColW),log10(Factor(In0)),Poly);
   %TemplateFactor = 10.^(polyval(Par,Template(:,ColW)));
   %ScaledTemplate = Template;
   %ScaledTemplate(:,ColI) = ScaledTemplate(:,ColI).*TemplateFactor;
   %ScalePar       = [ScaledTemplate(:,ColW), TemplateFactor];
end


%-------------------------------------------
% call_find_shift_scale_spec
% a driver function for find_shift_scale_spec
%-------------------------------------------
function [Shift,Scale] = call_find_shift_scale_spec(Spec,Template,ShiftMethod,ScaleMethod,SigmaClip);

Fit_ScaleOption = 8;
Fit_ShiftOption = 7;

switch ShiftMethod
    case 'fit'
        ScaleMethod = 'fit';
    otherwise
        % do nothing
end

[Shift,Scale] = find_shift_scale_spec(Spec,Template,ShiftMethod,ScaleMethod,SigmaClip);



%--------------------------------------------------------
% redisplay_spec
% set spectra data into UserData in redisplay spectra
%--------------------------------------------------------
function redisplay_spec(ReDisplayOption);

if (nargin==0),
   ReDisplayOption = 'Value';
end     

ColW = 1;
ColI = 2;

Vel  = [];

%--- Read Redshift ---
Redshift    = str2num(get(findobj(gcbf,'Tag','edit_Redshift'),'String'));

%--- Read Ebv ---
Ebv_Source    = str2num(get(findobj(gcbf,'Tag','edit_SourceE_BV'),'String'));
Ebv_Observer  = str2num(get(findobj(gcbf,'Tag','edit_ObserverE_BV'),'String'));

%--- Read R ---
R_Source    = str2num(get(findobj(gcbf,'Tag','edit_SourceR'),'String'));
R_Observer  = str2num(get(findobj(gcbf,'Tag','edit_ObserverR'),'String'));


%--- Read Redshift Template | Spec ---
RedshiftTemplate = get(findobj(gcbf,'Tag','radiobutton_RedshiftTemplate'),'Value');


%--- Get UserData ---
UserData = get(matchspec,'UserData');  
CurTemp  = UserData.CurTemp;
Template = UserData.Template{CurTemp};
Htemp    = UserData.Htemp{CurTemp};
Spec     = UserData.Spec;
Hspec    = UserData.Hspec;
CurTemp  = UserData.CurTemp;


%--- Scale template ---
if (isempty(Htemp)),
   % do nothing
   switch RedshiftTemplate
    case 1
       % No template - do nothing
    case 0
       Spec      = [Spec(:,ColW).*(1+Redshift), Spec(:,ColI)];
    otherwise
       error('Unknonw RedshiftTemplate Option');
   end    
else
   switch RedshiftTemplate
    case 1
       if (isempty(Htemp)),
          % do nothing
       else
          Template = get_spectra(Template,[Ebv_Source, Ebv_Observer],[R_Source, R_Observer],Redshift,Vel);
       end
    case 0
       if (isempty(Htemp)),
          % do nothing
       else
          Template = get_spectra(Template,[Ebv_Source, Ebv_Observer],[R_Source, R_Observer],0,[]);
       end

       Spec      = [Spec(:,ColW).*(1+Redshift), Spec(:,ColI)];
    otherwise
       error('Unknonw RedshiftTemplate Option');
   end

   
   %--- find Shift and Scale ---
   
   switch ReDisplayOption
    case 'Method'
       %--- Read ShiftMethod ---
       ShiftOption = get(findobj(gcbf,'Tag','popupmenu_ZeroShift'),'Value');
       ShiftMethod = index2shiftmethod(ShiftOption);

       %--- Read ScaleMethod ---
       ScaleOption = get(findobj(gcbf,'Tag','popupmenu_Scaling'),'Value');
       ScaleMethod = index2scalemethod(ScaleOption);

       SigmaClip            = 2.5;
       [ShiftVal, ScaleVal] = call_find_shift_scale_spec(Spec,Template,ShiftMethod,ScaleMethod,SigmaClip);
   
       %--- Update Shift and Scale edit windows ---
       set(findobj(gcbf,'Tag','edit_ShiftVal'),'String',num2str(ShiftVal));
       set(findobj(gcbf,'Tag','edit_ScaleVal'),'String',num2str(ScaleVal));
       Poly = 0;   % override polynomial
       
    case 'Value'
       %--- Read Shift ---
       ShiftVal    = str2num(get(findobj(gcbf,'Tag','edit_ShiftVal'),'String'));
       %--- Read Scale ---
       ScaleVal    = str2num(get(findobj(gcbf,'Tag','edit_ScaleVal'),'String'));
       %--- Read Poly ---
       Poly        = str2num(get(findobj(gcbf,'Tag','edit_PolyDivide'),'String'));
      
    otherwise
       error('Unknown ReDisplayOption');
   end
   
   %--- find Shift and Scale ---
   [Template,ScalePar] = scale_template(Template,Spec,ShiftVal,ScaleVal,Poly);
end




%--- Get template line properties ---
TemplateColor = get(findobj(gcbf,'Tag','pushbutton_Color'),'ForegroundColor');
TemplateLine  = get(findobj(gcbf,'Tag','popupmenu_LineType'),'Value');
TemplateLine = index2linetype(TemplateLine);

switch RedshiftTemplate
 case 1
    if (isempty(Htemp)),
       % do nothing
    else
       set(get(Htemp,'Children'),'XData',Template(:,ColW));
       set(get(Htemp,'Children'),'YData',Template(:,ColI));

       set(get(Htemp,'Children'),'Color',TemplateColor);
       set(get(Htemp,'Children'),'LineStyle',TemplateLine);
    end
 case 0
    set(get(Hspec,'Children'),'XData',Spec(:,ColW));
    set(get(Hspec,'Children'),'YData',Spec(:,ColI));
    if (isempty(Htemp)),
       % do nothing
    else
       set(get(Htemp,'Children'),'XData',Template(:,ColW));
       set(get(Htemp,'Children'),'YData',Template(:,ColI));
    end

 otherwise
    error('Unknonw RedshiftTemplate Option');
end














% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function pushbutton_Color_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_Color (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton_AccretionDisk.
function pushbutton_AccretionDisk_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_AccretionDisk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ColW = 1;
ColI = 2;
WaveVec = [1000:10:20000].';

%--- Read Accretion disk parameters ---
Mass     = str2num(get(findobj(gcbf,'Tag','edit_AD_Mass'),'String'));
MassRate = str2num(get(findobj(gcbf,'Tag','edit_AD_SM_Year'),'String'));
InnerRad = str2num(get(findobj(gcbf,'Tag','edit_InnerRad'),'String'));
OuterRad = str2num(get(findobj(gcbf,'Tag','edit_OuterRad'),'String'));

Template  = accretion_disk(Mass,MassRate,InnerRad,OuterRad,WaveVec);

UserData       = get(matchspec,'UserData');

UserData.Ntemp = UserData.Ntemp + 1;
UserData.CurTemp = UserData.Ntemp;
UserData.Template{UserData.Ntemp} = Template;

   
%---------------------
%--- Plot Template ---
%---------------------
TemplateColor = get(findobj(gcbf,'Tag','pushbutton_Color'),'ForegroundColor');


Hfig = get(get(UserData.Hspec,'Parent'),'Parent');
figure(Hfig);
Htemp         = stairs(Template(:,ColW),Template(:,ColI));

set(Htemp,'Color',TemplateColor);

UserData.Htemp{UserData.Ntemp} = Htemp;
set(matchspec,'UserData',UserData);





function edit_AD_Mass_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AD_Mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AD_Mass as text
%        str2double(get(hObject,'String')) returns contents of edit_AD_Mass as a double


% --- Executes during object creation, after setting all properties.
function edit_AD_Mass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AD_Mass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_AD_SM_Year_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AD_SM_Year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_AD_SM_Year as text
%        str2double(get(hObject,'String')) returns contents of edit_AD_SM_Year as a double


% --- Executes during object creation, after setting all properties.
function edit_AD_SM_Year_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AD_SM_Year (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_InnerRad_Callback(hObject, eventdata, handles)
% hObject    handle to edit_InnerRad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_InnerRad as text
%        str2double(get(hObject,'String')) returns contents of edit_InnerRad as a double


% --- Executes during object creation, after setting all properties.
function edit_InnerRad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_InnerRad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_OuterRad_Callback(hObject, eventdata, handles)
% hObject    handle to edit_OuterRad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_OuterRad as text
%        str2double(get(hObject,'String')) returns contents of edit_OuterRad as a double


% --- Executes during object creation, after setting all properties.
function edit_OuterRad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_OuterRad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton_RedshiftSpectrum.
function radiobutton_RedshiftSpectrum_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_RedshiftSpectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_RedshiftSpectrum



function edit_CCF_Step_Callback(hObject, eventdata, handles)
% hObject    handle to edit_CCF_Step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_CCF_Step as text
%        str2double(get(hObject,'String')) returns contents of edit_CCF_Step as a double


% --- Executes during object creation, after setting all properties.
function edit_CCF_Step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CCF_Step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_LoadTemplateFromFile.
function pushbutton_LoadTemplateFromFile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_LoadTemplateFromFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ColW = 1;
ColI = 2;

%--- Col X ---
TempColW    = str2num(get(findobj(gcbf,'Tag','edit_TempColW'),'String'));
%--- Col Y ---
TempColI    = str2num(get(findobj(gcbf,'Tag','edit_TempColI'),'String'));
%--- Col Err ---
TempColErr  = str2num(get(findobj(gcbf,'Tag','edit_TempColErr'),'String'));

[FileName,FilePath] = uigetfile;
Temp                = load([FilePath,FileName]);
if (size(Temp,2)>2),
   Template            = Temp(:,[TempColW, TempColI, TempColErr]);
else
   Template            = Temp(:,[TempColW, TempColI]);
end
   
UserData       = get(matchspec,'UserData');

UserData.Ntemp = UserData.Ntemp + 1;
UserData.CurTemp = UserData.Ntemp;
UserData.Template{UserData.Ntemp} = Template;

%---------------------
%--- Plot Template ---
%---------------------
TemplateColor = get(findobj(gcbf,'Tag','pushbutton_Color'),'ForegroundColor');


Hfig = get(get(UserData.Hspec,'Parent'),'Parent');
figure(Hfig);
Htemp         = stairs(Template(:,ColW),Template(:,ColI));

set(Htemp,'Color',TemplateColor);

UserData.Htemp{UserData.Ntemp} = Htemp;
set(matchspec,'UserData',UserData);





function edit_TempColW_Callback(hObject, eventdata, handles)
% hObject    handle to edit_TempColW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_TempColW as text
%        str2double(get(hObject,'String')) returns contents of edit_TempColW as a double


% --- Executes during object creation, after setting all properties.
function edit_TempColW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_TempColW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_TempColI_Callback(hObject, eventdata, handles)
% hObject    handle to edit_TempColI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_TempColI as text
%        str2double(get(hObject,'String')) returns contents of edit_TempColI as a double


% --- Executes during object creation, after setting all properties.
function edit_TempColI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_TempColI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_TempColErr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_TempColErr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_TempColErr as text
%        str2double(get(hObject,'String')) returns contents of edit_TempColErr as a double


% --- Executes during object creation, after setting all properties.
function edit_TempColErr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_TempColErr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_LoadSpecFromFile.
function pushbutton_LoadSpecFromFile_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_LoadSpecFromFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ColW = 1;
ColI = 2;

%--- Col X ---
SpecColW    = str2num(get(findobj(gcbf,'Tag','edit_SpecColW'),'String'));
%--- Col Y ---
SpecColI    = str2num(get(findobj(gcbf,'Tag','edit_SpecColI'),'String'));
%--- Col Err ---
SpecColErr  = str2num(get(findobj(gcbf,'Tag','edit_SpecColErr'),'String'));

[FileName,FilePath] = uigetfile;
Sp                  = load([FilePath,FileName]);
if (size(Sp,2)>2),
   Spec            = Sp(:,[SpecColW, SpecColI, SpecColErr]);
else
   Spec            = Sp(:,[SpecColW, SpecColI]);
end
   
UserData       = get(matchspec,'UserData');


UserData.Spec = Spec;

%---------------------
%--- Plot Template ---
%---------------------


Hfig = get(get(UserData.Hspec,'Parent'),'Parent');
figure(Hfig);
Hspec         = stairs(Spec(:,ColW),Spec(:,ColI));


UserData.Hspec = Hspec;
set(matchspec,'UserData',UserData);





function edit_SpecColW_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SpecColW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SpecColW as text
%        str2double(get(hObject,'String')) returns contents of edit_SpecColW as a double


% --- Executes during object creation, after setting all properties.
function edit_SpecColW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SpecColW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SpecColY_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SpecColY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SpecColY as text
%        str2double(get(hObject,'String')) returns contents of edit_SpecColY as a double


% --- Executes during object creation, after setting all properties.
function edit_SpecColY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SpecColY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_SpecColErr_Callback(hObject, eventdata, handles)
% hObject    handle to edit_SpecColErr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_SpecColErr as text
%        str2double(get(hObject,'String')) returns contents of edit_SpecColErr as a double


% --- Executes during object creation, after setting all properties.
function edit_SpecColErr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_SpecColErr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_MeasureLines.
function pushbutton_MeasureLines_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_MeasureLines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(matchspec,'Visible','off');
%--- Read Measure Lines Template | Spec ---
ML_Template = get(findobj(gcbf,'Tag','radiobutton_ML_Template'),'Value');

UserData = get(matchspec,'UserData');  
CurTemp  = UserData.CurTemp;
Template = UserData.Template{CurTemp};
Htemp    = UserData.Htemp{CurTemp};
Spec     = UserData.Spec;
Hspec    = UserData.Hspec;
CurTemp  = UserData.CurTemp;


switch ML_Template
    case 1
        % Measure lines in template
        SpecM = Template;
        SpecType = 'Template';
    case 0
        % measure lines in spectrum
        SpecM = Spec;
        SpecType = 'Spectrum';
    otherwise
        error('Unknown ML_Template Option');
end

disp(sprintf('--- Measure Lines in rest frame %s ---',SpecType))
disp(sprintf('--- OPen New figure ---'))
[LinesProp,ConPos]=eq_width(SpecM);
%            [Wavelength at extramum intensity,
%             Line center (first moment),
%             Line second moment,
%             Line FWHM, (not implemented)
%             Extramum value,
%             Flux of line,
%             Line eqwivalent width].


disp(sprintf('Line measurments results are available in session variables: '))
disp(sprintf('     matchspec_LinesProp'))
disp(sprintf('     matchspec_ConPos'))
disp(sprintf('Check out eq_width.m for more details'))

assignin('base','matchspec_LinesProp',LinesProp);
assignin('base','matchspec_ConPos',ConPos);
   




% --- Executes during object creation, after setting all properties.
function radiobutton_RedshiftTemplate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to radiobutton_RedshiftTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


