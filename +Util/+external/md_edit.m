function varargout = md_edit(varargin)
%MD_EDIT creates an edit window for interactive editing of variables
%      This functions allows you to edit variables in a spreadsheet like way.
%      To edit a variable 'in1', type:
%
%        out1 = md_edit(in1);
%
%      The variable 'out1' will contain the edited data. You can edit several
%      matrices of equal dimensions simultaneously by typing them as the first
%      parameters of MD_EDIT:
%
%        [out1,out2,out3,out4] = md_edit(in1, in2, in3, in4);
%
%      The elements of the matrices at specific indices can be edited in the
%      upper part of the window. The indices can be selected in the lower part.
%      If one matrix is edited, the editing takes place in the lower part.
%
%      When the input variable is a char matrix containing only text,
%      you can choose to edit it as a matrix of characters (default),
%      a set of single lines, or as a text (like a small text editor).
%
%      When the input variable is a cell array clicking on an element switches
%      the editing to the editing of that element. The popupmenu at the bottom
%      of the figure can be used to get back to a higher level. Similar for a
%      struct array of which the fields can be selected in the upper part of
%      the window.
%
%      When the array has dimension higher than two you are asked to select
%      a slice out of the array to display. This is done by selecting two
%      dimensions that you want to display completely and by specifying
%      indices for the other dimensions.
%
%      The array size can be changed by selecting columns or rows and deleting
%      and inserting them.
%
%      OPTIONS:
%
%        out1 = md_edit(in1, option1, option2, ...);
%
%      Options can be specified after the input arrays. Valid options are:
%
%        'specmode' followed by 'array','line','multiline' to start editing
%                   a char matrix in character, line, or text mode.
%
%        'dimfixed' followed by a vector containing 1 for dimensions that
%                   should have fixed size, e.g. 'dimfixed',[1 0] allows the
%                   insertion and deletion of columns, but not of rows.
%
%        'label'    followed by a cell array containing for each dimension
%                   one of the following entries:
%
%                   '123'             array indices numbered 1,2,3, etc. 
%                   'ABC'             array indices numbered A,B,C, etc. 
%                   'none'            blank buttons
%                   {'string1','string2', ...} user defined names
%
%         'dimlabel' followed by a cell array containing for each dimension
%                    a string indicating the name of that dimension.
%
%         'multiple' followed by a cell array containing for each input
%                    array a name to display as fieldname.
%
%      For instance:
%         xt = get(gca,'xtick')';
%         xtl = get(gca,'xticklabel');
%         [xt,xtl] = md_edit(xt,xtl,'multiple',{'tick','ticklabel'}, ...
%                             'label',{'123' 'none'},'dimlabel',{'tick number',''});
%         set(gca,'xtick',xt,'xticklabel',xtl);

%      Copyright (c) Aug 19, 1998 by H.R.A. Jagers
%                                    The Netherlands

%*-------------------------------------------------------------------------------------------------
%* CATCH STACKUDF CALLS
%*-------------------------------------------------------------------------------------------------
if (nargin==2) && strcmp(varargin{1},'MDEDIT-INTERNAL') && ...
      ~isempty(gcbf) && strcmp(get(gcbf,'tag'),'IDEAS - GUI'),
  LcStackUDF(gcbf,'CommandStack',varargin{2});
  return;
end;

%*-------------------------------------------------------------------------------------------------
%* DEFINE WINDOW SIZES
%*-------------------------------------------------------------------------------------------------

DisplayNumber.Fields=1;
DisplayNumber.Columns=14;
DisplayNumber.Rows=20;

%*-------------------------------------------------------------------------------------------------
%* USER INPUT SHOULD BE EVALUATED IN THE CALLING M-FILE (COULD BE SET TO BASE)
%*-------------------------------------------------------------------------------------------------

EvalWhere='caller';

%*-------------------------------------------------------------------------------------------------
%* PROCESS INPUT: CHECK FOR EDITING PARALLEL EDITING OF TWO OR MORE MATRICES
%* SINCE OPTIONS MIGHT BE SPECIFIED; THE NUMBER OF MATRICES CANNOT BE DETERMINED
%* ON THE NUMBER OF INPUT ARGUMENTS. THEREFORE IT HAS TO BE DETERMINED ON THE
%* NUMBER OF OUTPUT ARGUMENTS.
%*-------------------------------------------------------------------------------------------------

multiedit=max(nargout,1);

%*-------------------------------------------------------------------------------------------------
%* PROCESS INPUT: CREATE APPROPRIATE EMPTY OUTPUT WHILE EDITING IS NOT FINISHED
%*-------------------------------------------------------------------------------------------------

if nargout>0,
  varargout=cell(1,nargout);
  for i=1:nargout,
    varargout{i}=[];
  end;
end;
%*-------------------------------------------------------------------------------------------------
%* PROCESS INPUT: CHECK NUMBER OF INPUT AND OUTPUT ARGUMENTS
%*-------------------------------------------------------------------------------------------------

if (nargin==0),
  uiwait(msgbox('At least one input argument expected.'));
  return;
elseif (nargout>nargin),
  uiwait(msgbox('Unexpected number of output arguments.'));
  return;
end;

%*-------------------------------------------------------------------------------------------------
%* PROCESS INPUT: STORE DATA IN OriginalData AND OPTIONS IN options
%* ADD 'multiple' KEYWORD WHEN MORE THAN ONE MATRIX WAS
%* SPECIFIED.
%*-------------------------------------------------------------------------------------------------

if multiedit>1,
  OriginalData=varargin(1:multiedit);
  options={'multiple','123',varargin{(multiedit+1):nargin}};
else
  OriginalData=varargin{1};
  options=varargin((multiedit+1):nargin);
end;

%*-------------------------------------------------------------------------------------------------
%* CREATE THE EDITING WINDOW
%*-------------------------------------------------------------------------------------------------

[Handles,Color]=LcUiEdit(DisplayNumber);

%*-------------------------------------------------------------------------------------------------
%* MAKE A COPY OF THE MATRIX FOR EDITING AND KEEP THE ORIGINAL FOR RESETTING
%*-------------------------------------------------------------------------------------------------

EditedData=OriginalData;

%*-------------------------------------------------------------------------------------------------
%* INITIALIZE EDITING VARIABLES
%*-------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%
Selected.Columns=[];      % Contains the indices of the selected columns.
Selected.Rows=[];         % Contains the indices of the selected rows.
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Display.Dimensions=[];    % Contains the dimension numbers currently displayed
Display.TotalDataSize=[]; % Contains the size of the total data structure at the level of editing
                          % Display.TotalDataSize(Display.Dimensions) = size(Display.Data) 
Display.Data=[];          % Contains the displayed data and all data below
Display.Element=[1 1];    % Contains the coordinates of the displayed structure (when relevant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Display.Index.Ref.type='';                %
Display.Index.Ref.subs={};                %
Display.Index.Ref=Display.Index.Ref(1:0); % Display.Index contains indices to obtain displayed data
Display.PreviousSelectedDims={};          % Contains the subs part of the last Display.Index before redrawing
Display.MultipleField={};                 % Display.Field contains in the first column the index of
                                          % the dataset in multiple mode and in the second column
                                          % field name if it is a struct dataset.
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Offset.Row=0;             % Offset.Row is upper row number minus one;
                          % so rows Offset.Row+1 to Offset.Row+DisplayNumber.Rows are shown (if available)
Offset.Column=0;          % Offset.Column left column number minus one,
                          % so columns Offset.Column+1 to Offset.Column+show_numbercolumns are shown (if available)
Offset.Field=0;           % Offset of the element list of a structure (if relevant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Options.Multiple={};      % Labels of matrices in case of multiple matrices
Options.MultipleLines=0;  % One if value vectors and string vectors
Options.DimElmLabel={};   % Labels of dimensions
Options.DimFixed=[];      % Array of logical values determining whether certain dimensions are fixed i.e. size(*,n)==constant
Options.DimLabel={};      % Labels of indices
Options.EditModes={'array','line','multiline'};
Options.EditMode=Options.EditModes{1};    % Initial EditMode
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%*-------------------------------------------------------------------------------------------------
%* START OF EDITING LOOP 
%*-------------------------------------------------------------------------------------------------

gui_quit=0;               % Becomes one if the interface has to quit.
stack=[];                 % Contains the stack of commands; read from 'userdata' field of the figure

while ~gui_quit,

%%*************************************************************************************************
%%%% UPDATE SCREEN BEFORE WAITING FOR COMMAND
%%*************************************************************************************************

  drawnow;

%%*************************************************************************************************
%%%% WAIT UNTIL A COMMAND IS ON THE STACK IN THE USERDATA FIELD OF THE FIGURE
%%*************************************************************************************************

  if ishandle(Handles.Figure),
    if isempty(LcGetUDF(Handles.Figure,'CommandStack')),
      LcWaitForUDF(Handles.Figure,'CommandStack');
    end;
  end;

%%*************************************************************************************************
%%%% SET POINTER TO WATCH WHILE PROCESSING COMMANDS ON STACK
%%%% FIRST CHECK WHETHER FIGURE STILL EXISTS
%%*************************************************************************************************

  if ishandle(Handles.Figure),
    stack=LcGetUDF(Handles.Figure,'CommandStack');
    LcSetUDF(Handles.Figure,'CommandStack',{});
    set(Handles.Figure,'pointer','watch');
  else
    uiwait(msgbox('Unexpected removal of Edit window!','modal'));
    gui_quit=1;
  end;

%%*************************************************************************************************
%%%% START OF WHILE COMMANDS ON STACK LOOP
%%*************************************************************************************************

  while ~isempty(stack),
    cmd=stack{1};
    stack=stack(2:size(stack,1),:);

    if isnan(cmd(1)), % rebuild window (initialize or reset if isnan(cmd(2)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% UPDATE DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     [Selected,Display,Offset,Options]= ...
       LcUiEdit_update(DisplayNumber,Color,cmd,Handles,EditedData,Selected,Display,Offset,Options,options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% UPDATE DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif min(cmd)>0, % index element edited (or clicked)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ELEMENT CHANGED OR CLICKED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Determine type of data
      if ~isempty(Options.Multiple),
        DataClass='struct';
      else
        DataClass=class(Display.Data);
      end

      % Determine what to do given data type
      switch DataClass,
      case {'double','uint8','sparse'},
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% A NUMERIC ELEMENT WAS EDITED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(Handles.Index.Elements(cmd(1),cmd(2)),'foregroundcolor',Color.Message);
        lasterr('');
        %temp=evalin(EvalWhere,get(Handles.Index.Elements(cmd(1),cmd(2)),'string'),'NaN');
        temp=str2num(get(Handles.Index.Elements(cmd(1),cmd(2)),'string'));
        Str=lasterr;
        if isempty(Str),
          if isequal(size(temp),[1 1]) && isnumeric(temp),
            Display.Data(Offset.Row+cmd(1),Offset.Column+cmd(2))=temp;
            EditedData=LcSubsAsgn(EditedData,Display.Index(end).Ref,Display.Data);
          else
            Str='Input must be a single numeric value.';
            uiwait(msgbox(Str,'modal'));
            set(Handles.General.MessageWindow,'string',Str);
          end;
        else
          Str=multiline(Str);
          Str=Str(size(Str,1),:);
          uiwait(msgbox(Str,'modal'));
          set(Handles.General.MessageWindow,'string',Str);
        end;
        set( Handles.Index.Elements(cmd(1),cmd(2)) ,'string', sprintf('%g',double(Display.Data(Offset.Row+cmd(1),Offset.Column+cmd(2)))), ...
                                 'foregroundcolor',Color.Normal );

      case {'char'},
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% A STRING ELEMENT WAS EDITED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(Options.EditMode,Options.EditModes{1}), % array
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%% A CHARACTER WAS EDITED
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          set(Handles.Index.Elements(cmd(1),cmd(2)),'foregroundcolor',Color.Message);
          temp=get(Handles.Index.Elements(cmd(1),cmd(2)),'string');
          if length(temp)>1 && temp(1)=='#', % possibly ASCII code intended
            if temp(1)=='#' && all(abs(temp(2:length(temp)))>=48) & all(abs(temp(2:length(temp)))<=57),
              temp=sscanf(temp(2:length(temp)),'%i');
              if temp>65535,
                Str = 'Valid character codes are 0:65535.';
                uiwait(msgbox(Str,'modal'));
                set(Handles.General.MessageWindow,'string',Str);
              else
                Display.Data(Offset.Row+cmd(1),Offset.Column+cmd(2))=char(temp);
                EditedData=LcSubsAsgn(EditedData,Display.Index(end).Ref,Display.Data);
              end;
            else
              Str = 'Input must be a single character or # followed by a character code.';
              uiwait(msgbox(Str,'modal'));
              set(Handles.General.MessageWindow,'string',Str);
            end;
          elseif length(temp)==1,  
            Display.Data(Offset.Row+cmd(1),Offset.Column+cmd(2))=temp;
            EditedData=LcSubsAsgn(EditedData,Display.Index(end).Ref,Display.Data);
          elseif isempty(temp),
            temp=char(0);
            Display.Data(Offset.Row+cmd(1),Offset.Column+cmd(2))=char(0);
            EditedData=LcSubsAsgn(EditedData,Display.Index(end).Ref,Display.Data);
          else
            Str = 'Input must be a single character or # followed by a character code.';
            uiwait(msgbox(Str,'modal'));
            set(Handles.General.MessageWindow,'string',Str);
          end;
          if abs(Display.Data(Offset.Row+cmd(1),Offset.Column+cmd(2)))<32 || ...
             abs(Display.Data(Offset.Row+cmd(1),Offset.Column+cmd(2)))>255,
            Str=['#' num2str(abs(Display.Data(Offset.Row+cmd(1),Offset.Column+cmd(2))))];
            set(Handles.Menu.Char.Main,'enable','off');
          else
            Str=Display.Data(Offset.Row+cmd(1),Offset.Column+cmd(2));
          end;
          set( Handles.Index.Elements(cmd(1),cmd(2)) ,'string', Str, ...
                                   'foregroundcolor',Color.Normal );

        elseif strcmp(Options.EditMode,Options.EditModes{2}), % line
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%% A LINE WAS EDITED
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          temp=deblank(get(Handles.LineEdit.Lines(cmd(1)),'string'));
          if length(temp)>size(Display.Data,2),
            % Fit matrix to new line
            Display.Data(:,(size(Display.Data,2)+1):length(temp))=' ';
            Display.Data(Offset.Row+cmd(1),:)=temp;
          elseif length(temp)==size(Display.Data,2),
            % No fitting necessary
            if ~isempty(Display.Data),
              Display.Data(Offset.Row+cmd(1),:)=temp;
            end;
          else
            % Fit new text to matrix
            temp((length(temp)+1):size(Display.Data,2))=' ';
            Display.Data(Offset.Row+cmd(1),:)=temp;
          end;
          EditedData=LcSubsAsgn(EditedData,Display.Index(end).Ref(1:(end-1)),Display.Data);
          Display.TotalDataSize=size(LcSubsRef(EditedData,Display.Index(end).Ref(1:(end-1))));

        elseif strcmp(Options.EditMode,Options.EditModes{3}), % multiline
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%% A TEXT MATRIX WAS EDITED
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          Display.Data=deblank(get(Handles.TextEdit.Text,'string'));
          EditedData=LcSubsAsgn(EditedData,Display.Index(end).Ref(1:(end-1)),Display.Data);
          Display.TotalDataSize=size(LcSubsRef(EditedData,Display.Index(end).Ref(1:(end-1))));
        end;

      case {'cell'},
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% A CELL ELEMENT WAS CLICKED ON
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch get(Handles.Figure,'selectiontype'),
        case 'normal',
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%% A CELL ELEMENT WAS LEFT CLICKED ON
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          TempVar=Display.Index(end).Ref(end).subs;
          TempVar(Display.Dimensions)={Offset.Row+cmd(1) Offset.Column+cmd(2)};
          Display.Index(end+1)=Display.Index(end);
          Display.Index(end).Ref(end).type='{}';
          Display.Index(end).Ref(end).subs=TempVar;
          if isstruct(Display.Data{Offset.Row+cmd(1),Offset.Column+cmd(2)}),
            Display.Element=[1 1];
          end;
          LcStackUDF(Handles.Figure,'CommandStack',[NaN 1],'top');
        case 'alt',
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%% A CELL ELEMENT WAS RIGHT CLICKED ON
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          TempRef=Display.Index(end).Ref;
          TempRef(end).type='{}';
          TempRef(end).subs(Display.Dimensions)={Offset.Row+cmd(1) Offset.Column+cmd(2)};
          TempVar=LcSubsRef(EditedData,TempRef);
          switch class(TempVar),
          case {'double','uint8','char','sparse','cell','struct'},
            Str = 'Cell does not contain object.';
            uiwait(msgbox(Str,'modal'));
            set(Handles.General.MessageWindow,'string',Str);
          otherwise,
            Confirm=questdlg('Convert object to struct?', ...
                       'Conversion warning','Yes','No','No');
            if strcmp(Confirm,'Yes')
              EditedData=LcSubsAsgn(EditedData,TempRef,struct(TempVar));
              LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
            end
          end;
        end;
      case {'struct'}, % struct -> clicked on
                       % struct data or ~isempty(Options.Multiple),
        Display.Element=[Offset.Row+cmd(1) Offset.Column+cmd(2)];
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
      end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ELEMENT CHANGED OR CLICKED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif min(cmd)==0, % column or row selected

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% A ROW OR COLUMN WAS (UN)SELECTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if cmd(1)==0,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% A COLUMN WAS (UN)SELECTED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cmd(2)=cmd(2)+Offset.Column;
        TempVar=find( ~(Selected.Columns-cmd(2)) );
        if any(TempVar),
          Selected.Columns=Selected.Columns(logical((Selected.Columns-cmd(2))));
        else
          Selected.Columns=sort([Selected.Columns cmd(2)]);
        end
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');

      else % cmd(2)==0,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% A ROW WAS (UN)SELECTED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cmd(1)=cmd(1)+Offset.Row;
        TempVar=find( ~(Selected.Rows-cmd(1)) );
        if any(TempVar),
          Selected.Rows=Selected.Rows(logical((Selected.Rows-cmd(1))));
        else
          Selected.Rows=sort([Selected.Rows cmd(1)]);
        end;
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
      end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% A ROW OR COLUMN WAS (UN)SELECTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif cmd(1)==-1, % other button

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ANOTHER UICONTROL OBJECT WAS USED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if cmd(2)==1,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% THE ACCEPT BUTTON WAS PRESSED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gui_quit=1;

      elseif cmd(2)==2,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% THE RESET BUTTON WAS PRESSED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        EditedData=OriginalData;
        Display.Index(2:end)=[];
        Display.Index.Ref=Display.Index.Ref(1:0);
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 1],'top');

      elseif cmd(2)==3,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% THE COLUMN SLIDER WAS MOVED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TempVar=Offset.Column;
        Offset.Column=get(Handles.Index.Scrollbar.Column,'value');
        if (Offset.Column>TempVar),
          Offset.Column=max(round(Offset.Column),min(TempVar+1,get(Handles.Index.Scrollbar.Column,'max')));
        elseif(Offset.Column<TempVar),
          Offset.Column=min(round(Offset.Column),max(TempVar-1,0));
        end;
        if (Offset.Column~=TempVar),
          LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
        end;

      elseif cmd(2)==4,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% THE ROW SLIDER WAS MOVED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TempVar=Offset.Row;
        switch Options.EditMode,
        case Options.EditModes(1), % index editing
          Offset.Row=-get(Handles.Index.Scrollbar.Row,'value');
        case Options.EditModes(2), % line editing
          Offset.Row=-get(Handles.LineEdit.Scrollbar.Line,'value');
        end;
        if (Offset.Row>TempVar),
          switch Options.EditMode,
          case Options.EditModes(1), % index editing
            Offset.Row=max(round(Offset.Row),min(TempVar+1,-get(Handles.Index.Scrollbar.Row,'min')));
          case Options.EditModes(2), % line editing
            Offset.Row=max(round(Offset.Row),min(TempVar+1,-get(Handles.LineEdit.Scrollbar.Line,'min')));
          end;
        elseif(Offset.Row<TempVar),
          Offset.Row=min(round(Offset.Row),max(TempVar-1,0));
        end;
        if (Offset.Row~=TempVar),
          LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
        end;

      elseif cmd(2)==5,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% THE LIST POPUPMENU WAS USED: CHANGE BACK TO HIGHER EDITING LEVEL ?
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TempVar=get(Handles.General.List,'value');
        Str=get(Handles.General.List,'string');
        Display.Index=Display.Index(1:TempVar);
        if length(Display.Index(end).Ref(end).subs)>2, % multidimensional
          Display.PreviousSelectedDims=Display.Index(end).Ref(end).subs;
        end;
        Display.Index(end).Ref(end)=[];
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 1],'top');

      elseif cmd(2)==6,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% THE FIELD SLIDER OF THE STRUCT WINDOW PART WAS MOVED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TempVar=Offset.Field;
        Offset.Field=-get(Handles.Field.Scrollbar.Field,'value');
        if (Offset.Field>TempVar),
          Offset.Field=max(round(Offset.Field),min(TempVar+1,-get(Handles.Field.Scrollbar.Field,'min')));
        elseif(Offset.Field<TempVar),
          Offset.Field=min(round(Offset.Field),max(TempVar-1,0));
        end;
        if (Offset.Field~=TempVar),
          LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
        end;

      else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% AN UNKNOWN UICONTROL (BUTTON/SLIDER/POPUPMENU) WAS USED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Str='Unknown uicontrol command';
        fprintf(1,[Str ':\n']);
        disp(cmd);
        set(Handles.General.MessageWindow,'string',[Str '.']);
      end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ANOTHER UICONTROL OBJECT WAS USED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif cmd(1)==-2,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% THE COLUMN MENU WAS USED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if cmd(2)==1,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% THE DELETE COLUMN(S) OPTION WAS SELECTED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dataset=1;
        NDatasets=length(Display.Data);
        First=1;
        while (~isempty(Options.Multiple) && (Dataset<=NDatasets)) || First,
          First=0;
          TempSubsIndex=Display.Index(end).Ref;
          if ~isempty(Options.Multiple),
            TempSubsIndex(end+1)=TempSubsIndex(end);
            TempSubsIndex(end-1).type='{}';
            TempSubsIndex(end-1).subs={Dataset};
            Dataset=Dataset+1;
          end;
          TempSubsIndex(end).subs(:)={':'};
          TempSubsIndex(end).subs{Display.Dimensions(2)}=Selected.Columns;

          TempVar=LcSubsRef(EditedData,TempSubsIndex(1:(end-1)));
          switch class(TempVar),
          case {'struct'},
            % To correct bugs in Matlab 5.2
            % A(5,1).a=1; A(:,1)=[]; returns a 1x0 struct
            % B=ones(5,1); B(:,1)=[]; returns a 5x0 matrix
            % A a 5x0 struct with field a; A(2,:)=[]; returns: ??? Improper null index in null matrix.
            % B=ones(5,0); B(2,:)=[]; returns a 4x0 matrix

            if ~isempty(TempVar),
              TempVar2=size(LcSubsRef(EditedData,TempSubsIndex(1:(end-1))));
              TempVar2(Display.Dimensions(2))=TempVar2(Display.Dimensions(2))-length(Selected.Columns);
              EditedData=subsasgn(EditedData,TempSubsIndex,[]);
              TempVar3=size(LcSubsRef(EditedData,TempSubsIndex(1:(end-1))));
              if ~isequal(TempVar2,TempVar3),
                EditedData=LcSubsAsgn(EditedData,TempSubsIndex(1:(end-1)),reshape(LcSubsRef(EditedData,TempSubsIndex(1:(end-1))),TempVar2));
              end;
            else
              TempVar2=size(LcSubsRef(EditedData,TempSubsIndex(1:(end-1))));
              TempVar2(Display.Dimensions(2))=TempVar2(Display.Dimensions(2))-length(Selected.Columns);
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex(1:(end-1)),reshape(LcSubsRef(EditedData,TempSubsIndex(1:(end-1))),TempVar2));
            end;
          otherwise,
            EditedData=subsasgn(EditedData,TempSubsIndex,[]);
          end;
        end;

        Display.Element(2)=max(1,Display.Element(2)-sum(Selected.Columns<=Display.Element(2)));
        Selected.Columns=[];
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');

      elseif cmd(2)==2, % insert column left
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% THE INSERT COLUMN LEFT OPTION WAS SELECTED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dataset=1;
        NDatasets=length(Display.Data);
        First=1;
        while (~isempty(Options.Multiple) && (Dataset<=NDatasets)) || First,
          First=0;
          TempSubsIndex=Display.Index(end).Ref;
          if ~isempty(Options.Multiple),
            TempSubsIndex(end+1)=TempSubsIndex(end);
            TempSubsIndex(end-1).type='{}';
            TempSubsIndex(end-1).subs={Dataset};
            Dataset=Dataset+1;
          end;
          TempSubsIndex(end).subs(:)={':'};

          TempVar=LcSubsRef(EditedData,TempSubsIndex(1:(end-1)));
          if isempty(TempVar),
            TempVar2=size(TempVar);
            TempVar2(Display.Dimensions(2))=TempVar2(Display.Dimensions(2))+length(Selected.Columns);
            EditedData=LcSubsAsgn(EditedData,TempSubsIndex(1:(end-1)),reshape(TempVar,TempVar2));
          else
            TempSubsIndex(end).subs{Display.Dimensions(2)}=setdiff(1:(size(TempVar,Display.Dimensions(2))+length(Selected.Columns)),Selected.Columns+(1:length(Selected.Columns))-1);
            EditedData=LcSubsAsgn(EditedData,TempSubsIndex,TempVar);
            TempSubsIndex(end).subs{Display.Dimensions(2)}=Selected.Columns+(1:length(Selected.Columns))-1;
            switch class(TempVar),
            case {'char'},
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,' ');
            case {'cell'},
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,{[]});
            case {'double','uint8','sparse'},
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,0);
            case {'struct'},
              TempVar2=[];
              Fields=fieldnames(TempVar);
              for j=1:length(Fields),
                TempVar3.type='.';
                TempVar3.subs=Fields{j};
                TempVar2=LcSubsAsgn(TempVar2,TempVar3,[]);
              end;
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,TempVar2);
            otherwise,
            end;
          end;
        end;

        Display.Element(2)=Display.Element(2)+sum(Selected.Columns<=Display.Element(2));
        Selected.Columns=Selected.Columns+(1:length(Selected.Columns));
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');

      elseif cmd(2)==3,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% THE INSERT COLUMN RIGHT OPTION WAS SELECTED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dataset=1;
        NDatasets=length(Display.Data);
        First=1;
        while (~isempty(Options.Multiple) && (Dataset<=NDatasets)) || First,
          First=0;
          TempSubsIndex=Display.Index(end).Ref;
          if ~isempty(Options.Multiple),
            TempSubsIndex(end+1)=TempSubsIndex(end);
            TempSubsIndex(end-1).type='{}';
            TempSubsIndex(end-1).subs={Dataset};
            Dataset=Dataset+1;
          end;
          TempSubsIndex(end).subs(:)={':'};

          TempVar=LcSubsRef(EditedData,TempSubsIndex(1:(end-1)));
          if isempty(TempVar),
            TempVar2=size(TempVar);
            TempVar2(Display.Dimensions(2))=TempVar2(Display.Dimensions(2))+length(Selected.Columns);
            EditedData=LcSubsAsgn(EditedData,TempSubsIndex(1:(end-1)),reshape(TempVar,TempVar2));
          else
            TempSubsIndex(end).subs{Display.Dimensions(2)}=setdiff(1:(size(TempVar,Display.Dimensions(2))+length(Selected.Columns)),Selected.Columns+(1:length(Selected.Columns)));
            EditedData=LcSubsAsgn(EditedData,TempSubsIndex,TempVar);
            TempSubsIndex(end).subs{Display.Dimensions(2)}=Selected.Columns+(1:length(Selected.Columns));
            switch class(TempVar),
            case {'char'},
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,' ');
            case {'cell'},
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,{[]});
            case {'double','uint8','sparse'},
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,0);
            case {'struct'},
              TempVar2=[];
              Fields=fieldnames(TempVar);
              for j=1:length(Fields),
                TempVar3.type='.';
                TempVar3.subs=Fields{j};
                TempVar2=LcSubsAsgn(TempVar2,TempVar3,[]);
              end;
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,TempVar2);
            otherwise,
            end;
          end;
        end;

        Selected.Columns=Selected.Columns+(1:length(Selected.Columns))-1;
        Display.Element(2)=Display.Element(2)+sum(Selected.Columns<Display.Element(2));
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');

      elseif cmd(2)==4,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% THE CREATE COLUMN OPTION WAS SELECTED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dataset=1;
        NDatasets=length(Display.Data);
        First=1;
        while (~isempty(Options.Multiple) && (Dataset<=NDatasets)) || First,
          First=0;
          TempSubsIndex=Display.Index(end).Ref;
          if ~isempty(Options.Multiple),
            TempSubsIndex(end+1)=TempSubsIndex(end);
            TempSubsIndex(end-1).type='{}';
            TempSubsIndex(end-1).subs={Dataset};
            Dataset=Dataset+1;
          end;
          TempSubsIndex(end).subs(:)={':'};
          TempVar=LcSubsRef(EditedData,TempSubsIndex(1:(end-1)));
          [TempSubsIndex(end).subs{size(TempVar)==0}]=deal(1);
          switch class(TempVar),
          case {'char'},
            EditedData=LcSubsAsgn(EditedData,TempSubsIndex,' ');
          case {'cell'},
            EditedData=LcSubsAsgn(EditedData,TempSubsIndex,{[]});
          case {'double','uint8','sparse'},
            EditedData=LcSubsAsgn(EditedData,TempSubsIndex,0);
          case {'struct'},
            TempVar2=size(LcSubsRef(EditedData,TempSubsIndex(1:(end-1))));
            TempVar3=TempSubsIndex(end).subs; % backup the subscript index
            for j=1:length(TempSubsIndex(end).subs), % replace ':' by dimension
              if strcmp(TempSubsIndex(end).subs{j},':'),
                if j>length(TempVar2),
                  TempSubsIndex(end).subs{j}=1;
                else
                  TempSubsIndex(end).subs{j}=TempVar2(j);
                end;
              end;
            end;
            Fields=fieldnames(TempVar);
            for j=1:length(Fields),
              F.type='.';
              F.subs=Fields{j};
              EditedData=LcSubsAsgn(EditedData,[TempSubsIndex F],[]);
            end;
            TempSubsIndex(end).subs=TempVar3; % reset back to ':'
          otherwise,
          end;
          OtherZeroDims=(size(TempVar)==0); % reset other dimensions that are empty back to emptiness
          OtherZeroDims(Display.Dimensions(2))=false;  %logical(0);
          if any(OtherZeroDims),
            [TempSubsIndex(end).subs{OtherZeroDims}]=deal([]);
          end;
          EditedData=LcSubsAsgn(EditedData,TempSubsIndex(1:(end-1)),LcSubsRef(EditedData,TempSubsIndex));
        end;
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');

      else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% UNKNOWN OPTION OF THE COLUMN MENU
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Str='Unknown column menu command';
        fprintf(1,[Str ':\n']);
        disp(cmd);
        uiwait(msgbox(Str,'modal'));
        set(Handles.General.MessageWindow,'string',[Str '.']);
      end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% THE COLUMN MENU WAS USED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif cmd(1)==-3,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% THE ROW MENU WAS USED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if cmd(2)==1,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% THE DELETE ROW(S) OPTION WAS SELECTED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dataset=1;
        NDatasets=length(Display.Data);
        First=1;
        while (~isempty(Options.Multiple) && (Dataset<=NDatasets)) || First,
          First=0;
          TempSubsIndex=Display.Index(end).Ref;
          if ~isempty(Options.Multiple),
            TempSubsIndex(end+1)=TempSubsIndex(end);
            TempSubsIndex(end-1).type='{}';
            TempSubsIndex(end-1).subs={Dataset};
            Dataset=Dataset+1;
          end;
          TempSubsIndex(end).subs(:)={':'};
          TempSubsIndex(end).subs{Display.Dimensions(1)}=Selected.Rows;

          TempVar=LcSubsRef(EditedData,TempSubsIndex(1:(end-1)));
          switch class(TempVar),
          case {'struct'},
            % To correct bugs in Matlab 5.2
            % A(5,1).a=1; A(:,1)=[]; returns a 1x0 struct
            % B=ones(5,1); B(:,1)=[]; returns a 5x0 matrix
            % A a 5x0 struct with field a; A(2,:)=[]; returns: ??? Improper null index in null matrix.
            % B=ones(5,0); B(2,:)=[]; returns a 4x0 matrix

            if ~isempty(TempVar),
              TempVar2=size(LcSubsRef(EditedData,TempSubsIndex(1:(end-1))));
              TempVar2(Display.Dimensions(1))=TempVar2(Display.Dimensions(1))-length(Selected.Rows);
              EditedData=subsasgn(EditedData,TempSubsIndex,[]);
              TempVar3=LcSubsRef(EditedData,TempSubsIndex(1:(end-1)));
              if ~isequal(TempVar2,TempVar3),
                EditedData=LcSubsAsgn(EditedData,TempSubsIndex(1:(end-1)),reshape(LcSubsRef(EditedData,TempSubsIndex(1:(end-1))),TempVar2));
              end;
            else
              TempVar2=size(LcSubsRef(EditedData,TempSubsIndex(1:(end-1))));
              TempVar2(Display.Dimensions(1))=TempVar2(Display.Dimensions(1))-length(Selected.Rows);
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex(1:(end-1)),reshape(LcSubsRef(EditedData,TempSubsIndex(1:(end-1))),TempVar2));
            end;
          otherwise,
            EditedData=subsasgn(EditedData,TempSubsIndex,[]);
          end;
        end;

        Display.Element(1)=max(1,Display.Element(1)-sum(Selected.Rows<=Display.Element(1)));
        Selected.Rows=[];
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');

      elseif cmd(2)==2,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% THE INSERT ROW ABOVE OPTION WAS SELECTED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dataset=1;
        NDatasets=length(Display.Data);
        First=1;
        while (~isempty(Options.Multiple) && (Dataset<=NDatasets)) || First,
          First=0;
          TempSubsIndex=Display.Index(end).Ref;
          if ~isempty(Options.Multiple),
            TempSubsIndex(end+1)=TempSubsIndex(end);
            TempSubsIndex(end-1).type='{}';
            TempSubsIndex(end-1).subs={Dataset};
            Dataset=Dataset+1;
          end;
          TempSubsIndex(end).subs(:)={':'};

          TempVar=LcSubsRef(EditedData,TempSubsIndex(1:(end-1)));
          if isempty(TempVar),
            TempVar2=size(TempVar);
            TempVar2(Display.Dimensions(1))=TempVar2(Display.Dimensions(1))+length(Selected.Rows);
            EditedData=LcSubsAsgn(EditedData,TempSubsIndex(1:(end-1)),reshape(TempVar,TempVar2));
          else
            TempSubsIndex(end).subs{Display.Dimensions(1)}=setdiff(1:(size(TempVar,Display.Dimensions(1))+length(Selected.Rows)),Selected.Rows+(1:length(Selected.Rows))-1);
            EditedData=LcSubsAsgn(EditedData,TempSubsIndex,TempVar);
            TempSubsIndex(end).subs{Display.Dimensions(1)}=Selected.Rows+(1:length(Selected.Rows))-1;
            switch class(TempVar),
            case {'char'},
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,' ');
            case {'cell'},
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,{[]});
            case {'double','uint8','sparse'},
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,0);
            case {'struct'},
              TempVar2=[];
              Fields=fieldnames(TempVar);
              for j=1:length(Fields),
                TempVar3.type='.';
                TempVar3.subs=Fields{j};
                TempVar2=LcSubsAsgn(TempVar2,TempVar3,[]);
              end;
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,TempVar2);
            otherwise,
            end;
          end;
        end;

        Display.Element(1)=Display.Element(1)+sum(Selected.Rows<=Display.Element(1));
        Selected.Rows=Selected.Rows+(1:length(Selected.Rows));
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');

      elseif cmd(2)==3,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% THE INSERT ROW BELOW OPTION WAS SELECTED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dataset=1;
        NDatasets=length(Display.Data);
        First=1;
        while (~isempty(Options.Multiple) && (Dataset<=NDatasets)) || First,
          First=0;
          TempSubsIndex=Display.Index(end).Ref;
          if ~isempty(Options.Multiple),
            TempSubsIndex(end+1)=TempSubsIndex(end);
            TempSubsIndex(end-1).type='{}';
            TempSubsIndex(end-1).subs={Dataset};
            Dataset=Dataset+1;
          end;
          TempSubsIndex(end).subs(:)={':'};

          TempVar=LcSubsRef(EditedData,TempSubsIndex(1:(end-1)));
          if isempty(TempVar),
            TempVar2=size(TempVar);
            TempVar2(Display.Dimensions(1))=TempVar2(Display.Dimensions(1))+length(Selected.Rows);
            EditedData=LcSubsAsgn(EditedData,TempSubsIndex(1:(end-1)),reshape(TempVar,TempVar2));
          else
            TempSubsIndex(end).subs{Display.Dimensions(1)}=setdiff(1:(size(TempVar,Display.Dimensions(1))+length(Selected.Rows)),Selected.Rows+(1:length(Selected.Rows)));
            EditedData=LcSubsAsgn(EditedData,TempSubsIndex,TempVar);
            TempSubsIndex(end).subs{Display.Dimensions(1)}=Selected.Rows+(1:length(Selected.Rows));
            switch class(TempVar),
            case {'char'},
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,' ');
            case {'cell'},
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,{[]});
            case {'double','uint8','sparse'},
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,0);
            case {'struct'},
              TempVar2=[];
              Fields=fieldnames(TempVar);
              for j=1:length(Fields),
                TempVar3.type='.';
                TempVar3.subs=Fields{j};
                TempVar2=LcSubsAsgn(TempVar2,TempVar3,[]);
              end;
              EditedData=LcSubsAsgn(EditedData,TempSubsIndex,TempVar2);
            otherwise,
            end;
          end;
        end;

        Display.Element(1)=Display.Element(1)+sum(Selected.Rows<Display.Element(1));
        Selected.Rows=Selected.Rows+(1:length(Selected.Rows))-1;
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');

      elseif cmd(2)==4,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% THE CREATE ROW OPTION WAS SELECTED
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Dataset=1;
        NDatasets=length(Display.Data);
        First=1;
        while (~isempty(Options.Multiple) && (Dataset<=NDatasets)) || First,
          First=0;
          TempSubsIndex=Display.Index(end).Ref;
          if ~isempty(Options.Multiple),
            TempSubsIndex(end+1)=TempSubsIndex(end);
            TempSubsIndex(end-1).type='{}';
            TempSubsIndex(end-1).subs={Dataset};
            Dataset=Dataset+1;
          end;
          TempSubsIndex(end).subs(:)={':'};
          TempVar=LcSubsRef(EditedData,TempSubsIndex(1:(end-1)));
          [TempSubsIndex(end).subs{size(TempVar)==0}]=deal(1);
          switch class(TempVar),
          case {'char'},
            EditedData=LcSubsAsgn(EditedData,TempSubsIndex,' ');
          case {'cell'},
            EditedData=LcSubsAsgn(EditedData,TempSubsIndex,{[]});
          case {'double','uint8','sparse'},
            EditedData=LcSubsAsgn(EditedData,TempSubsIndex,0);
          case {'struct'},
            TempVar2=size(LcSubsRef(EditedData,TempSubsIndex(1:(end-1))));
            TempVar3=TempSubsIndex(end).subs; % backup the subscript index
            for j=1:length(TempSubsIndex(end).subs), % replace ':' by dimension
              if strcmp(TempSubsIndex(end).subs{j},':'),
                if j>length(TempVar2),
                  TempSubsIndex(end).subs{j}=1;
                else
                  TempSubsIndex(end).subs{j}=TempVar2(j);
                end;
              end;
            end;
            Fields=fieldnames(TempVar);
            for j=1:length(Fields),
              F.type='.';
              F.subs=Fields{j};
              EditedData=LcSubsAsgn(EditedData,[TempSubsIndex F],[]);
            end;
            TempSubsIndex(end).subs=TempVar3; % reset back to ':'
          otherwise,
          end;
          OtherZeroDims=(size(TempVar)==0); % reset other dimensions that are empty back to emptiness
          OtherZeroDims(Display.Dimensions(1))=false; %logical(0);
          if any(OtherZeroDims),
            [TempSubsIndex(end).subs{OtherZeroDims}]=deal([]);
            EditedData=LcSubsAsgn(EditedData,TempSubsIndex(1:(end-1)),LcSubsRef(EditedData,TempSubsIndex));
          end;
        end;

        LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');

      else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% UNKNOWN OPTION OF THE ROW MENU
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Str='Unknown row menu command';
        fprintf(1,[Str ':\n']);
        disp(cmd);
        set(Handles.General.MessageWindow,'string',[Str '.']);
      end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% THE ROW MENU WAS USED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif cmd(1)==-4,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% THE MULTIDIMENSIONAL MENU WAS USED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if cmd(2)==1, % select dimensions and indices
        Display.PreviousSelectedDims=Display.Index(end).Ref(end).subs;
        Display.Index(end).Ref(end)=[];
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 2],'top');
      end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% THE MULTIDIMENSIONAL MENU WAS USED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif cmd(1)==-5,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% THE CHARACTER MENU WAS USED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if cmd(2)==1,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% SWITCH TO CHARACTER MATRIX EDITING MODE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(Handles.Menu.Column.Main,'enable','on');
        set(Handles.Menu.Row.Main,'enable','on');
        Options.EditMode=Options.EditModes{1};
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');

      elseif cmd(2)==2,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% SWITCH TO LINE (ONE LINE STRINGS) EDITING MODE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(Handles.Menu.Column.Main,'enable','off');
        set(Handles.Menu.Row.Main,'enable','on');
        Options.EditMode=Options.EditModes{2};
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');

      elseif cmd(2)==3,
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% SWITCH TO TEXT (MULTI LINE STRING) EDITING MODE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        set(Handles.Menu.Column.Main,'enable','off');
        set(Handles.Menu.Row.Main,'enable','off');
        Options.EditMode=Options.EditModes{3};
        LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
      end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% THE CHARACTER MENU WAS USED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif cmd(1)==-6,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% A FIELD OF THE STRUCTURE / OR MULTIPLE EDITING MODE WAS EDITED OR CLICKED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if ~isempty(Options.Multiple),
        %%---------------------------------------------------------------------
        %% MULTIPLE editing mode
        %%---------------------------------------------------------------------
        cmd2=Display.MultipleField{Offset.Field+cmd(2),1};
        switch class(Display.Data{cmd2}),
        case {'char'},
          if ~Options.MultipleLines,
            set(Handles.Field.Field(cmd(2)),'foregroundcolor',Color.Message);
            temp=get(Handles.Field.Field(cmd(2)),'string');
            TempRef=Display.Index(end).Ref(end).subs;
            TempRef(Display.Dimensions)={Display.Element(1) Display.Element(2)};
            if length(temp)>1 && temp(1)=='#', % possibly ASCII code intended
              if temp(1)=='#' && all(abs(temp(2:length(temp)))>=48) & all(abs(temp(2:length(temp)))<=57),
                temp=sscanf(temp(2:length(temp)),'%i');
                if temp>65535,
                  Str=['Valid character codes are 0:65535.'];
                  uiwait(msgbox(Str,'modal'));
                  set(Handles.General.MessageWindow,'string',Str);
                else
                  Display.Data{cmd2}(TempRef{:})=char(temp);
                  EditedData=LcSubsAsgn(EditedData,Display.Index(1:(end-1)).Ref,Display.Data);
                end;
              else
                Str = 'Input must be a single character or # followed by a character code.';
                uiwait(msgbox(Str,'modal'));
                set(Handles.General.MessageWindow,'string',Str);
              end;
            elseif length(temp)==1,
              TempRef=Display.Index(end).Ref(end).subs;
              TempRef(Display.Dimensions)={Display.Element(1) Display.Element(2)};
              Display.Data{cmd2}(TempRef{:})=temp;
              EditedData=LcSubsAsgn(EditedData,Display.Index(end).Ref(1:(end-1)),Display.Data);
            elseif isempty(temp),
              TempRef=Display.Index(end).Ref(end).subs;
              TempRef(Display.Dimensions)={Display.Element(1) Display.Element(2)};
              Display.Data{cmd2}(TempRef{:})=char(0);
              EditedData=LcSubsAsgn(EditedData,Display.Index(end).Ref(1:(end-1)),Display.Data);
            else
              Str=['Input must be a single character or # followed by a character code.'];
              uiwait(msgbox(Str,'modal'));
              set(Handles.General.MessageWindow,'string',Str);
            end;
            if abs(Display.Data{cmd2}(TempRef{:}))<32 || ...
               abs(Display.Data{cmd2}(TempRef{:}))>255,
              Str=['#' num2str(abs(Display.Data{cmd2}(TempRef{:})))];
              set(Handles.Menu.Char.Main,'enable','off');
            else
              Str=Display.Data{cmd2}(TempRef{:});
            end;
            set( Handles.Field.Field(cmd(2)) ,'string', Str, ...
                                       'foregroundcolor',Color.Normal );
          else
            % Options.MultipleLines
            temp=deblank(get(Handles.Field.Field(cmd(2)),'string'));
            if length(temp)>size(Display.Data{cmd2},2),
              % Fit matrix to new line
              Display.Data{cmd2}(:,(size(Display.Data{cmd2},2)+1):length(temp))=' ';
              Display.Data{cmd2}(Display.Element(1),:)=temp;
            elseif length(temp)==size(Display.Data{cmd2},2),
              % No fitting necessary
              Display.Data{cmd2}(Display.Element(1),:)=temp;
            else,
              % Fit new text to matrix
              temp((length(temp)+1):size(Display.Data{cmd2},2))=' ';
              Display.Data{cmd2}(Display.Element(1),:)=temp;
            end;
            EditedData=LcSubsAsgn(EditedData,Display.Index(end).Ref(1:(end-1)),Display.Data);
          end;
        case {'double','uint8','sparse'},
          set(Handles.Field.Field(cmd(2)),'foregroundcolor',Color.Message);
          lasterr('');
          %temp=evalin(EvalWhere,get(Handles.Field.Field(cmd(2)),'string'),'NaN');
          temp=str2num(get(Handles.Field.Field(cmd(2)),'string'));
          Str=lasterr;
          if isempty(Str),
            if isequal(size(temp),[1 1]) && isnumeric(temp),
              Display.Data{cmd2}(Display.Element(1),Display.Element(2))=temp;
              EditedData=LcSubsAsgn(EditedData,Display.Index(1:(end-1)).Ref,Display.Data);
            else
              Str='Input must be a single numeric value.';
              uiwait(msgbox(Str,'modal'));
              set(Handles.General.MessageWindow,'string',Str);
            end;
          else
            Str=multiline(Str);
            Str=Str(size(Str,1),:);
            uiwait(msgbox(Str,'modal'));
            set(Handles.General.MessageWindow,'string',Str);
          end;
          set(Handles.Field.Field(cmd(2)),'string', sprintf('%g',double(Display.Data{cmd2}(Display.Element(1),Display.Element(2)))), ...
                                   'foregroundcolor',Color.Normal );
        case 'cell',
          switch get(Handles.Figure,'selectiontype'),
          case 'normal',
            TempVar=Display.Index(end).Ref(end).subs;
            TempVar(Display.Dimensions)={Display.Element(1) Display.Element(2)};
            Display.Index(end+1)=Display.Index(end);
            Display.Index(end).Ref(end).type='{}';
            Display.Index(end).Ref(end).subs={cmd2};
            Display.Index(end).Ref(end+1).type='{}';
            Display.Index(end).Ref(end).subs=TempVar;
            if isstruct(LcSubsRef(EditedData,Display.Index(end).Ref)),
              Display.Element=[1 1];
            end;
            LcStackUDF(Handles.Figure,'CommandStack',[NaN 1],'top');
          case 'alt',
            TempRef=Display.Index(end).Ref;
            TempVar=TempRef(end).subs;
            TempVar(Display.Dimensions)={Display.Element(1) Display.Element(2)};
            TempRef(end).type='{}';
            TempRef(end).subs={cmd2};
            TempRef(end+1).type='{}';
            TempRef(end).subs=TempVar;
            TempVar=LcSubsRef(EditedData,TempRef);
            switch class(TempVar),
            case {'double','uint8','char','sparse','cell','struct'},
              Str = 'Cell does not contain object.';
              uiwait(msgbox(Str,'modal'));
              set(Handles.General.MessageWindow,'string',Str);
            otherwise,
              Confirm=questdlg('Convert object to struct?', ...
                         'Conversion warning','Yes','No','No');
              if strcmp(Confirm,'Yes')
                EditedData=LcSubsAsgn(EditedData,TempRef,struct(TempVar));
                LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
              end
            end;
          end;
        case 'struct',
          switch get(Handles.Figure,'selectiontype'),
          case 'normal',
            TempVar=Display.Index(end).Ref(end).subs;
            TempVar(Display.Dimensions)={Display.Element(1) Display.Element(2)};
            Display.Index(end+1)=Display.Index(end);
            Display.Index(end).Ref(end).type='{}';
            Display.Index(end).Ref(end).subs={cmd2};
            Display.Index(end).Ref(end+1).type='()';
            Display.Index(end).Ref(end).subs=TempVar;
            Display.Index(end).Ref(end+1).type='.';
            Display.Index(end).Ref(end).subs=Display.MultipleField{Offset.Field+cmd(2),2};
            if isstruct(LcSubsRef(EditedData,Display.Index(end).Ref)),
              Display.Element=[1 1];
            end;
            LcStackUDF(Handles.Figure,'CommandStack',[NaN 1],'top');
          case 'alt',
            TempRef=Display.Index(end).Ref;
            TempVar=TempRef(end).subs;
            TempVar(Display.Dimensions)={Display.Element(1) Display.Element(2)};
            TempRef(end).type='{}';
            TempRef(end).subs={cmd2};
            TempRef(end+1).type='()';
            TempRef(end).subs=TempVar;
            TempRef(end+1).type='.';
            TempRef(end).subs=Display.MultipleField{Offset.Field+cmd(2),2};
            TempVar=LcSubsRef(EditedData,TempRef);
            switch class(TempVar),
            case {'double','uint8','char','sparse','cell','struct'},
              Str = 'Cell does not contain object.';
              uiwait(msgbox(Str,'modal'));
              set(Handles.General.MessageWindow,'string',Str);
            otherwise,
              Confirm=questdlg('Convert object to struct?', ...
                         'Conversion warning','Yes','No','No');
              if strcmp(Confirm,'Yes')
                EditedData=LcSubsAsgn(EditedData,TempRef,struct(TempVar));
                LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
              end
            end;
          end;
        otherwise,
          Str=['Editing ',class(Display.Data{Offset.Field+cmd(2)}),' in multiple mode not supported.'];
          uiwait(msgbox(Str,'modal'));
          set(Handles.General.MessageWindow,'string',Str);
        end;

      else
        %%---------------------------------------------------------------------
        %% STRUCTURE editing mode
        %%---------------------------------------------------------------------
        switch get(Handles.Figure,'selectiontype'),
        case 'normal',
          TempVar=Display.Index(end).Ref(end).subs;
          TempVar(Display.Dimensions)={Display.Element(1) Display.Element(2)};
          Display.Index(end+1)=Display.Index(end);
          Display.Index(end).Ref(end).type='()';
          Display.Index(end).Ref(end).subs=TempVar;
          Display.Index(end).Ref(end+1).type='.';
          TempVar=fieldnames(Display.Data);
          TempVar=TempVar{Offset.Field+cmd(2)};
          Display.Index(end).Ref(end).subs=TempVar;
          if isstruct(getfield(Display.Data(Display.Element(1),Display.Element(2)),TempVar)),
            Display.Element=[1 1];
          end;
          LcStackUDF(Handles.Figure,'CommandStack',[NaN 1],'top');
        case 'alt',
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%% A STRUCTURE FIELD WAS RIGHT CLICKED ON
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          TempRef=Display.Index(end).Ref;
          TempRef(end).type='()';
          TempRef(end).subs(Display.Dimensions)={Display.Element(1) Display.Element(2)};
          TempRef(end+1).type='.';
          TempVar=fieldnames(Display.Data);
          TempVar=TempVar{Offset.Field+cmd(2)};
          TempRef(end).subs=TempVar;
          TempVar=LcSubsRef(EditedData,TempRef);
          switch class(TempVar),
          case {'double','uint8','char','sparse','cell','struct'},
            Str = 'Cell does not contain object.';
            uiwait(msgbox(Str,'modal'));
            set(Handles.General.MessageWindow,'string',Str);
          otherwise,
            Confirm=questdlg('Convert object to struct?', ...
                       'Conversion warning','Yes','No','No');
            if strcmp(Confirm,'Yes')
              EditedData=LcSubsAsgn(EditedData,TempRef,struct(TempVar));
              LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
            end
          end;
        end;
      end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% A FIELD OF THE STRUCTURE WAS EDITED OR CLICKED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif cmd(1)==-7,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% A STRUCTURE MENU WAS CLICKED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if cmd(2)==1,
        prompt={'Enter new field name:'};
        NewField=inputdlg(prompt,'',1,{''});

        if ~isempty(NewField),
          NewField=NewField{1};
          CharN=abs(NewField);
          NotValid=isempty(CharN);
          if ~NotValid,
            NotValid= (CharN(1)<65) | ...
                   any(CharN<48) | any((CharN>57) & (CharN<65)) | any(CharN>122) | ...
                   any((CharN>90) & (CharN<97) & (CharN~=95));
          end;

          if ~NotValid,
            ExistingFields=fieldnames(Display.Data);

            if isempty(strmatch(NewField,ExistingFields,'exact')),
              TempRef=Display.Index(end).Ref;
              TempRef(end)=[];
              if isempty(LcSubsRef(EditedData,TempRef)),
                TempVar2=fieldnames(Display.Data);
                TempVar2{end+1}=NewField;
                TempVar=reshape({},[length(TempVar2) size(LcSubsRef(EditedData,TempRef))]);
                EditedData=LcSubsAsgn(EditedData,TempRef,cell2struct(TempVar,TempVar2,1));
              else
                TempVar=Display.Index(end).Ref(end).subs;
                TempVar(Display.Dimensions)={Display.Element(1) Display.Element(2)};

                TempRef=Display.Index(end).Ref;
                TempRef(end).type='()';
                TempRef(end).subs=TempVar;
                TempRef(end+1).type='.';
                TempRef(end).subs=NewField;

                EditedData=LcSubsAsgn(EditedData,TempRef,[]);
              end;
              LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
            else
              Str = 'Field already exists.';
              uiwait(msgbox(Str,'modal'));
              set(Handles.General.MessageWindow,'string',Str);
            end;

          else
            Str = 'Invalid field name specified.';
            uiwait(msgbox(Str,'modal'));
            set(Handles.General.MessageWindow,'string',Str);
          end;
        end;
      end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% A STRUCTURE MENU WAS CLICKED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif cmd(1)==-8,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% A STRUCTURE FIELD LABEL WAS CLICKED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      TempVar=fieldnames(Display.Data);
      NoClear=length(TempVar)==1;
      TempVar=TempVar{Offset.Field+cmd(2)};
      OldField=TempVar;
      if NoClear,
        prompt={'Enter new field name:'};
      else
        prompt={'Enter new field name (clear to remove field):'};
      end;
      NewField=inputdlg(prompt,'',1,{TempVar});

      if ~isempty(NewField),
        NewField=NewField{1};
        CharN=abs(NewField);
        if (~NoClear) && isempty(CharN), % field name remove -> remove field
          TempRef=Display.Index(end).Ref;
          TempRef(end)=[];

          TempVar=struct2cell(LcSubsRef(EditedData,TempRef));
          if isempty(TempVar),
            TempVar=reshape(TempVar,[length(TempVar2) size(LcSubsRef(EditedData,TempRef))]);
          end;
          TempVar2=size(TempVar);
          TempVar(Offset.Field+cmd(2),:)=[];
          TempVar2(1)=TempVar2(1)-1;
          TempVar=reshape(TempVar,TempVar2);

          TempVar2=fieldnames(Display.Data);
          TempVar2(Offset.Field+cmd(2))=[];
          EditedData=LcSubsAsgn(EditedData,TempRef,cell2struct(TempVar,TempVar2,1));
          LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
        else % change field name
          NotValid=isempty(CharN);
          if ~NotValid,
            NotValid= (CharN(1)<65) | ...
                   any(CharN<48) | any((CharN>57) & (CharN<65)) | any(CharN>122) | ...
                   any((CharN>90) & (CharN<97) & (CharN~=95));
          end;

          if ~NotValid,
            ExistingFields=fieldnames(Display.Data);

            if isempty(strmatch(NewField,ExistingFields,'exact')),
              TempRef=Display.Index(end).Ref;
              TempRef(end)=[];

              TempVar=struct2cell(LcSubsRef(EditedData,TempRef));
              TempVar2=fieldnames(Display.Data);
              if isempty(TempVar),
                TempVar=reshape(TempVar,[length(TempVar2) size(LcSubsRef(EditedData,TempRef))]);
              end;
              TempVar2{Offset.Field+cmd(2)}=NewField;
              EditedData=LcSubsAsgn(EditedData,TempRef,cell2struct(TempVar,TempVar2,1));
              LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
            elseif strmatch(NewField,OldField),
              % no renaming: OldFieldName=NewFieldName
              % nor error message
            else
              Str = 'Field already exists.';
              uiwait(msgbox(Str,'modal'));
              set(Handles.General.MessageWindow,'string',Str);
            end;

          else
            Str = 'Invalid field name specified.';
            uiwait(msgbox(Str,'modal'));
            set(Handles.General.MessageWindow,'string',Str);
          end;
        end;
      end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% A STRUCTURE FIELD LABEL WAS CLICKED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif cmd(1)==-9,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CONVERT MENU WAS SELECTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      TempRef=Display.Index(end).Ref;
      TempRef(end)=[];
      TempVar=LcSubsRef(EditedData,TempRef);
      switch cmd(2), % to ...
      case 1, % double
        switch class(TempVar),
        case {'uint8','char'},
          EditedData=LcSubsAsgn(EditedData,TempRef,double(TempVar));
          LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
        end;
      case 2, % uint8
        switch class(TempVar),
        case {'double','char'},
          EditedData=LcSubsAsgn(EditedData,TempRef,uint8(TempVar));
          LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
        end;
      case 3, % sparse
      case 4, % char
        switch class(TempVar),
        case 'double',
          EditedData=LcSubsAsgn(EditedData,TempRef,char(TempVar));
          LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
        case 'uint8',
          EditedData=LcSubsAsgn(EditedData,TempRef,char(double(TempVar)));
          LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
        end;
      case 5, % cell
        switch class(TempVar),
        case {'double','uint8','char'},
          EditedData=LcSubsAsgn(EditedData,TempRef,num2cell(TempVar));
          LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
        end;
      case 6, % struct
      end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% A CONVERT MENU WAS SELECTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    else

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% UNKNOWN COMMAND ENCOUNTERED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Str='Unknown command';
      fprintf(1,[Str ':\n']);
      disp(cmd);
      set(Handles.General.MessageWindow,'string',[Str '.']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% UNKNOWN COMMAND ENCOUNTERED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end;

  end;

%%*************************************************************************************************
%%%% END OF WHILE COMMANDS ON STACK LOOP
%%*************************************************************************************************

%%*************************************************************************************************
%%%% RESET POINTER
%%*************************************************************************************************

  if ishandle(Handles.Figure),
    set(Handles.Figure,'pointer','arrow');
  end;

end;

%*-------------------------------------------------------------------------------------------------
%* END OF EDITING LOOP 
%*-------------------------------------------------------------------------------------------------

%*-------------------------------------------------------------------------------------------------
%* DELETE FIGURE IF IT STILL EXISTS
%*-------------------------------------------------------------------------------------------------

if ishandle(Handles.Figure),
  delete(Handles.Figure);
end;

%*-------------------------------------------------------------------------------------------------
%* CREATE APPROPRIATE OUTPUT
%*-------------------------------------------------------------------------------------------------

if nargout>1,
  varargout=EditedData;
else
  varargout{1}=EditedData;
end;


function [Handles,Color]=LcUiEdit(DisplayNumber)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INITIALIZE WINDOW 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Fld_NumFld=DisplayNumber.Fields;
Fld_LabelHeight=23;
Fld_Height=Fld_NumFld*Fld_LabelHeight;
Fld_Scrollbar=Fld_LabelHeight;
Fld_Button=28;

Lst_VertOffset=0;
Lst_HorzOffset=0;
Lst_ButtonHeight=Fld_LabelHeight;
Lst_Height=3*Lst_ButtonHeight;
Lst_TypeWidth=Lst_ButtonHeight;

Ind_HorzNumElement=DisplayNumber.Columns;
Ind_VertNumElement=DisplayNumber.Rows;
Ind_ElementWidth=80;
Ind_ElementHeight=Fld_LabelHeight;
Ind_Scrollbar=Fld_LabelHeight;
Ind_Label=Fld_LabelHeight;
Ind_Width =(Ind_HorzNumElement+1)*Ind_ElementWidth+Ind_Scrollbar+Ind_Label;
Ind_Height=(Ind_VertNumElement+1)*Ind_ElementHeight+Ind_Scrollbar+Ind_Label;
Ind_VertOffset=Lst_Height+5;
Ind_HorzOffset=0;

Fld_HorzOffset=0;
Fld_VertOffset=Ind_VertOffset+Ind_Height+5;
Fld_Width=Ind_Width;
Fld_LabelWidth=2*Ind_ElementWidth+Ind_Label-Fld_Button;

Lst_Width =Ind_Width;

Fig_Width =Ind_HorzOffset+Ind_Width;
Fig_Height=Fld_VertOffset+Fld_Height;

Color.Background=[1 1 1]*192/255;
Color.Edit=[1 1 1];
Color.Message=[1 0 0];
Color.Normal=[0 0 0];
Color.Selected=[166 202 240]/255;
Color.Struct=[1 1 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FIGURE INITIALIZATION %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ss = get(0,'ScreenSize');
swidth = ss(3);
sheight = ss(4);
left = (swidth-Fig_Width)/2;
bottom = (sheight-Fig_Height)/2;
rect = [left bottom Fig_Width Fig_Height];

Handles.Figure=figure('visible','off', ...
           'menu','none', ...
           'units','pixels', ...
           'color',Color.Background, ...
           'renderer','zbuffer', ...
           'inverthardcopy','off', ...
           'closerequestfcn','', ...
           'resize','off', ...
           'integerhandle','off', ...
           'numbertitle','off', ...
           'handlevisibility','off', ...
           'name','', ...
           'tag','IDEAS - GUI', ...
           'defaultuicontrolbackgroundcolor',Color.Background, ...
           'position',rect, ...
           'name','Edit window', ...
           'closerequestfcn','closereq');

ax1=axes(...
  'parent',Handles.Figure, ...
  'units','pixels', ...
  'position',[1 1 Fig_Width Fig_Height], ...
  'xlim',[0 Fig_Width-1],'ylim',[0 Fig_Height-1], ...
  'visible','off');

l1=LcHLine3D(0,Fld_VertOffset-3,Fig_Width,'parent',ax1);
l2=LcHLine3D(0,Ind_VertOffset-3,Fig_Width,'parent',ax1);

%%%%%%%%%%%%%%%%%%%%%
%%% UICONTEXTMENU %%%
%%%%%%%%%%%%%%%%%%%%%

Handles.UIContextMenu.Main=uicontextmenu('parent',Handles.Figure);

Handles.UIContextMenu.Convert=uimenu('parent',Handles.UIContextMenu.Main, ...
                                     'label','convert');

%%%%%%%%%%%%%%%%%%%%%%%
%%% INDEX SUBWINDOW %%%
%%%%%%%%%%%%%%%%%%%%%%%

% class type display
Handles.Index.ClassName=uicontrol(...
  'parent',Handles.Figure, ...
  'units','pixels', ...
  'position',[Ind_HorzOffset+Ind_Label Ind_VertOffset+Ind_Height-2*Ind_ElementHeight Ind_ElementWidth Ind_ElementHeight], ...
  'string','class', ...
  'buttondownfcn','', ...
  'style','pushbutton', ...
  'backgroundcolor',Color.Background, ...
  'enable','inactive', ...
  'userdata','2D Edit');

% column labels
Handles.Index.Label.Column(Ind_HorzNumElement)=0;
for j=1:Ind_HorzNumElement,
  str=LcColLabel(j);
  Handles.Index.Label.Column(j)=uicontrol(...
    'parent',Handles.Figure, ...
    'units','pixels', ...
    'callback',['md_edit(''MDEDIT-INTERNAL'',[0 ',num2str(j),'])'], ...
    'buttondownfcn','', ...
    'position',[Ind_HorzOffset+Ind_Label+j*Ind_ElementWidth Ind_VertOffset+Ind_Scrollbar+Ind_VertNumElement*Ind_ElementHeight Ind_ElementWidth Ind_ElementHeight], ...
    'string',str, ...
    'enable','off', ...
    'userdata','2D Edit');
end;

% row labels
Handles.Index.Label.Row(Ind_VertNumElement)=0;
for i=1:Ind_VertNumElement,
  str=num2str(i);
  Handles.Index.Label.Row(i)=uicontrol(...
    'parent',Handles.Figure, ...
    'units','pixels', ...
    'callback',['md_edit(''MDEDIT-INTERNAL'',[',str,' 0])'], ...
    'buttondownfcn','',...
    'position',[Ind_HorzOffset+Ind_Label Ind_VertOffset+Ind_Height-(2+i)*Ind_ElementHeight Ind_ElementWidth Ind_ElementHeight], ...
    'string',str, ...
    'enable','off', ...
    'userdata','2D Edit');
end;

%labels
xx=Handles.Index.Label.Column(1);
b1=LcBorder3D(Ind_HorzOffset,Ind_VertOffset+Ind_Scrollbar,Ind_Label,Ind_VertNumElement*Ind_ElementHeight,'parent',ax1,'userdata','2D Edit');
Handles.Index.Label.RowMain=text(Ind_HorzOffset+Ind_Label/2,Ind_VertOffset+Ind_Scrollbar+Ind_VertNumElement*Ind_ElementHeight/2,'-','parent',ax1);
set(Handles.Index.Label.RowMain,'fontangle',get(xx,'fontangle'), ...
       'fontsize',get(xx,'fontsize'), ...
       'fontweight',get(xx,'fontweight'), ...
       'color',get(xx,'foregroundcolor'), ...
       'fontunits',get(xx,'fontunits'), ...
       'verticalalignment','middle', ...
       'horizontalalignment','center', ...
       'rotation',90, ...
       'userdata','2D Edit');
b2=LcBorder3D(Ind_HorzOffset+Ind_Label+Ind_ElementWidth,Ind_VertOffset+Ind_Height-Ind_Label,Ind_HorzNumElement*Ind_ElementWidth,Ind_Label,'parent',ax1,'userdata','2D Edit');
Handles.Index.Label.ColumnMain=text(Ind_HorzOffset+Ind_Label+Ind_ElementWidth+Ind_HorzNumElement*Ind_ElementWidth/2,Ind_VertOffset+Ind_Height-Ind_Label+Ind_Label/2,'-','parent',ax1);
set(Handles.Index.Label.ColumnMain,'fontangle',get(xx,'fontangle'), ...
       'fontsize',get(xx,'fontsize'), ...
       'fontweight',get(xx,'fontweight'), ...
       'color',get(xx,'foregroundcolor'), ...
       'fontunits',get(xx,'fontunits'), ...
       'verticalalignment','middle', ...
       'horizontalalignment','center', ...
       'userdata','2D Edit');

% vertical (row) scroll bar
Handles.Index.Scrollbar.Row=uicontrol(...
  'parent',Handles.Figure, ...
  'units','pixels', ...
  'position',[Ind_HorzOffset+Ind_Width-Ind_Scrollbar Ind_VertOffset+Ind_Scrollbar Ind_Scrollbar Ind_VertNumElement*Ind_ElementHeight], ...
  'style','slider', ...
  'enable','off', ...
  'min',-1,'max',0, ...
  'userdata','2D Edit');
set(Handles.Index.Scrollbar.Row,'value',0, ...
 'callback',['md_edit(''MDEDIT-INTERNAL'',[-1 4])']);

% horizontal (column) scroll bar
Handles.Index.Scrollbar.Column=uicontrol(...
  'parent',Handles.Figure, ...
  'units','pixels', ...
  'position',[Ind_HorzOffset+Ind_Label+Ind_ElementWidth Ind_VertOffset Ind_ElementWidth*Ind_HorzNumElement Ind_Scrollbar], ...
  'style','slider', ...
  'enable','off', ...
  'min',0,'max',1, ...
  'userdata','2D Edit');
set(Handles.Index.Scrollbar.Column,'value',0, ...
 'callback',['md_edit(''MDEDIT-INTERNAL'',[-1 3])']);

% element items (editable or clickable)
Handles.Index.Elements(Ind_VertNumElement,Ind_HorzNumElement)=0;
for i=1:Ind_VertNumElement,
  for j=1:Ind_HorzNumElement,
    Handles.Index.Elements(i,j)=uicontrol(...
      'parent',Handles.Figure, ...
      'units','pixels', ...
      'position',[Ind_HorzOffset+Ind_Label+j*Ind_ElementWidth Ind_VertOffset+Ind_Height-(i+2)*Ind_ElementHeight Ind_ElementWidth Ind_ElementHeight], ...
      'style','edit', ...
      'backgroundcolor',Color.Background, ...
      'enable','off', ...
      'callback',['md_edit(''MDEDIT-INTERNAL'',[',num2str(i),' ',num2str(j),'])'], ...
      'userdata','2D Edit');
  end;
end;

Ind_Handles=findobj(Handles.Figure,'userdata','2D Edit');
set(Ind_Handles,'visible','off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DIMENSION EDIT SUBWINDOW %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dimension labels
for i=1:(Ind_VertNumElement+3),
  str=num2str(i);
  Handles.DimEdit.Label.Dim(i)=uicontrol(...
    'parent',Handles.Figure, ...
    'units','pixels', ...
    'position',[Ind_HorzOffset-2*Ind_ElementWidth+Ind_Width-Ind_Scrollbar Ind_VertOffset+Ind_Height-i*Ind_ElementHeight Ind_ElementWidth Ind_ElementHeight], ...
    'string',str, ...
    'enable','on', ...
    'buttondownfcn',['md_edit(''MDEDIT-INTERNAL'',[',str,' 0])'], ...
    'callback',['md_edit(''MDEDIT-INTERNAL'',[',str,' 0])'], ...
    'userdata','Dim Edit');
end;

% dimensions scroll bar
Handles.DimEdit.Scrollbar.Dim=uicontrol(...
  'parent',Handles.Figure, ...
  'units','pixels', ...
  'position',[Ind_HorzOffset+Ind_Width-Ind_Scrollbar Ind_VertOffset Ind_Scrollbar (Ind_VertNumElement+3)*Ind_ElementHeight], ...
  'style','slider', ...
  'enable','off', ...
  'min',-1,'max',0, ...
  'userdata','Dim Edit');
set(Handles.DimEdit.Scrollbar.Dim,'value',0, ...
 'callback',['md_edit(''MDEDIT-INTERNAL'',[-1 4])']);

% dimensions
Handles.DimEdit.Dims(Ind_VertNumElement+3)=0;
for i=1:(Ind_VertNumElement+3),
  Handles.DimEdit.Dim(i)=uicontrol(...
    'parent',Handles.Figure, ...
    'units','pixels', ...
    'position',[Ind_HorzOffset-Ind_ElementWidth+Ind_Width-Ind_Scrollbar Ind_VertOffset+Ind_Height-i*Ind_ElementHeight Ind_ElementWidth Ind_ElementHeight], ...
    'string','1', ...
    'style','edit', ...
    'horizontalalignment','right', ...
    'backgroundcolor',Color.Edit, ...
    'enable','on', ...
    'callback',['md_edit(''MDEDIT-INTERNAL'',[',num2str(i),' 1])'], ...
    'userdata','Dim Edit');
end;

%label
xx=Handles.DimEdit.Label.Dim(1);
b3=LcBorder3D(Ind_HorzOffset-2*Ind_ElementWidth+Ind_Width-Ind_Scrollbar-Ind_Label,Ind_VertOffset,Ind_Label,(Ind_VertNumElement+3)*Ind_ElementHeight,'parent',ax1);
set(b3,'userdata','Dim Edit');
Handles.DimEdit.Label.DimMain=text(Ind_HorzOffset-2*Ind_ElementWidth+Ind_Width-Ind_Scrollbar-Ind_Label+Ind_Label/2,Ind_VertOffset+(Ind_VertNumElement+3)*Ind_ElementHeight/2,'Dimensions','parent',ax1);
set(Handles.DimEdit.Label.DimMain,'fontangle',get(xx,'fontangle'), ...
       'fontsize',get(xx,'fontsize'), ...
       'fontweight',get(xx,'fontweight'), ...
       'color',get(xx,'foregroundcolor'), ...
       'fontunits',get(xx,'fontunits'), ...
       'verticalalignment','middle', ...
       'horizontalalignment','center', ...
       'rotation',90, ...
       'userdata','Dim Edit');

Dim_Handles=findobj(Handles.Figure,'userdata','Dim Edit');
set(Dim_Handles,'visible','off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SMALLEDIT SUBWINDOW %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% line labels
for i=1:(Ind_VertNumElement+3),
  str=num2str(i);
  Handles.LineEdit.Label.Line(i)=uicontrol(...
    'parent',Handles.Figure, ...
    'units','pixels', ...
    'position',[Ind_HorzOffset+Ind_Label Ind_VertOffset+Ind_Height-i*Ind_ElementHeight Ind_ElementWidth Ind_ElementHeight], ...
    'string',str, ...
    'enable','off', ...
    'buttondownfcn',['md_edit(''MDEDIT-INTERNAL'',[',str,' 0])'], ...
    'callback',['md_edit(''MDEDIT-INTERNAL'',[',str,' 0])'], ...
    'userdata','1D Edit');
end;

% vertical (lines) scroll bar
Handles.LineEdit.Scrollbar.Line=uicontrol(...
  'parent',Handles.Figure, ...
  'units','pixels', ...
  'position',[Ind_HorzOffset+Ind_Width-Ind_Scrollbar Ind_VertOffset Ind_Scrollbar (Ind_VertNumElement+3)*Ind_ElementHeight], ...
  'style','slider', ...
  'enable','off', ...
  'min',-1,'max',0, ...
  'userdata','1D Edit');
set(Handles.LineEdit.Scrollbar.Line,'value',0, ...
 'callback',['md_edit(''MDEDIT-INTERNAL'',[-1 4])']);

% string lines
Handles.LineEdit.Lines(Ind_VertNumElement+3)=0;
for i=1:(Ind_VertNumElement+3),
  Handles.LineEdit.Lines(i)=uicontrol(...
    'parent',Handles.Figure, ...
    'units','pixels', ...
    'position',[Ind_HorzOffset+Ind_Label+Ind_ElementWidth Ind_VertOffset+Ind_Height-i*Ind_ElementHeight Ind_Width-Ind_ElementWidth-Ind_Scrollbar-Ind_Label Ind_ElementHeight], ...
    'style','edit', ...
    'horizontalalignment','left', ...
    'backgroundcolor',Color.Background, ...
    'enable','off', ...
    'callback',['md_edit(''MDEDIT-INTERNAL'',[',num2str(i),' 1])'], ...
    'userdata','1D Edit');
end;

%label
xx=Handles.LineEdit.Label.Line(1);
b3=LcBorder3D(Ind_HorzOffset,Ind_VertOffset,Ind_Label,(Ind_VertNumElement+3)*Ind_ElementHeight,'parent',ax1);
set(b3,'userdata','1D Edit');
Handles.LineEdit.Label.LineMain=text(Ind_HorzOffset+Ind_Label/2,Ind_VertOffset+(Ind_VertNumElement+3)*Ind_ElementHeight/2,'-','parent',ax1);
set(Handles.LineEdit.Label.LineMain,'fontangle',get(xx,'fontangle'), ...
       'fontsize',get(xx,'fontsize'), ...
       'fontweight',get(xx,'fontweight'), ...
       'color',get(xx,'foregroundcolor'), ...
       'fontunits',get(xx,'fontunits'), ...
       'verticalalignment','middle', ...
       'horizontalalignment','center', ...
       'rotation',90, ...
       'userdata','1D Edit');

SmEd_Handles=findobj(Handles.Figure,'userdata','1D Edit');
set(SmEd_Handles,'visible','off');

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% BIGEDIT SUBWINDOW %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

Handles.TextEdit.Text=uicontrol(...
      'parent',Handles.Figure, ...
      'units','pixels', ...
      'position',[Ind_HorzOffset Ind_VertOffset Ind_Width Ind_Height], ...
      'style','edit', ...
      'backgroundcolor',Color.Background, ...
      'horizontalalignment','left', ...
      'max',2, ...
      'enable','off', ...
      'callback',['md_edit(''MDEDIT-INTERNAL'',[1 1])'], ...
      'userdata','0D Edit');

BigEd_Handles=Handles.TextEdit.Text;
set(BigEd_Handles,'visible','off');

Handles.Index.All={Ind_Handles,SmEd_Handles,BigEd_Handles};

%%%%%%%%%%%%%%%%%%%%%%%%
%%% STRUCT SUBWINDOW %%%
%%%%%%%%%%%%%%%%%%%%%%%%

% field labels (in lower half)
for i=1:Fld_NumFld,
  Handles.Field.Label.FieldName(i)=uicontrol(...
    'parent',Handles.Figure, ...
    'units','pixels', ...
    'position',[Fld_HorzOffset Fld_VertOffset+Fld_Height-i*Ind_ElementHeight Fld_LabelWidth Fld_LabelHeight], ...
    'enable','off', ...
    'string','', ...
    'callback',['md_edit(''MDEDIT-INTERNAL'',[-8 ',num2str(i),'])'], ...
    'userdata','Field');
  Handles.Field.Label.FieldType(i)=uicontrol(...
    'parent',Handles.Figure, ...
    'units','pixels', ...
    'position',[Fld_HorzOffset+Fld_LabelWidth Fld_VertOffset+Fld_Height-i*Ind_ElementHeight Fld_Button Fld_LabelHeight], ...
    'enable','off', ...
    'string',' ', ...
    'userdata','Field');
  Handles.Field.Field(i)=uicontrol(...
    'parent',Handles.Figure, ...
    'units','pixels', ...
    'position',[Fld_HorzOffset+Fld_LabelWidth+Fld_Button Fld_VertOffset+Fld_Height-i*Ind_ElementHeight Fld_Width-Fld_LabelWidth-Fld_Button-Fld_Scrollbar Fld_LabelHeight], ...
    'style','edit', ...
    'enable','inactive', ...
    'backgroundcolor',Color.Background, ...
    'userdata','Field');
end;

% vertical (field) scroll bar
Handles.Field.Scrollbar.Field=uicontrol(...
  'parent',Handles.Figure, ...
  'units','pixels', ...
  'position',[Fld_HorzOffset+Fld_Width-Fld_Scrollbar Fld_VertOffset Fld_Scrollbar Fld_NumFld*Fld_LabelHeight], ...
  'style','slider', ...
  'enable','off', ...
  'min',-1,'max',0, ...
  'userdata','Field');
set(Handles.Field.Scrollbar.Field,'value',0, ...
  'callback',['md_edit(''MDEDIT-INTERNAL'',[-1 6])']);

Handles.Field.All=findobj(Handles.Figure,'userdata','Field');
set(Handles.Field.All,'visible','off');

%%%%%%%%%%%%%%%%%%%%%%
%%% LIST SUBWINDOW %%%
%%%%%%%%%%%%%%%%%%%%%%

% reset and accept buttons
Handles.General.Reset=uicontrol(...
  'parent',Handles.Figure, ...
  'units','pixels', ...
  'position',[Lst_HorzOffset Lst_VertOffset Lst_Width/2 Lst_ButtonHeight], ...
  'enable','on', ...
  'string','reset', ...
  'callback',['md_edit(''MDEDIT-INTERNAL'',[-1 2])'], ...
  'userdata','List');

Handles.General.Accept=uicontrol(...
  'parent',Handles.Figure, ...
  'units','pixels', ...
  'position',[Lst_HorzOffset+Lst_Width/2 Lst_VertOffset Lst_Width/2 Lst_ButtonHeight], ...
  'enable','on', ...
  'string','accept', ...
  'callback',['md_edit(''MDEDIT-INTERNAL'',[-1 1])'], ...
  'userdata','List');

Str='(c) April 1999, H.R.A.Jagers';
Handles.General.MessageWindow=uicontrol(...
  'parent',Handles.Figure, ...
  'units','pixels', ...
  'position',[Lst_HorzOffset Lst_VertOffset+Lst_ButtonHeight Lst_Width Lst_ButtonHeight], ...
  'style','edit', ...
  'horizontalalignment','left', ...
  'enable','inactive', ...
  'backgroundcolor',Color.Edit, ...
  'foregroundcolor',Color.Message, ...
  'string',Str, ...
  'userdata','List');

% list
Handles.General.List=uicontrol(...
  'parent',Handles.Figure, ...
  'units','pixels', ...
  'position',[Lst_HorzOffset Lst_VertOffset+2*Lst_ButtonHeight Lst_Width Lst_ButtonHeight], ...
  'style','popupmenu', ...
  'string',' ', ...
  'enable','on', ...
  'horizontalalignment','left', ...
  'backgroundcolor',Color.Background, ...
  'callback',['md_edit(''MDEDIT-INTERNAL'',[-1 5])'], ...
  'userdata','List');

Lst_Handles=findobj(Handles.Figure,'userdata','List');
%set(Lst_Handles,'visible','off');

%%%%%%%%%%%%%%%%%%
%%% MENU ITEMS %%%
%%%%%%%%%%%%%%%%%%

Handles.Menu.Column.Main=uimenu(Handles.Figure, ...
   'label','&column', ...
   'enable','off');
Handles.Menu.Column.Delete=uimenu(Handles.Menu.Column.Main, ...
   'label','&delete column(s)', ...
   'enable','off', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-2 1])']);
Handles.Menu.Column.CreateOne=uimenu(Handles.Menu.Column.Main, ...
   'label','create &column', ...
   'enable','off', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-2 4])']);
Handles.Menu.Column.InsLeft=uimenu(Handles.Menu.Column.Main, ...
   'label','insert column &left', ...
   'enable','off', ...
   'separator','on', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-2 2])']);
Handles.Menu.Column.InsRight=uimenu(Handles.Menu.Column.Main, ...
   'label','insert column &right', ...
   'enable','off', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-2 3])']);

Handles.Menu.Row.Main=uimenu(Handles.Figure, ...
   'label','&row', ...
   'enable','off');
Handles.Menu.Row.Delete=uimenu(Handles.Menu.Row.Main, ...
   'label','&delete row(s)', ...
   'enable','off', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-3 1])']);
Handles.Menu.Row.CreateOne=uimenu(Handles.Menu.Row.Main, ...
   'label','create r&ow', ...
   'enable','off', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-3 4])']);
Handles.Menu.Row.InsAbove=uimenu(Handles.Menu.Row.Main, ...
   'label','insert row &above', ...
   'enable','off', ...
   'separator','on', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-3 2])']);
Handles.Menu.Row.InsBelow=uimenu(Handles.Menu.Row.Main, ...
   'label','insert row &below', ...
   'enable','off', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-3 3])']);

Handles.Menu.Data.Main=uimenu(Handles.Figure, ...
   'label','&data', ...
   'enable','on');

Handles.Menu.Dimen.Main=uimenu(Handles.Menu.Data.Main, ...
   'label','&multidimensional', ...
   'enable','off');
seldim=uimenu(Handles.Menu.Dimen.Main, ...
   'label','&select dimensions and indices', ...
   'enable','on', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-4 1])']);

Handles.Menu.Convert.Main=uimenu(Handles.Menu.Data.Main, ...
   'label','c&onvert to ...', ...
   'enable','on');
Handles.Menu.Convert.Double=uimenu(Handles.Menu.Convert.Main, ...
   'label','&double', ...
   'enable','off', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-9 1])']);
Handles.Menu.Convert.Uint8=uimenu(Handles.Menu.Convert.Main, ...
   'label','&uint8', ...
   'enable','off', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-9 2])']);
Handles.Menu.Convert.Sparse=uimenu(Handles.Menu.Convert.Main, ...
   'label','&sparse', ...
   'enable','off', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-9 3])']);
Handles.Menu.Convert.Char=uimenu(Handles.Menu.Convert.Main, ...
   'label','c&har', ...
   'enable','off', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-9 4])']);
Handles.Menu.Convert.Cell=uimenu(Handles.Menu.Convert.Main, ...
   'label','&cell', ...
   'enable','off', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-9 5])']);
Handles.Menu.Convert.Struct=uimenu(Handles.Menu.Convert.Main, ...
   'label','&struct', ...
   'enable','off', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-9 6])']);

Handles.Menu.Char.Main=uimenu(Handles.Menu.Data.Main, ...
   'label','c&haracter', ...
   'separator','on', ...
   'enable','off');
castr=uimenu(Handles.Menu.Char.Main, ...
   'label','&character array', ...
   'enable','on', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-5 1])']);
slstr=uimenu(Handles.Menu.Char.Main, ...
   'label','&single line strings', ...
   'enable','on', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-5 2])']);
mlstr=uimenu(Handles.Menu.Char.Main, ...
   'label','&multiline string', ...
   'enable','on', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-5 3])']);

Handles.Menu.Structure.Main=uimenu(Handles.Menu.Data.Main, ...
   'label','&structure', ...
   'enable','off');
Handles.Menu.Structure.NewField=uimenu(Handles.Menu.Structure.Main, ...
   'label','&new field', ...
   'enable','on', ...
   'callback',['md_edit(''MDEDIT-INTERNAL'',[-7 1])']);

%%%%%%%%%%%%%%%%%%%
%%% SHOW WINDOW %%%
%%%%%%%%%%%%%%%%%%%

set(Handles.Figure, ...
  'pointer','arrow', ...
  'visible','on');

qqq=findobj(Handles.Figure);
set(qqq,'interruptible','on');
drawnow;

LcSetUDF(Handles.Figure,'CommandStack',{[NaN 1]});


function [Selected,Display,Offset,Options]=LcUiEdit_update(DisplayNumber,Color,cmd,Handles,EditedData,Selected,Display,Offset,Options,options);
% LcUiEdit_UPDATE
%        cmd(2) should be 0, 1, or 2

if (cmd(2)==1),
  Selected.Columns=[];
  Selected.Rows=[];
  Offset.Row=0;
  Offset.Column=0;
  Offset.Field=0;
  Display.Element=[1 1];
end;
Display.Data=LcSubsRef(EditedData,Display.Index(end).Ref);
Display.TotalDataSize=size(Display.Data);
if ~ischar(Display.Data) | (cmd(2)~=0),
  Options.EditMode=Options.EditModes{1};
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DETERMINE CLASS OF DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TempClass=class(Display.Data);
switch TempClass,
case {'double','uint8','char','cell','struct'}, % 'sparse' gives problems for index e.g. b{1}(:,:) where b{1} is sparse

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% A DEFAULT CLASS OF MATLAB: NUMERIC, CHARACTER, CELL OR STRUCT
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%---------------------------------------------------------------------------
  %% Select the mode of reinitialization.
  %%---------------------------------------------------------------------------
  switch cmd(2),
  case {1,2},
    %%-------------------------------------------------------------------------
    %% Complete reinitialization.
    %% cmd(2)==2 is used to force multidimensional selection.
    %%-------------------------------------------------------------------------

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Process options
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Set default options

    Options.Multiple={};
    Options.MultipleLines=0;
    Options.DimElmLabel={};
    Options.DimFixed=zeros(1,ndims(Display.Data));
    Options.DimLabel={};

    %% Check for options specified by the user

    if length(Display.Index(end).Ref)==0,
      TempOpt=options;
    else, % options for second level have not yet been programmed
      TempOpt={};
    end;
    Options=LcProcOpt(Display,Options,TempOpt);

    %% When editing vectors and string vectors (character matrix interpreted
    %% as a set of strings) you can't add columns; since there would be no
    %% difference between the number of columns and the length of the strings
    %% of the vector.

    if Options.MultipleLines,
      Options.DimFixed(2)=1;
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% End processing options
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%-------------------------------------------------------------------------
    %% Get the subscript reference to extract and display the data.
    %%-------------------------------------------------------------------------
    if ~isempty(Options.Multiple),
      MultiDim=ndims(Display.Data{1})>2;
      Display.TotalDataSize=size(Display.Data{1});
      if Options.MultipleLines,
        Display.TotalDataSize(2)=1;
      end;
    else,
      MultiDim=ndims(Display.Data)>2;
      Display.TotalDataSize=size(Display.Data);
    end;

    if (cmd(2)==2) | MultiDim,

      Display.Element=[1 1];
      %% a dataset with more than two dimensions.
      LengthSI=length(Display.Index(end).Ref)+1;
      Display.Index(end).Ref(LengthSI).type='()';
      Display.Index(end).Ref(LengthSI).subs={};
      set(Handles.Figure,'pointer','arrow');
      if sum(Display.TotalDataSize==0)>2,
        Str=sprintf('Too many empty dimensions to edit.');
        uiwait(msgbox(Str,'modal'));
        set(Handles.General.MessageWindow,'string',Str);
        Display.Index(end).Ref=Display.Index(end).Ref(1:(end-1));
        if length(Display.Index(end).Ref)==0,
          LcStackUDF(Handles.Figure,'CommandStack',[-1 1],'top');
        else, % Go back to previous level
          Display.Index(end)=[];
          LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
        end;
        return;
      end;
      [Display.Dimensions,Display.Index(end).Ref(LengthSI).subs]=LcSeldim(2,length(Display.TotalDataSize),Display.TotalDataSize,Display.PreviousSelectedDims);
      Display.PreviousSelectedDims={};
      set(Handles.Figure,'pointer','watch');
      set(Handles.Menu.Dimen.Main,'enable','on');

    else,

      %% a dataset with two dimensions.
      Display.Dimensions=[1 2];
      LengthSI=length(Display.Index(end).Ref)+1;
      Display.Index(end).Ref(LengthSI).type='()';
      Display.Index(end).Ref(LengthSI).subs={':',':'};
      set(Handles.Menu.Dimen.Main,'enable','on');

    end;

    if ~isempty(Options.Multiple),
      Display.Data=LcSubsRef(EditedData,Display.Index(end).Ref(1:(end-1)));
      % No squeezing necessary for in case of MULTIPLE editing mode
      % since only cell vectors are supported.
    else,
      Display.Data=LcSubsRef(EditedData,Display.Index(end).Ref);

      %%-------------------------------------------------------------------------
      %% Squeeze the data into a 2D matrix.
      %%-------------------------------------------------------------------------
      TempDim=size(Display.Data,Display.Dimensions(1));
      Display.Data=squeeze(Display.Data);
      if size(Display.Data,1)~=TempDim,
        Display.Data=transpose(Display.Data);
      end;
    end;

    %%-------------------------------------------------------------------------
    %% Enable the row menu when the number of rows is not fixed
    %%-------------------------------------------------------------------------
    if Display.Dimensions(1)>length(Options.DimFixed),
      set(Handles.Menu.Row.Main,'enable','on');
    elseif Options.DimFixed(Display.Dimensions(1)),
      set(Handles.Menu.Row.Main,'enable','off');
    else,
      set(Handles.Menu.Row.Main,'enable','on');
    end;

    %%-------------------------------------------------------------------------
    %% Enable the column menu when the number of columns is not fixed.
    %%-------------------------------------------------------------------------
    if Display.Dimensions(2)>length(Options.DimFixed),
      set(Handles.Menu.Column.Main,'enable','on');
    elseif Options.DimFixed(Display.Dimensions(2)),
      set(Handles.Menu.Column.Main,'enable','off');
    else,
      set(Handles.Menu.Column.Main,'enable','on');
    end;

  otherwise, % cmd(2) is member {0}

    %%-------------------------------------------------------------------------
    %% Incomplete reinitialization.
    %%-------------------------------------------------------------------------
    if ~isempty(Options.Multiple),
      Display.Data=LcSubsRef(EditedData,Display.Index(end).Ref(1:(end-1)));
      % No squeezing necessary for in case of MULTIPLE editing mode
      % since only cell vectors are supported.
    else,
      Display.Data=LcSubsRef(EditedData,Display.Index(end).Ref);

      %%-------------------------------------------------------------------------
      %% Squeeze the data into a 2D matrix.
      %%-------------------------------------------------------------------------
      TempDim=size(Display.Data,Display.Dimensions(1));
      Display.Data=squeeze(LcSubsRef(EditedData,Display.Index(end).Ref));
      if size(Display.Data,1)~=TempDim,
        Display.Data=transpose(Display.Data);
      end;
    end;
  end;

  %%-------------------------------------------------------------------------
  %% Enable/disable convert menu items
  %%-------------------------------------------------------------------------
  TempHandles=struct2cell(Handles.Menu.Convert);
  TempHandles=[TempHandles{:}];
  set(TempHandles,'enable','off');
  set(Handles.Menu.Convert.Main,'enable','on');
  switch class(Display.Data),
  case 'double',
    set(Handles.Menu.Convert.Uint8,'enable','on');
    set(Handles.Menu.Convert.Char,'enable','on');
    set(Handles.Menu.Convert.Cell,'enable','on');
  case 'uint8',
    set(Handles.Menu.Convert.Double,'enable','on');
    set(Handles.Menu.Convert.Char,'enable','on');
    set(Handles.Menu.Convert.Cell,'enable','on');
  case 'char',
    set(Handles.Menu.Convert.Double,'enable','on');
    set(Handles.Menu.Convert.Uint8,'enable','on');
    set(Handles.Menu.Convert.Cell,'enable','on');
  case 'sparse',
  case 'cell',
  case 'struct',
  otherwise,
    set(Handles.Menu.Convert.Struct,'enable','on');
  end;

  %%-------------------------------------------------------------------------
  %% Store Display.Data in TempContent and change Display.Data
  %% to one of the datasets edited in MULTIPLE mode.
  %% MULTIPLE mode is always combined with EditModes{1} % array
  %% Display.Data is restored in the block updating the
  %% 2D specific elements of the interface.
  %%-------------------------------------------------------------------------
  if ~isempty(Options.Multiple),
    TempContent=Display.Data;
    Display.Data=LcSubsRef(Display.Data{1},Display.Index(end).Ref(end));

    TempDim=size(Display.Data,Display.Dimensions(1));
    Display.Data=squeeze(Display.Data);
    if size(Display.Data,1)~=TempDim,
      Display.Data=transpose(Display.Data);
    end;

    if Options.MultipleLines & ~isempty(Display.Data),
      Display.Data=Display.Data(:,1);
    end;
  end;

  %%---------------------------------------------------------------------------
  %% Switch to EditMode 1 if the data is character and if it has more
  %% than 2 dimensions and/or contains non-printing characters
  %%---------------------------------------------------------------------------
  if ischar(Display.Data) & ( length(Display.Index(end).Ref(end).subs)>2 | ...
                              any((Display.Data(:)<32) | (Display.Data(:)>255)) ),
    Options.EditMode=Options.EditModes{1};
  end;

  %%---------------------------------------------------------------------------
  %% If the data doesn't contain a row, enable the menu to create one.
  %%---------------------------------------------------------------------------
  if size(Display.Data,1)==0,
    set(Handles.Menu.Row.CreateOne,'enable','on');
  else,
    set(Handles.Menu.Row.CreateOne,'enable','off');
  end;

  %%---------------------------------------------------------------------------
  %% Enable the menu to insert other rows above and below.
  %%---------------------------------------------------------------------------
  if length(Selected.Rows)>0,
    set(Handles.Menu.Row.InsAbove,'enable','on');
    set(Handles.Menu.Row.InsBelow,'enable','on');
  else,
    set(Handles.Menu.Row.InsAbove,'enable','off');
    set(Handles.Menu.Row.InsBelow,'enable','off');
  end;

  %%---------------------------------------------------------------------------
  %% If one or more rows are selected, enable the menu to delete them.
  %%---------------------------------------------------------------------------
  if length(Selected.Rows)>0,
    set(Handles.Menu.Row.Delete,'enable','on');
    if length(Selected.Rows)>1,
      set(Handles.Menu.Row.Delete,'label','&delete rows');
    else,
      set(Handles.Menu.Row.Delete,'label','&delete row');
    end;
  else,
    set(Handles.Menu.Row.Delete,'enable','off');
  end;

  %%---------------------------------------------------------------------------
  %% Fill the listbox containing the settings of previous layers.
  %%---------------------------------------------------------------------------
  for i=1:length(Display.Index),
    if i==1,
      Str=LcRef2Str(Display.Index(1).Ref);
    else,
      Str=str2mat(Str,LcRef2Str(Display.Index(i).Ref));
    end;
  end;  
  set(Handles.General.List,'string',Str,'value',length(Display.Index));

  %%---------------------------------------------------------------------------
  %% Remove all information about previous editing in Multiple mode.
  %%---------------------------------------------------------------------------
  Display.MultipleField={};

  switch Options.EditMode,
  case Options.EditModes(1), % array
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Now update the 2D specific elements of the interface.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%-------------------------------------------------------------------------
    %% If the data doesn't contain a column, enable the menu to create one.
    %%-------------------------------------------------------------------------
    if size(Display.Data,2)==0,
      set(Handles.Menu.Column.CreateOne,'enable','on');
    else,
      set(Handles.Menu.Column.CreateOne,'enable','off');
    end;

    %%-------------------------------------------------------------------------
    %% Enable the menu to insert other columns to the left and right.
    %%-------------------------------------------------------------------------
    if length(Selected.Columns)>0,
      set(Handles.Menu.Column.InsRight,'enable','on');
      set(Handles.Menu.Column.InsLeft,'enable','on');
    else,
      set(Handles.Menu.Column.InsRight,'enable','off');
      set(Handles.Menu.Column.InsLeft,'enable','off');
    end;

    %%-------------------------------------------------------------------------
    %% If one or more columns are selected, enable the menu to delete them.
    %%-------------------------------------------------------------------------
    if length(Selected.Columns)>0,
      set(Handles.Menu.Column.Delete,'enable','on');
      if length(Selected.Columns)>1,
        set(Handles.Menu.Column.Delete,'label','&delete columns');
      else,
        set(Handles.Menu.Column.Delete,'label','&delete column');
      end;
    else,
      set(Handles.Menu.Column.Delete,'enable','off');
    end;

    %%-------------------------------------------------------------------------
    %% Set the class name in the upper left corner of the editing window.
    %%-------------------------------------------------------------------------
    if ~isempty(Options.Multiple),
      set(Handles.Index.ClassName,'string','*');
    else,
      set(Handles.Index.ClassName,'string',class(Display.Data));
    end;

    %%-------------------------------------------------------------------------
    %% Set the label of the rows.
    %%-------------------------------------------------------------------------
    if isempty(Options.DimElmLabel),
      Str=sprintf('dimension %i',Display.Dimensions(1));
    else,
      TempVar=Options.DimElmLabel{Display.Dimensions(1)};
      if ischar(TempVar),
        Str=TempVar;
      else,
        Str=sprintf('dimension %i',Display.Dimensions(1));
      end;
    end;
    set(Handles.Index.Label.RowMain,'string',Str);

    %%-------------------------------------------------------------------------
    %% Set the label of the columns.
    %%-------------------------------------------------------------------------
    if isempty(Options.DimElmLabel),
      Str=sprintf('dimension %i',Display.Dimensions(2));
    else,
      TempVar=Options.DimElmLabel{Display.Dimensions(2)};
      if ischar(TempVar),
        Str=TempVar;
      else,
        Str=sprintf('dimension %i',Display.Dimensions(2));
      end;
    end;
    set(Handles.Index.Label.ColumnMain,'string',Str);

    %%-------------------------------------------------------------------------
    %% Set the fields of the editing controls.
    %%-------------------------------------------------------------------------
    Offset.Row=min(Offset.Row,max(size(Display.Data,1)-DisplayNumber.Rows,0));
    Offset.Column=min(Offset.Column,max(size(Display.Data,2)-DisplayNumber.Columns,0));
    for i=1:DisplayNumber.Rows,
      for j=1:DisplayNumber.Columns,
        if all([i j]<=size(Display.Data)),
          if isempty(find(~(Selected.Rows-i-Offset.Row))) & isempty(find(~(Selected.Columns-j-Offset.Column))),
            TempColor=Color.Edit;
          else,
            TempColor=Color.Selected;
          end;
          if isstruct(Display.Data) | ~isempty(Options.Multiple),
            Str=sprintf('[%i,%i]',Offset.Row+i,Offset.Column+j);
          elseif ischar(Display.Data),
            if abs(Display.Data(Offset.Row+i,Offset.Column+j))<32 | ...
               abs(Display.Data(Offset.Row+i,Offset.Column+j))>255,
              Str=sprintf('#%i',abs(Display.Data(Offset.Row+i,Offset.Column+j)));
            else,
              Str=Display.Data(Offset.Row+i,Offset.Column+j);
            end;
          elseif iscell(Display.Data),
            TempVar=Display.Data{Offset.Row+i,Offset.Column+j};
            Str=LcVal2Str(TempVar);
          else,
            Str=sprintf('%g',double(Display.Data(Offset.Row+i,Offset.Column+j)));
          end;
          if iscell(Display.Data) | isstruct(Display.Data) | ~isempty(Options.Multiple),
            if (isstruct(Display.Data) | ~isempty(Options.Multiple)) & all([Offset.Row+i Offset.Column+j]==Display.Element), 
              set(Handles.Index.Elements(i,j),'string',Str, ...
                           'backgroundcolor',Color.Struct, ...
                           'buttondownfcn',['md_edit(''MDEDIT-INTERNAL'',[',num2str(i),' ',num2str(j),'])'], ...
                           'uicontextmenu',Handles.UIContextMenu.Main, ...
                           'enable','inactive');
            else,
              set(Handles.Index.Elements(i,j),'string',Str, ...
                           'backgroundcolor',TempColor, ...
                           'buttondownfcn',['md_edit(''MDEDIT-INTERNAL'',[',num2str(i),' ',num2str(j),'])'], ...
                           'uicontextmenu',Handles.UIContextMenu.Main, ...
                           'enable','inactive');
            end;
          else,
            set(Handles.Index.Elements(i,j),'string',Str, ...
                         'backgroundcolor',TempColor, ...
                         'buttondownfcn','', ...
                         'enable','on');
          end;
        else,
          set(Handles.Index.Elements(i,j),'string','', ...
                       'backgroundcolor',Color.Background, ...
                       'buttondownfcn','', ...
                       'enable','off');
        end;
      end;
    end;

    %%---------------------------------------------------------------------------
    %% Set the row labels.
    %%---------------------------------------------------------------------------
    for i=1:DisplayNumber.Rows,
      if isempty(Options.DimLabel),
        Str=num2str(Offset.Row+i);
      else,
        TempVar=Options.DimLabel{Display.Dimensions(1)};
        if ischar(TempVar),
          if strcmp(TempVar,'123'),
            Str=num2str(Offset.Row+i);
          elseif strcmp(TempVar,'ABC'),
            Str=LcColLabel(Offset.Row+i);
          elseif strcmp(TempVar,'none'),
            Str='';
          else,
            Str=num2str(Offset.Row+i);
          end;
        else, % cell
          if (Offset.Row+i)>length(TempVar),
            Str=num2str(Offset.Row+i);
          else,
            Str=TempVar{Offset.Row+i};
          end;
        end;
      end;
      if (i<=size(Display.Data,1)),
        if isempty(find(~(Selected.Rows-i-Offset.Row))),
          set(Handles.Index.Label.Row(i),'string',Str, ...
                     'style','togglebutton', ...
                     'value',0, ...
                     'buttondownfcn','', ...
                     'enable','on');
        else,
          set(Handles.Index.Label.Row(i),'string',Str, ...
                     'style','togglebutton', ...
                     'value',1, ...
                     'enable','on');
        end;
      else,
        set(Handles.Index.Label.Row(i),'string',Str, ...
                     'style','pushbutton', ...
                     'value',0, ...
                     'enable','off');
      end;
    end;

    %%-------------------------------------------------------------------------
    %% Set the column labels.
    %%-------------------------------------------------------------------------
    for j=1:DisplayNumber.Columns,
      if isempty(Options.DimLabel),
        Str=LcColLabel(Offset.Column+j);
      else,
        TempVar=Options.DimLabel{Display.Dimensions(2)};
        if ischar(TempVar),
          if strcmp(TempVar,'123'),
            Str=num2str(Offset.Column+j);
          elseif strcmp(TempVar,'ABC'),
            Str=LcColLabel(Offset.Column+j);
          elseif strcmp(TempVar,'none'),
            Str='';
          else,
            Str=LcColLabel(Offset.Column+j);
          end;
        else, % cell
          if (Offset.Column+j)>length(TempVar),
            Str=LcColLabel(Offset.Column+j);
          else,
            Str=TempVar{Offset.Column+j};
          end;
        end;
      end;
      if (j<=size(Display.Data,2)),
        if isempty(find(~(Selected.Columns-j-Offset.Column))),
          set(Handles.Index.Label.Column(j),'string',Str, ...
                     'style','togglebutton', ...
                     'value',0, ...
                     'enable','on');
        else,
          set(Handles.Index.Label.Column(j),'string',Str, ...
                     'style','togglebutton', ...
                     'value',1, ...
                     'enable','on');
        end;
      else,
        set(Handles.Index.Label.Column(j),'string',Str, ...
                     'style','pushbutton', ...
                     'value',0, ...
                     'enable','off');
      end;
    end;

    %%-------------------------------------------------------------------------
    %% Enable the row scrollbar when necessary.
    %%-------------------------------------------------------------------------
    if DisplayNumber.Rows<size(Display.Data,1),
      set(Handles.Index.Scrollbar.Row,'min',-size(Display.Data,1)+DisplayNumber.Rows, ...
              'max',0, ...
              'value',-Offset.Row, ...
              'enable','on');
    else,
      set(Handles.Index.Scrollbar.Row,'min',-1, ...
              'max',0, ...
              'value',0, ...
              'enable','off');
    end;

    %%-------------------------------------------------------------------------
    %% Enable the column scrollbar when necessary.
    %%-------------------------------------------------------------------------
    if DisplayNumber.Columns<size(Display.Data,2),
      set(Handles.Index.Scrollbar.Column,'max',size(Display.Data,2)-DisplayNumber.Columns, ...
              'min',0, ...
              'value',Offset.Column, ...
              'enable','on');
    else,
      set(Handles.Index.Scrollbar.Column,'min',0, ...
              'max',1, ...
              'value',0, ...
              'enable','off');
    end;

    %%-------------------------------------------------------------------------
    %% Restore the Display.Data from TempContent.
    %%-------------------------------------------------------------------------
    if ~isempty(Options.Multiple),
      Display.Data=TempContent;
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STRUCTURE and MULTIPLE OPTION specific elements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(Handles.Menu.Structure.Main,'enable','off');
    if isstruct(Display.Data) || ~isempty(Options.Multiple),

      %%-----------------------------------------------------------------------
      %% Set the field labels.
      %%-----------------------------------------------------------------------
      if isstruct(Display.Data),
        set(Handles.Menu.Structure.Main,'enable','on');
        if ~isempty(Display.Data),
          TempElm=Display.Data(Display.Element(1),Display.Element(2));
        else
          TempElm=Display.Data;
        end;
        TempFld=fieldnames(TempElm);
      else  % ~isempty(Options.Multiple) % MULTIPLE mode
        if ~isempty(Display.Data{1}),
          TempRef=Display.Index(end).Ref(end);
          TempRef.subs(Display.Dimensions)={Display.Element(1) Display.Element(2)};
          TempElm=LcSubsRef(Display.Data{1},TempRef);
        else
          TempElm=Display.Data{1};
        end;

        %%---------------------------------------------------------------------
        %% Check for struct datasets among the MULTIPLE datasets.
        %%---------------------------------------------------------------------
        for i=1:length(Display.Data),
          if strcmp(Options.Multiple,'123'),
            TempStr=num2str(i);
          elseif strcmp(Options.Multiple,'ABC'),
            TempStr=LcColLabel(i);
          else,
            TempStr=Options.Multiple{i};
          end;
          switch class(Display.Data{i}),
          case 'struct',
            FNames=fieldnames(Display.Data{i});
            for j=1:length(FNames),
              if isempty(Display.MultipleField),
                Display.MultipleField={i FNames{j}};
                TempFld={[TempStr '.' FNames{j}]};
              else
                Display.MultipleField(end+1,1:2)={i FNames{j}};
                TempFld(end+1,1)={[TempStr '.' FNames{j}]};
              end;
            end;
          otherwise,
            if isempty(Display.MultipleField),
              Display.MultipleField={i []};
              TempFld={TempStr};
            else
              Display.MultipleField(end+1,1:2)={i []};
              TempFld(end+1,1)={TempStr};
            end;
          end;
        end;
      end;

      %%-----------------------------------------------------------------------
      %% Load icons: NUMERIC, SPARSE, CHAR, CELL, STRUCT and other.
      %%-----------------------------------------------------------------------
      Icon=Icons;

      %%-----------------------------------------------------------------------
      %% Set the value controls.
      %%-----------------------------------------------------------------------
      for i=1:DisplayNumber.Fields,
        if Offset.Field+i<=length(TempFld),
          set(Handles.Field.Label.FieldName(i),'string',TempFld{Offset.Field+i});
          FieldCData=[];
          if isempty(TempElm),
            set(Handles.Field.Field(i),'string','', ...
                                'enable','off', ...
                                'backgroundcolor',Color.Background);
            set(Handles.Field.Label.FieldName(i),'enable','on');
          else

            %%-----------------------------------------------------------------
            %% Get field of structure to display when editing a structure.
            %%-----------------------------------------------------------------
            if isstruct(Display.Data),
              set(Handles.Field.Label.FieldName(i),'enable','on');

              TempVar=getfield(TempElm,TempFld{Offset.Field+i});
              Str=LcVal2Str(TempVar);
              switch class(TempVar),
              case {'double','uint8'},
                FieldCData=Icon.Matrix;
              case {'sparse'},
                FieldCData=Icon.Sparse;
              case {'char'},
                FieldCData=Icon.Char;
              case {'struct'},
                FieldCData=Icon.Struct;
              case {'cell'},
                FieldCData=Icon.Cell;
              otherwise, % userdefined class
                FieldCData=Icon.CustomClass;
              end;
              set(Handles.Field.Field(i),'string',Str, ...
                                  'buttondownfcn',['md_edit(''MDEDIT-INTERNAL'',[-6 ',num2str(i),'])'], ...
                                  'callback','', ...
                                  'enable','inactive', ...
                                  'backgroundcolor',Color.Edit);

            %%-----------------------------------------------------------------
            %% Editing in MULTIPLE mode, get value from appropriate dataset.
            %%-----------------------------------------------------------------
            else  % ~isstruct(Display.Data)
              set(Handles.Field.Label.FieldName(i),'enable','inactive');
              FieldEnable='inactive';

              %%---------------------------------------------------------------
              %% Get the appropriate data from the matrices or vectors.
              %%---------------------------------------------------------------
              if ischar(Display.MultipleField{i,2}), % part of a structure
                TempRef=Display.Index(end).Ref(end);
                TempRef.subs(Display.Dimensions)={Display.Element(1) Display.Element(2)};
                TempVar=LcSubsRef(Display.Data{Offset.Field+Display.MultipleField{i,1}},TempRef);
                TempVar=getfield(TempVar,Display.MultipleField{i,2});
                Str=LcVal2Str(TempVar);
              else
                if Options.MultipleLines,
                  if ischar(Display.Data{Offset.Field+Display.MultipleField{i,1}}),
                    TempVar=Display.Data{Offset.Field+Display.MultipleField{i,1}}(Display.Element(1),:);
                  else
                    TempVar=Display.Data{Offset.Field+Display.MultipleField{i,1}}(Display.Element(1));
                  end;
                else
                  TempRef=Display.Index(end).Ref(end);
                  TempRef.subs(Display.Dimensions)={Display.Element(1) Display.Element(2)};
                  TempVar=LcSubsRef(Display.Data{Offset.Field+Display.MultipleField{i,1}},TempRef);
                end;

                %%---------------------------------------------------------------
                %% Convert the data into a string.
                %%---------------------------------------------------------------
                if size(TempVar,2)>1, % 'char', % part of a STRING vector [i.e. CHAR matrix among vectors]
                  Str=TempVar;
                  FieldEnable='on';
                else  % size(TempVar)==[1 1]
                  switch class(TempVar),
                  case {'double','uint8','sparse'},
                    Str=sprintf('%g',double(TempVar));
                    FieldEnable='on';
                  case {'char'},
                    if abs(TempVar)<32 || abs(TempVar)>255,
                      Str=sprintf('#%i',abs(TempVar));
                      set(Handles.Menu.Char.Main,'enable','off');
                    else
                      Str=TempVar;
                    end;
                    FieldEnable='on';
                  case {'cell'},
                    TempVar=TempVar{1};
                    Str=LcVal2Str(TempVar);
                  otherwise, % clickable classes
                    Str=['[' class(TempVar) ']'];
                  end;
                end;
              end;

              %%---------------------------------------------------------------
              %% Set the appropriate icon.
              %%---------------------------------------------------------------
              switch class(TempVar),
              case {'double','uint8'},
                FieldCData=Icon.Matrix;
              case {'sparse'},
                FieldCData=Icon.Sparse;
              case {'char'},
                FieldCData=Icon.Char;
              case {'struct'},
                FieldCData=Icon.Struct;
              case {'cell'},
                FieldCData=Icon.Cell;
              otherwise, % userdefined class
                FieldCData=Icon.CustomClass;
              end;

              %%---------------------------------------------------------------
              %% Finally set the value control.
              %%----------------------------------------------------------------
              set(Handles.Field.Field(i),'string',Str, ...
                                  'buttondownfcn',['md_edit(''MDEDIT-INTERNAL'',[-6 ',num2str(i),'])'], ...
                                  'callback',['md_edit(''MDEDIT-INTERNAL'',[-6 ',num2str(i),'])'], ...
                                  'enable',FieldEnable, ...
                                  'backgroundcolor',Color.Edit);
            end;
          end;
          set(Handles.Field.Label.FieldType(i),'enable','inactive','string','','cdata',FieldCData);
        else
          set(Handles.Field.Label.FieldName(i),'string','','enable','off');
          set(Handles.Field.Label.FieldType(i),'enable','inactive','string','','cdata',[]);
          set(Handles.Field.Field(i),'string','','enable','inactive','backgroundcolor',Color.Background);
        end;
      end;

      %%-----------------------------------------------------------------------
      %% Enable the scrollbar when necessary.
      %%-----------------------------------------------------------------------
      if DisplayNumber.Fields<length(TempFld),
        set(Handles.Field.Scrollbar.Field,'max',0, ...
                'min',-length(TempFld)+DisplayNumber.Fields, ...
                'value',-Offset.Field, ...
                'enable','on');
      else
        set(Handles.Field.Scrollbar.Field,'min',0, ...
                'max',1, ...
                'value',0, ...
                'enable','off');
      end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% End of STRUCTURE and MULTIPLE OPTION specific elements
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% End of the 2D specific elements of the interface.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case Options.EditModes(2), % line

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Now update the 1D specific elements of the interface.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%-------------------------------------------------------------------------
    %% Set the main label of the strings.
    %%-------------------------------------------------------------------------
    if isempty(Options.DimElmLabel),
      Str=sprintf('dimension %i',Display.Dimensions(1));
    else
      TempVar=Options.DimElmLabel{Display.Dimensions(1)};
      if ischar(TempVar),
        Str=TempVar;
      else
        Str=sprintf('dimension %i',Display.Dimensions(1));
      end;
    end;
    set(Handles.LineEdit.Label.LineMain,'string',Str);

    %%-------------------------------------------------------------------------
    %% Set the string edit controls.
    %%-------------------------------------------------------------------------
    Offset.Row=min(Offset.Row,max(size(Display.Data,1)-(DisplayNumber.Rows+3),0));
    Offset.Column=0;
    for i=1:(DisplayNumber.Rows+3),
      if i<=size(Display.Data,1),
        if isempty(find(~(Selected.Rows-i-Offset.Row))),
          TempColor=Color.Edit;
        else
          TempColor=Color.Selected;
        end;
        Str=deblank(Display.Data(Offset.Row+i,:));
        set(Handles.LineEdit.Lines(i),'string',Str, ...
                    'backgroundcolor',TempColor, ...
                    'enable','on');
      else
        set(Handles.LineEdit.Lines(i),'string','', ...
                     'backgroundcolor',Color.Background, ...
                     'enable','off');
      end;
    end;

    %%-------------------------------------------------------------------------
    %% Set the string labels.
    %%-------------------------------------------------------------------------
    for i=1:(DisplayNumber.Rows+3),
      if isempty(Options.DimLabel),
        Str=num2str(Offset.Row+i);
      else
        TempVar=Options.DimLabel{Display.Dimensions(1)};
        if ischar(TempVar),
          if strcmp(TempVar,'123'),
            Str=num2str(Offset.Row+i);
          elseif strcmp(TempVar,'ABC'),
            Str=LcColLabel(Offset.Row+i);
          elseif strcmp(TempVar,'none'),
            Str='';
          else
            Str=num2str(Offset.Row+i);
          end;
        else  % cell
          if (Offset.Row+i)>length(TempVar),
            Str=num2str(Offset.Row+i);
          else
            Str=TempVar{Offset.Row+i};
          end;
        end;
      end;
      if (i<=size(Display.Data,1)),
        if isempty(find(~(Selected.Rows-i-Offset.Row))),
          set(Handles.LineEdit.Label.Line(i),'string',Str, ...
                     'style','togglebutton', ...
                     'value',0, ...
                     'enable','on');
        else
          set(Handles.LineEdit.Label.Line(i),'string',Str, ...
                     'style','togglebutton', ...
                     'value',1, ...
                     'enable','on');
        end;
      else
        set(Handles.LineEdit.Label.Line(i),'string',Str, ...
                     'style','pushbutton', ...
                     'enable','off');
      end;
    end;

    %%-------------------------------------------------------------------------
    %% Enable the string scrollbar when appropriate.
    %%-------------------------------------------------------------------------
    if (DisplayNumber.Rows+3)<size(Display.Data,1),
      set(Handles.LineEdit.Scrollbar.Line,'min',-size(Display.Data,1)+DisplayNumber.Rows+3, ...
              'max',0, ...
              'value',-Offset.Row, ...
              'enable','on');
    else
      set(Handles.LineEdit.Scrollbar.Line,'min',-1, ...
              'max',0, ...
              'value',0, ...
              'enable','off');
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% End of the 1D specific elements of the interface.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  case Options.EditModes(3), % multilline

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Now update the 0D specific elements of the interface.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Offset.Row=0;
    Offset.Column=0;
    set(Handles.TextEdit.Text,'string',Display.Data, ...
                'backgroundcolor',Color.Edit, ...
                'enable','on');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% End of the 0D specific elements of the interface.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end;

  %%---------------------------------------------------------------------------
  %% Show the field elements when appropriate.
  %%---------------------------------------------------------------------------
  if isstruct(Display.Data) || ~isempty(Options.Multiple),
    set(Handles.Field.All,'visible','on');
  else
    set(Handles.Field.All,'visible','off');
  end;

  %%---------------------------------------------------------------------------
  %% Enable the char menu when a char matrix is edited that can be edited
  %% in the other editing modes.
  %%---------------------------------------------------------------------------
  if ischar(Display.Data) & ndims(Display.Data)==2 & (isempty(Display.Data) | all(all((Display.Data>=32) & (Display.Data<=255)))),
    set(Handles.Menu.Char.Main,'enable','on');
  else
    set(Handles.Menu.Char.Main,'enable','off');
  end;

  %%---------------------------------------------------------------------------
  %% Show the appropriate elements.
  %%---------------------------------------------------------------------------
  for i=1:length(Options.EditModes),
    if strcmp(Options.EditModes{i},Options.EditMode),
      set(Handles.Index.All{i},'visible','on');
    else
      set(Handles.Index.All{i},'visible','off');
    end;
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% A DEFAULT CLASS OF MATLAB: NUMERIC, CHARACTER, CELL OR STRUCT
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

otherwise, % other (probably user defined) class

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% A NOT-SUPPORTED CLASS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % although struct() can be used to extract all data from
  % such a structure, converting it back is not possible
  % anymore. I.e. no editing posible. Possible extention:
  % call a user defined editor.

  MsgBoxStr=sprintf('Class ''%s'' not supported',TempClass);
  uiwait(msgbox(MsgBoxStr,'modal'));
  set(Handles.General.MessageWindow,'string',MsgBoxStr);
  if length(Display.Index(end).Ref)==0,
    LcStackUDF(Handles.Figure,'CommandStack',[-1 1],'top');
  else % Go back to previous level
    Display.Index(end)=[];
    LcStackUDF(Handles.Figure,'CommandStack',[NaN 0],'top');
  end;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% A NOT-SUPPORTED CLASS
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LcSeldim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sdim,indx]=LcSeldim(nsel,ndim,nsiz,init),

sdim=[];
if nargout>1,
  indx={};
end;

if nargin<3,
  fprintf(1,'* Not enough input arguments.\n');
  return;
elseif nargin==3,
  init={};
elseif nargin>4,
  fprintf(1,'* Too many input arguments.\n');
  return;
end;

if ndim<nsel,
  nsel=ndim;
end;

if nsel==1,
  dimens='dimension';
else
  dimens='dimensions';
end;

[fig,nd,rb]=LcUiSeldim(nsel,ndim,nsiz,init);

while length(sdim)~=nsel,
  set(fig,'userdata',[]);
  waitfor(fig,'userdata');
  switch get(fig,'userdata'),
  case 0,
    for dim=1:ndim,
      if get(rb(dim),'value')==1,
        init{dim}=':';
      else
        init{dim}=get(nd(dim),'string');
      end;
    end;
    ndim=ndim+1;
    nsiz(ndim)=1;
    init{ndim}=1;
    delete(fig);
    [fig,nd,rb]=LcUiSeldim(nsel,ndim,nsiz,init);
  case 1,
    for dim=1:ndim,
      if get(rb(dim),'value')==1,
        sdim=[sdim dim];
      end;
    end;
    if isempty(sdim),
    elseif length(sdim)~=nsel,
      Str=sprintf('Please select exactly %i %s.',nsel,dimens);
      uiwait(msgbox(Str,'modal'));
      sdim=[];
    else
      for dim=1:ndim,
        if all(dim~=sdim),
          TempVar=get(nd(dim),'string');
          [TempVal,TempCnt,TempErr,NextIndex] = sscanf(TempVar,'%i',1);
          if ~(TempCnt==1) || ~isempty(deblank(TempVar(NextIndex:end))),
            Str=sprintf('Specified index for dimension %i is not valid.',dim);
            uiwait(msgbox(Str,'modal'));
            sdim=[];
            break; % for dim=1:ndim,
          else
            set(nd(dim),'string',num2str(TempVal));
          end;
        end;
      end;
    end;
  end;
end;
if nargout>1,
  for dim=1:ndim,
    if all(dim~=sdim),
      indx{dim} = sscanf(get(nd(dim),'string'),'%i',1);
    else
      indx{dim}=':';
    end;
  end;
end;
delete(fig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LcUiSeldim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fig,nd,rb]=LcUiSeldim(nsel,ndim,nsiz,init)

LightGray=[1 1 1]*192/255;

if nsel==1,
  dimens='dimension';
else
  dimens='dimensions';
end;

Fig_Width=280;
Num_Field=40;
Field_Height=20;
Margin=20;
Fig_Height=Field_Height*(ndim+3);

ss = get(0,'ScreenSize');
swidth = ss(3);
sheight = ss(4);
left = (swidth-Fig_Width)/2;
bottom = (sheight-Fig_Height)/2;
rect = [left bottom Fig_Width Fig_Height];

fig=figure('menu','none', ...
        'units','pixels', ...
        'position',rect, ...
        'color',LightGray, ...
        'inverthardcopy','off', ...
        'closerequestfcn','', ...
        'resize','off', ...
        'numbertitle','off', ...
        'handlevisibility','off', ...
        'name','Select ...', ...
        'tag','Select Window');

for dim=1:ndim,
  Str=sprintf('dimension %i  (indices 1:%i)',dim,nsiz(dim));
  nd(dim)=uicontrol('style','edit', ...
            'position',[Fig_Width-Margin-Num_Field Fig_Height-(1+dim)*Field_Height Num_Field Field_Height], ...
            'parent',fig, ...
            'enable','on', ...
            'string','1');
  rb(dim)=uicontrol('style','radiobutton', ...
            'position',[Margin Fig_Height-(1+dim)*Field_Height Fig_Width-2*Margin-Num_Field Field_Height], ...
            'callback', 'if get(gcbo,''value''), set(get(gcbo,''userdata''),''visible'',''off''); else, set(get(gcbo,''userdata''),''visible'',''on''); end;', ...
            'parent',fig, ...
            'string',Str, ...
            'userdata',nd(dim));
  if ~isempty(init),
    if ischar(init{dim}),
      if strcmp(init{dim},':'), % init{dim}==':'
        set(nd(dim),'visible','off');
        set(rb(dim),'value',1);
      else
        set(nd(dim),'string',init{dim});
      end;
    else
      set(nd(dim),'string',num2str(init{dim}));
    end;
  end;
end;

Str=sprintf('Select %i %s and remaining indices ...',nsel,dimens);
uicontrol('style','text', ...
          'position',[Margin Fig_Height-Field_Height Fig_Width-2*Margin Field_Height], ...
          'string',Str, ...
          'parent',fig);
uicontrol('style','pushbutton', ...
          'position',[Margin Field_Height Fig_Width-2*Margin Field_Height], ...
          'string','add dimension', ...
          'parent',fig, ...
          'callback','set(gcbf,''userdata'',0)');
uicontrol('style','pushbutton', ...
          'position',[Margin 0 Fig_Width-2*Margin Field_Height], ...
          'string','continue', ...
          'parent',fig, ...
          'callback','set(gcbf,''userdata'',1)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INTERPRET OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Options=LcProcOpt(Display,DefaultOptions,TempOpt)

Options=DefaultOptions;

for opt=1:2:length(TempOpt),
  %%%%%%%%%%%%%%%%
  switch TempOpt{opt},
  case 'multiple', % matrix labels
  %%%%%%%%%%%%%%%%
    if ~iscell(Display.Data) || ndims(Display.Data)>2 || min(size(Display.Data))>1 || isempty(Display.Data),
      uiwait(msgbox('Option ''multiple'' only valid for 1-D cell arrays.','modal'));
    elseif max(size(Display.Data))==1, % one element
      uiwait(msgbox('Option ''multiple'' does not match single input.','modal'));
    else
      TempVar=size(Display.Data{1});
      NonCharSize=[];
      CharSize=[];
      Options.Multiple=1;
      Options.MultipleLines=0;
      for matr=1:length(Display.Data),
        if ischar(Display.Data{matr}),
          if isempty(CharSize),
            CharSize=size(Display.Data{matr});
          elseif ~all(CharSize==size(Display.Data{matr})),
            if CharSize(1)~=size(Display.Data{matr},1),
              Options.Multiple=0;
            else
              Options.MultipleLines=1;
            end;
          end;
        else
          if isempty(NonCharSize),
            NonCharSize=size(Display.Data{matr});
          elseif ~all(NonCharSize==size(Display.Data{matr})),
            Options.Multiple=0;
          end;
        end;
      end;
      if ~isempty(CharSize) && ~isempty(NonCharSize),
        if ~isequal(CharSize,NonCharSize)
          if (NonCharSize(2)~=1) || (CharSize(1)~=NonCharSize(1)),
            Options.Multiple=0;
          elseif (NonCharSize(2)==1) && (CharSize(1)==NonCharSize(1)),
            Options.MultipleLines=1;
          end;
        end;
      end;
      if Options.MultipleLines,
        for matr=1:length(Display.Data),
          if ischar(Display.Data{matr})
            if (~isempty(Display.Data{matr}) && ~all(all((Display.Data{matr}>=32 & Display.Data{matr}<=255) | Display.Data{matr}==0))),
              Options.Multiple=0;
            end;
          end;
        end;
      end;
      if Options.Multiple,
        Options.Multiple={};
        TempVar=TempOpt{opt+1};
        if ischar(TempVar),
          if strcmp(TempVar,'123') || strcmp(TempVar,'ABC'),
            Options.Multiple=TempVar;
          else
            if size(TempVar,1)==length(Display.Data),
              Options.Multiple=cellstr(TempVar);
            else
              uiwait(msgbox('Number of names does not equal number of matrices.','Option ''multiple''','modal'));
            end;
          end;
        elseif iscell(TempVar),
          if length(TempVar)==length(Display.Data),
            Options.Multiple=TempVar;
          else
            uiwait(msgbox('Number of names does not equal number of matrices.','Option ''multiple''','modal'));
          end;
        else
          uiwait(msgbox('Value should be a cell of appropriate length or a string.','Option ''multiple''','modal'));
        end
      else
        Options.Multiple={};
        uiwait(msgbox('Matrix dimensions do not match.','Option ''multiple''','modal'));
      end;
    end;
  %%%%%%%%%%%%%%%%
  case 'dimlabel', % dimension labels
  %%%%%%%%%%%%%%%%
    TempVar=TempOpt{opt+1};
    if ischar(TempVar),
      if size(TempVar,1)==ndims(Display.Data),
        Options.DimElmLabel=cellstr(TempVar);
      else
        uiwait(msgbox('Number of names does not equal number of dimensions.','Option ''dimlabel''','modal'));
      end;
    elseif iscell(TempVar),
      if length(TempVar)==ndims(Display.Data),
        Options.DimElmLabel=TempVar;
      else
        uiwait(msgbox('Number of names does not equal number of dimensions.','Option ''dimlabel''','modal'));
      end;
    else
      uiwait(msgbox('Value should be a cell of appropriate length or a string.','Option ''dimlabel''','modal'));
    end;
  %%%%%%%%%%%%%%%%
  case 'dimfixed', % dimensions fixed indicators
  %%%%%%%%%%%%%%%%
    TempVar=TempOpt{opt+1};
    if isnumeric(TempVar) && ndims(TempVar)==2 && min(size(TempVar,1))==1 && length(TempVar)==ndims(Display.Data),
      Options.DimFixed=TempVar;
    else
      uiwait(msgbox('Value should be a numeric array of length equal to the number of dimensions.','Option ''dimfixed''','modal'));
    end;
  %%%%%%%%%%%%%%%%
  case 'label', % dimension indices labels
  %%%%%%%%%%%%%%%%
    TempVar=TempOpt{opt+1};
    if iscell(TempVar),
      if length(TempVar)==ndims(Display.Data),
        Options.DimLabel=TempVar;
      else
        uiwait(msgbox('Number of names does not equal number of dimensions.','Option ''label''','modal'));
      end;
    else
      uiwait(msgbox('Value should be a 1-D cell array.','Option ''label''','modal'));
    end;
  %%%%%%%%%%%%%%%%
  case 'specmode', % char edit mode
  %%%%%%%%%%%%%%%%
    if ischar(Display.Data) && ndims(Display.Data)==2,
      TempVar=TempOpt{opt+1};
      switch TempVar,
      case Options.EditModes,
        Options.EditMode=Options.EditModes{1};
      otherwise,
        uiwait(msgbox('Valid values are: ''array'', ''line'', ''multiline''.','Option ''specmode''','modal'));
      end;
    else,
      uiwait(msgbox('Option ''specmode'' valid only for 2-D char arrays.','Option ''specmode''','modal'));
    end;
  %%%%%%%%%%%%%%%%
  otherwise,
  %%%%%%%%%%%%%%%%
    fprintf(1,'Unknown option:');
    disp(TempOpt{opt});
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Convert a value into a string.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Str=LcVal2Str(Val),
if ndims(Val)>2 | min(size(Val))>1,
  Str=['[' num2str(ndims(Val)) 'D ' class(Val) ']'];
elseif min(size(Val))==1,
  if max(size(Val))>1,
    if ischar(Val) & size(Val,1)==1,
      Str=[Val ' [char]'];
    else,
      Str=['[1D ' class(Val) ']'];
    end;
  else, % size(Val)==[1 1]
    switch class(Val),
    case {'double','uint8','sparse'},
      Str=[num2str(double(Val)) ' [' class(Val) ']'];
    case {'char'},
      Str=[Val ' [char]'];
    case {'struct'},
      Str='[struct]';
    case {'cell'},
      Str='[cell]';
    otherwise, % userdefined class
      Str=['[' class(Val) ']'];
    end;
  end;
else, % min(size(Val))==0,
  Str=['[empty ' class(Val) ']'];
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LcBorder3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function l=LcBorder3D(x,y,dx,dy,varargin),
HShift=-1;
VShift=-1;
DarkGray=[1 1 1]*128/255;
White=[1 1 1];
Black=[0 0 0];
LightGray=[1 1 1]*223/255;
MidGray=[1 1 1]*192/255;
l(5)=patch(HShift+[x+dx-1 x x x+dx-1],VShift+[y+dy-1 y+dy-1 y y],-[1 1 1 1],1,'facecolor',MidGray,'edgecolor','none',varargin{:});
l(1)=line(HShift+[x+dx-1 x x],VShift+[y+dy-1 y+dy-1 y],-0.5*[1 1 1],'color',White,varargin{:});
l(2)=line(HShift+[x+dx-2 x+1 x+1],VShift+[y+dy-2 y+dy-2 y+1],-0.5*[1 1 1],'color',LightGray,varargin{:});
l(3)=line(HShift+[x+1 x+dx-2 x+dx-2],VShift+[y+1 y+1 y+dy-2],-0.5*[1 1 1],'color',DarkGray,varargin{:});
l(4)=line(HShift+[x x+dx-1 x+dx-1],VShift+[y y y+dy-1],-0.5*[1 1 1],'color',Black,varargin{:});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LcHLine3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function l=LcHLine3D(x,y,dx,varargin),
DarkGray=[.502 .502 .502];
White=[1 1 1];
l(1)=line([x x+dx],[y y],[0 0],'color',DarkGray,varargin{:});
l(2)=line([x x+dx],[y-1 y-1],[0 0],'color',White,varargin{:});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LcSubsRef - like subsref but dealing with empty subscript references
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B=LcSubsRef(A,S),
if isempty(S),
  B=A;
else,
  B=subsref(A,S);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LcSubsAsgn - like subsasgn but dealing with empty subscript references
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Aout=LcSubsAsgn(Ain,S,B),
if isempty(S),
  Aout=B;
else
  Aout=subsasgn(Ain,S,B);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LcRef2Str - creates a string from a reference list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str=LcRef2Str(ref)
% REF2STR creates a string from a reference list
%
%     See also SUBSINDEX

%     Copyright (c)  H.R.A. Jagers  12-05-1996

if nargin>1,
  fprintf(1,' * Too many input arguments\n');
elseif nargin==1,
  if isempty(ref) || isstruct(ref),
    str='';
    for k=1:length(ref),
      if strcmp(ref(k).type,'.'),
        str=[str '.' ref(k).subs];
      else % ref(k).type equals '()' or '{}'
        str=[str ref(k).type(1)];
        for  l=1:length(ref(k).subs),
          if l~=1,
            str=[str ','];
          end;
          if ischar(ref(k).subs{l}), % catch ':'
            str=[str sprintf('%s',ref(k).subs{l})];
          else
            str=[str sprintf('%i',ref(k).subs{l})];
          end;
        end;
        str=[str ref(k).type(2)];      
      end;
    end;
  else
    fprintf(1,' * Expected a reference list as input.\n');
  end;
else
  fprintf(1,' * Too few input arguments\n');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LcColLabel - generates the label of a column as used by spreadsheets, i.e. A-Z,AA-AZ,BA-BZ,...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function label=LcColLabel(i);
% COLLABEL generates the label of a column
%      as used by spreadsheets, i.e. A-Z,AA-AZ,BA-BZ, etc.

% Copyright (c) H.R.A. Jagers

if nargin~=1,
  fprintf(1,'* exactly one input argument expected');
  return;
end;

% if i is a matrix it must be a column vector
if min(size(i,2),size(i,1))~=1,
  fprintf(1,'* input argument must be a scalar or vector');
  return;
end;

if (size(i,2)>size(i,1)),
  i=i';
end;

label=[];
first=1;
while any(i~=0),
  t=fix((i-1)/26);
  j=i-26*t;
  j=j+32+(i~=0)*32;
  i=t;
  label=[j label];
  first=0;
end;
label=char(label);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LcGetUDF - Get the value of the field in a userdata structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Value=LcGetUDF(Handle,FieldName),
%GETUDF Get the value of the field in a userdata structure.
%
%    V = GETUDF(H,'PropertyName') returns the value of the specified
%    field 'PropertyName' in the structure stored in the UserData
%    property of the graphics object with handle H. If H is a
%    vector of handles, then LcGetUDF will return an M-by-1 cell array
%    of values where M is equal to length(H).
%
%    See also SETUDF, ISUDF, RMUDF, WAITFORUDF, GET.

% Copyright (c) 1999, H.R.A. Jagers, WL | delft hydraulics, The Netherlands

Value=cell(length(Handle(:)),1);
for H=1:length(Handle(:)),
  UserData=get(Handle(H),'userdata');
  if isfield(UserData,FieldName),
    Value{H}=getfield(UserData,FieldName);
  else,
    error(['invalid property: ' FieldName '.']);
    Value{H}=[];
  end;
end;

if length(Handle)==1,
  Value=Value{1};
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LcLcSetUDF - Set the value of the field in a userdata structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LcSetUDF(Handle,FieldName,Value),
% SETUDF Set field in a userdata structure.
%
%    SETUDF(H,'PropertyName',PropertyValue) sets the value of
%    the specified field 'PropertyName' in the structure stored
%    in the UserData property of the graphics object with handle H.
%    H can be a vector of handles, in which case SETUDF sets the
%    properties' values for all the objects.
%
%    See also GETUDF, ISUDF, RMUDF, WAITFORUDF, SET.

% Copyright (c) 1999, H.R.A. Jagers, WL | delft hydraulics, The Netherlands

for H=1:length(Handle(:)),
  UserData=get(Handle(H),'userdata');
  if ~isstruct(UserData),
    if ~isempty(UserData),
      error('nonstructure userdata.');
    end;
    UserData=[];
  end;
  UserData=setfield(UserData,FieldName,Value);
  set(Handle(H),'userdata',UserData);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LcStackUDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LcStackUDF(Handle,StackName,Value,Opt),

UserData=get(Handle,'userdata');
if ~isstruct(UserData),
  if ~isempty(UserData),
    uiwait(msgbox(['Overwriting nonstructure userdata.'],'modal'));
  end;
  UserData=[];
  Stack={Value};
elseif isfield(UserData,StackName),
  Stack=getfield(UserData,StackName);
  if (nargin==4) && strcmp(Opt,'top'),
    Stack={Value,Stack{:}};
  else,
    Stack{length(Stack)+1}=Value;
  end;
else,
  Stack={Value};
end;
UserData=setfield(UserData,StackName,Stack);
set(Handle,'userdata',UserData);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LcWaitForUDF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LcWaitForUDF(Handle,FieldName,WaitforValue),
% WAITFORUDF Block execution and wait for change of UserData field.
%
%    WAITFORUDF(H,'PropertyName'), returns when the value of
%    'PropertyName' for the graphics object handle changes.
%    If 'PropertyName' is not a valid field of the structure
%    stored in the userdata property of the object, or if the
%    handle does not exist, waitfor returns immediately without
%    processing any events.
% 
%    WAITFORUDF(H,'PropertyName',PropertyValue), returns when
%    the value of the 'PropertyName' field of the userdata of
%    the graphics object handle changes to PropertyValue. If
%    'PropertyName' is set to PropertyValue, if 'PropertyName'
%    is not a valid field of the structure stored in the userdata
%    property of the object, or if the handle does not exist,
%    waitfor returns immediately without processing any events.
% 
%    While LcWaitForUDF blocks an execution stream, it processes
%    events as would drawnow, allowing callbacks to execute. Nested
%    calls to LcWaitForUDF are supported, and earlier calls to
%    LcWaitForUDF will not return until all later calls have returned,
%    even if the condition upon which the earlier call is blocking
%    has been met.
% 
%    See also SETUDF, WAITFOR.

% Copyright (c) 1999, H.R.A. Jagers, WL | delft hydraulics, The Netherlands

if nargin<2,
  error('Not enough input arguments.');
end;
if ~ishandle(Handle),
  return;
end;
UserData=get(Handle,'userdata');
if ~isfield(UserData,FieldName),
  return;
else,
  StartValue=getfield(UserData,FieldName);
  if (nargin==3) && isequal(StartValue,WaitforValue),
    return;
  end;
  while 1,
    waitfor(Handle,'userdata'); % waitfor a change of the userdata field
    if ~ishandle(Handle), % graphics object deleted
      return;
    end;
    UserData=get(Handle,'userdata');
    if ~isfield(UserData,FieldName), % field deleted
      return;
    else,
      CurrentValue=getfield(UserData,FieldName);
      if (nargin==2), % no WaitforValue specified
        if ~isequal(CurrentValue,StartValue), % changed
          return;
        end;
      elseif isequal(CurrentValue,WaitforValue), % changed and equal to WaitforValue
        return;
      end;
    end;
  end;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Icons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Icon=Icons,
Icon.CustomClass=[131071 130111 126991 124903 122867 ...
122651 113677 113133 111597 112605 112413 110637 119787 ...
118747 124831 129151 131071];
Icon.Cell=[131071 65793 65793 80185 80185 80185 65793 ...
65793 131071 65793 65793 80185 80185 80185 65793 65793 ...
131071];
Icon.Struct=[131071 130687 130687 130687 130687 130687 ...
98305 98305 104857 104857 104857 104857 104857 104857 ...
104857 104857 131071];
Icon.Char=[130879 130591 130783 130783 130591 131039 ...
122911 130783 130783 130591 131071 130879 130591 130783 ...
130783 130783 131071];
Icon.Sparse=[131071 131071 100351 100351 100351 100351 ...
131071 130111 130111 130111 130111 131071 131041 131041 ...
131041 131041 131071];
Icon.Matrix=[131071 131071 99361 99361 99361 99361 ...
131071 99361 99361 99361 99361 131071 99361 99361 99361 ...
99361 131071];
nms=fieldnames(Icon);
Vals=uint8([0 240]);
for i=1:length(nms)
  tmp=abs(dec2base(getfield(Icon,nms{i}),2))'-47;
  Icon=setfield(Icon,nms{i},repmat(Vals(tmp),[1 1 3]));
end
