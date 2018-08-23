function Res=plot_int1(Hf,Type,WaitFor)
%------------------------------------------------------------------------------
% plot_int1 function                                                  plotting
% Description: Given a plot wait for the use to click a keyboard key or
%              a mouse click on the plot, and return the key and click
%              position.
% Input  : - Handle for current figure in which the data points are ploted.
%            If empty matrix, then use gca. Default is [].
%          - Type of key to wait for:
%            'key'   - wait for a keyboard click on the plot (default).
%            'mouse' - wait for a single mouse click on the plot.
%            'rect'  - wait for a mouse selection of a rectanguar region.
%            'mousem'- wait for a multiple left mouse clicks, and use
%                      right click to abort.
%          - Waitfor action {'y'|'n'}. Will return only after the user
%            clicked a key/mouse. Default is 'y'.
% Output : - Result structure containing the following fields.
%            .Key    - Keyborad key entered.
%            .Pos    - Mouse position [X,Y] or rectangule position
%                      [xmin ymin width height].
%            .MB     - Mouse botton (1-left; 2-middle; 3-right).
%                      Note that if Type is 'mousem', the last .MB value
%                      will be always "3" corresponding to the exit click.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    Jan 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Res = plot_int1(gcf,'mouse');
%          Res = plot_int([],'key');
% Reliable: 2
%------------------------------------------------------------------------------
global Res

Def.Hf      = [];
Def.Type    = 'key';
Def.WaitFor = 'y';
if (nargin==0),
    Hf      = Def.Hf;
    Type    = Def.Type;
    WaitFor = Def.WaitFor;
elseif (nargin==1),
    Type    = Def.Type;
    WaitFor = Def.WaitFor;
elseif (nargin==2),
    WaitFor = Def.WaitFor;
elseif (nargin==3),
    % do nothing
else
    error('Illegal number of input arguments');
end

if (isempty(Hf)),
    Hf = gca;
end

Res = [];
if (ischar(Hf)),
    % get gcf from gca
    GCF = get(gca,'Parent');
    
    CC  = get(GCF,'CurrentCharacter');
    %Pos = get(GCF,'CurrentPoint');

    Res.Key = CC;
    Res.Pos = [];
    Res.MB  = [];

    set(gcf,'WindowButtonDownFcn','',...
 	                'KeyPressFcn','');
    %fprintf('Interactive mode terminated\n');
else
    switch lower(Type)
        case 'key'
            %fprintf('Interactive mode\n');
            set(Hf,'WindowButtonDownFcn','',...
                   'KeyPressFcn','[Res]=plot_int1(''key_press'');');
        case 'mouse'
            [X,Y,MB] = ginput(1);
            Res.Key  = [];
            Res.Pos  = [X,Y];
            Res.MB   = MB;
        case 'rect'
            Rect     = getrect(Hf);
            Res.Key  = [];
            Res.Pos  = Rect;  % [xmin ymin width height]
            Res.MB   = [];
        case 'mousem'
            [X,Y,MB] = ginput(1);
            Res.Key  = [];
            Res.Pos  = [X,Y];
            Res.MB   = MB;
            while (MB==1),
                [X,Y,MB] = ginput(1);
                Res.Pos  = [Res.Pos; [X,Y]];
                Res.MB   = [Res.MB; MB];
            end
        otherwise
    end

    switch lower(WaitFor)
       case 'y'
           waitfor(Hf, 'KeyPressFcn', '');
       otherwise
           % do nothing
    end

end
       


