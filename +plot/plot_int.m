function Res=plot_int(Hf,IndRM,CallFun,FunPar,varargin)
%------------------------------------------------------------------------------
% plot_int function                                                   plotting
% Description: Plot a 2-dimensional graph and interactively remove and
%              restore data points from the plot.
% Input  : - Handle for current figure in which the data points are ploted.
%            If empty matrix, then use gca. Default is [].
%            The program will read the X and Y of the data points from
%            the figure handle.
%            Alternatively this can be a cell array of the data to plot.
%            The cell array should contain the following information:
%            {X,Y,PlotPars,GCApar,XLabel,YLabel}
%            where X and Y are vectors of data to plot.
%            PlotPars: is a cell array of additional parameters to pass to
%            the plot function (e.g., {'^','MarkerSize',12}, or a string
%            containing the marker type. Default is 'ko'.
%            GCApar: is a cell array of keyword, value,... to pass to
%            the set(gca,...) command. Default is {}.
%            XLabel and YLabel are strings containing the X and Y labels.      
%            Default is ''.
%            Alteranatively this can be a character for internal recursive
%            use.
%          - Indices of data points to mark as deleted when the data
%            is first presented on the screen. Default is empty matrix.
%            Alternativey if this is the string 'all' then all the points
%            are marked as deleted.
%            This is useful if you like to mark some or all the points
%            as deleted and let the user return the points to the sample.
%          - Optional function to call when the 'f' option is
%            selected (see menu). The Function has the form:
%            [...] = Fun(X,Y,FunPar{:}); where X and Y are the current
%            X and Y which are not deleted and not marked by red cross.
%            Default is empty matrix (i.e., no function).
%          - Cell array of parameters to pass to the function (FunPar).
%          * Pairs of ...,key,val,... input arguments.
%            The following keywords are avialble:
%            'FunParInd' - The index of the parameter in FunPar that
%                          can be modified by typing 'g' on the plot.
%                          If not given, or empty, then this parameter
%                          will be ignored. Default is []. For example,
%                          if CallFun is 'fitgenpoly' then setting this
%                          parameter to 2 will enable the user to modify
%                          the degree of the polynomial fit.
% 	         'FunBehav'  - One of the following behaviours:
%                          'c' - run the function only when the user
%                                click 'f'.
%                          'i' - run the function immidetly after each
%                                user click with the exception of
%                                'q' (default).
%            'DispFit'   - Display best info when calling the function
%                          {'y'|'n'}. This option works only when the
%                          call function is fitgenpoly.m. Default is 'n'.
% Output : - Result structure containing the following fields.
%            .OrigX   - Original X in the list, including artificial data.
%            .OrigY   - Original Y in the list, including artificial data.
%            .X       - Final X left in the list.
%            .Y       - Final Y left in the list.
%            .Ind     - Indices of the final X/Y left in the list.
%            .IndRM   - Indices of the X/Y taken out of the list.
%            .FlagAll - A flag indicating if the point was
%                       deleted (0) or not (1).
%            .CallFun - The call function.
%            .FunPar  - Cell array of additional parameters passed to
%                       the call function.
%            .FunOut  - Cell array of arbitrary number of output parameters
%                       from the last function call.
% Tested : Matlab 7.0
%     By : Eran O. Ofek                    May 2005
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% See also: plot_int1.m
% Example: X = [1:1:20].';  Y=1+2.*X+0.5.*X.^2+randn(size(X)).*0.5;
%          Res=plot_int({X,Y,'bo',{},'X','Y'},[],'polyfit',{2});
%          Res=plot_int({X,Y,'bo',{},'X','Y'},[],'fitgenpoly',{1,2,'Plot','fitonly'},...
%                        'FunParInd',2,'FunBehav','i');
%          Res=plot_int({X,Y,'bo',{},'X','Y'},[],'fitgenpoly',...
%               {1,2,'Plot','fitonly','Xplot',[0:1:25]'},'FunParInd',2,'FunBehav','i');
%          % or
%          plot(X,Y,'bo');
%          Res=plot_int([],[],'fitgenpoly',{1,2,'Plot','fitonly','Xplot',[0:1:25]'},...
%                        'FunParInd',2,'FunBehav','i');
%          % or
%          plot(X,Y,'bo');
%          Res=plot_int([],'all','fitgenpoly',{1,2,'Plot','fitonly','Xplot',[0:1:25]'},...
%                        'FunParInd',2,'FunBehav','i');
%          % in order to wait until the function will return use:
%          waitfor(gcf,'KeyPressFcn','');
% Reliable: 2
%------------------------------------------------------------------------------
import Util.Geom.*

CrossedMarkerSize = 10;


Def.Hf        = [];
Def.IndRM     = [];
Def.CallFun   = [];
Def.FunPar    = {};
if (nargin==0),
   Hf        = Def.Hf;
   IndRM     = Def.IndRM;
   CallFun   = Def.CallFun;
   FunPar    = Def.FunPar;
elseif (nargin==1),
   IndRM     = Def.IndRM;
   CallFun   = Def.CallFun;
   FunPar    = Def.FunPar;
elseif (nargin==2),
   CallFun   = Def.CallFun;
   FunPar    = Def.FunPar;
elseif (nargin==3),
   FunPar    = Def.FunPar;
else
   % do nothing
end

DefV.FunParInd = [];
DefV.FunBehav  = 'i';
DefV.DispFit   = 'n';
InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

if (~ischar(Hf)),
   % not a recursive call
   if (iscell(Hf)),
      % user sent data to plot

      Def.PlotPars = 'ko';
      Def.GCApar   = {};
      Def.XLabel   = '';
      Def.YLabel   = '';

      switch length(Hf)
       case 2
          PlotPars = Def.PlotPars;
          GCApar   = Def.GCApar;
          XLabel   = Def.XLabel;
          YLabel   = Def.YLabel;
       case 3
          GCApar   = Def.GCApar;
          XLabel   = Def.XLabel;
          YLabel   = Def.YLabel;
       case 4
          XLabel   = Def.XLabel;
          YLabel   = Def.YLabel;
       case 5
          YLabel   = Def.YLabel;
       case 6
          % do nothing
       otherwise
          error('Illegal number of input arguments in first input cell array');
      end
      
%      GUIpar.GraphicH = [];
      X        = Hf{1};
      Y        = Hf{2};
      PlotPars = Hf{3};
      GCApar   = Hf{4};
      XLabel   = Hf{5};
      YLabel   = Hf{6};
      if (~iscell(PlotPars)),
         PlotPars = {PlotPars};
      end
      H2  = plot(X,Y,PlotPars{:});
      hold on;
      H1  = plot(X,Y,'rx','MarkerSize',CrossedMarkerSize);
      set(H1,'XData',[]);
      set(H1,'YData',[]);

      Hxl = xlabel(XLabel);     
      set(Hxl,'FontSize',16);
      Hyl = ylabel(YLabel);
      set(Hyl,'FontSize',16);
      if (isempty(GCApar)),
          % do nothing
      else
         set(gca,GCApar{:}); 
      end
   else

      if (isempty(Hf)),
         % work on existing plot - assuming only one set of symbols
         Hs = get(gca,'Children');
          
         if (length(Hs)>1),
             error('Existing plot have more than one sets of symbols');
         else
             % assume length(Hs)==1
             H2 = Hs(1);
             %H1 = [];

             X = get(H2,'XData');
             Y = get(H2,'YData');
             hold on;
             H1  = plot(X,Y,'rx','MarkerSize',CrossedMarkerSize);
             set(H1,'XData',[]);
             set(H1,'YData',[]);
         end
      else
         if (isnumeric(Hf)),
            % assuming Hf contains the figure handle
	    figure(Hf);
            Hs = get(gca,'Children');

            % work on existing plot - assuming only one set of symbols
            if (length(Hs)>1),
               error('Existing plot have more than one sets of symbols');
            else
               % assume length(Hs)==1
               H2 = Hs(1);
               %H1 = [];

               X = get(H2,'XData');
               Y = get(H2,'YData');
               hold on;
               H1  = plot(X,Y,'rx','MarkerSize',CrossedMarkerSize);
               set(H1,'XData',[]);
               set(H1,'YData',[]);
            end
         else
            error('Illegal type for first input parameter (Hf)');
         end
      end
   end


   HoldStatus = get(gcf,'NextPlot');
   hold on;
   
   X  = get(H2,'XData');
   Y  = get(H2,'YData');
   %if (isempty(X)),
   %   X  = get(H1,'XData');
   %   Y  = get(H1,'YData');
   %end
   N  = length(X);

   GUIpar.GraphicH = [];

   % store information in UserData
   Res.OrigX     = X.';
   Res.OrigY     = Y.';
   Res.X         = X.';
   Res.Y         = Y.';
   Res.Ind       = (1:1:N).';
   Res.IndRM     = IndRM;
   Res.FlagAll   = ones(N,1);   % 1 in; 0 crossed
   Res.CallFun   = CallFun;
   Res.FunPar    = FunPar;
   Res.FunParInd = InPar.FunParInd;
   Res.FunBehav  = InPar.FunBehav;
   Res.DispFit   = InPar.DispFit;
   Res.FunOut    = [];

   if (ischar(IndRM)),
      switch lower(IndRM)
       case 'all'
           Res.Ind = [];
           Res.IndRM = (1:1:N).';
           Res.FlagAll = zeros(N,1);
           set(H2,'XData',[]);
           set(H2,'YData',[]);
           H1 = plot(Res.X(Res.FlagAll==0),Res.Y(Res.FlagAll==0),'rx','MarkerSize',CrossedMarkerSize);

       otherwise
	      error('Unknwon IndRM option');
      end
   end

   % initilalization
   GUIpar.Res  = Res;
   
   GUIpar.H1 = H1;
   GUIpar.H2 = H2;
   set(gcf,'UserData',GUIpar);

   fprintf('GUI initilalization\n');
   set(gcf,'WindowButtonDownFcn','',...
           'KeyPressFcn','[Res]=plot_int(''key_press'');');
   plot_rm_menu;

   switch lower(GUIpar.Res.FunBehav)
    case 'i'
        %Res=plot_int('f');
        Res=plot_int('f');

        GUIpar.Res  = Res;
        set(gcf,'UserData',GUIpar);

   end

end


if (ischar(Hf)),
    % a recursive call

    GUIpar = get(gcf,'UserData');
    if (isempty(GUIpar)),
        error('Illegal call to function - first argument is a char and UserData not defined');
    end


    switch lower(Hf)
     case 'key_press'    
        Hf = get(gcf,'CurrentCharacter');
     otherwise
        % do nothing - user input parameter
    end
    
    switch lower(Hf)
        case {'donothing','key_press',''}
            % do nothing
        otherwise
            GUIpar = get(gcf,'UserData');
            
            if (isempty(GUIpar)),
                error('Illegal call to function - first argument is a char and UserData not defined');
            end
    end
    
    GUIpar.Exit = 0;

    switch lower(Hf)
        case {'donothing','key_press',''}
            % do nothing
        case {'q','off'}
            fprintf('GUI termination\n');
            set(gcf,'WindowButtonDownFcn','',...
 	                'KeyPressFcn','');
            fprintf('Interactive mode terminated\n');

            GUIpar.Exit = 1;
        case 'x'
           fprintf('Delete a point nearest to the mouse click - mark with red cross\n')
           [Xm,Ym] = ginput(1);
           Dist = plane_dist(Xm,Ym, GUIpar.Res.X, GUIpar.Res.Y);
           [Min,MinInd] = min(Dist);
           Dist = plane_dist(GUIpar.Res.OrigX,GUIpar.Res.OrigY,GUIpar.Res.X(MinInd),GUIpar.Res.Y(MinInd));
           [Min,MinInd] = min(Dist);

           GUIpar.Res.IndRM = unique([GUIpar.Res.IndRM; MinInd]);
           GUIpar.Res.FlagAll(GUIpar.Res.IndRM) = 0;

           set(GUIpar.H2,'XData',GUIpar.Res.OrigX(GUIpar.Res.FlagAll==1),...
                         'YData',GUIpar.Res.OrigY(GUIpar.Res.FlagAll==1));

           % red cross
           set(GUIpar.H1,'XData',GUIpar.Res.OrigX(GUIpar.Res.FlagAll==0),...
                         'YData',GUIpar.Res.OrigY(GUIpar.Res.FlagAll==0));
           
            switch lower(GUIpar.Res.FunBehav)
             case 'i'
                GUIpar = call_pi_fun(GUIpar);
            end           
         case 'r'
            fprintf('Delete all points within a rectangular region - mark points with red cross\n')
            Rect = getrect;
            Xmin = Rect(1); Xmax = Rect(1)+Rect(3);
            Ymin = Rect(2); Ymax = Rect(2)+Rect(4);

            MinInd   = find(GUIpar.Res.OrigX>Xmin & GUIpar.Res.OrigX<Xmax & ...
                            GUIpar.Res.OrigY>Ymin & GUIpar.Res.OrigY<Ymax & GUIpar.Res.FlagAll==1);

            GUIpar.Res.IndRM = unique([GUIpar.Res.IndRM; MinInd]);
            GUIpar.Res.FlagAll(GUIpar.Res.IndRM) = 0;

            set(GUIpar.H2,'XData',GUIpar.Res.OrigX(GUIpar.Res.FlagAll==1),...
                          'YData',GUIpar.Res.OrigY(GUIpar.Res.FlagAll==1));

            % red cros
            set(GUIpar.H1,'XData',GUIpar.Res.OrigX(GUIpar.Res.FlagAll==0),...
                          'YData',GUIpar.Res.OrigY(GUIpar.Res.FlagAll==0));
            %H1 = [H1; plot(Res.X(Res.FlagAll==0),Res.Y(Res.FlagAll==0),'rx','MarkerSize',CrossedMarkerSize)];

            switch lower(GUIpar.Res.FunBehav)
             case 'i'
                GUIpar = call_pi_fun(GUIpar);
            end

         case 'd'
            % delete all points
            fprintf('Delete all points - mark all points with red cross\n');
            
            N = length(GUIpar.Res.OrigX);
            GUIpar.Res.IndRM = [1:1:N].';
            GUIpar.Res.FlagAll(GUIpar.Res.IndRM) = 0;

            set(GUIpar.H2,'XData',GUIpar.Res.OrigX(GUIpar.Res.FlagAll==1),...
                          'YData',GUIpar.Res.OrigY(GUIpar.Res.FlagAll==1));

            % red cros
            set(GUIpar.H1,'XData',GUIpar.Res.OrigX(GUIpar.Res.FlagAll==0),...
                          'YData',GUIpar.Res.OrigY(GUIpar.Res.FlagAll==0));

            switch lower(GUIpar.Res.FunBehav)
             case 'i'
                GUIpar = call_pi_fun(GUIpar);
            end

         case 'y'

            fprintf('return the nearest deleted point to the mouse click\n')
            if (isempty(find(GUIpar.Res.FlagAll==0))),
                fprintf(' plot does not contain red croses - skip\n')
            else
                [Xm,Ym] = ginput(1);
                Dist = plane_dist(Xm,Ym,...
                                  GUIpar.Res.OrigX(GUIpar.Res.FlagAll==0),...
                                  GUIpar.Res.OrigY(GUIpar.Res.FlagAll==0));
                [Min,MinInd] = min(Dist);
                F0 = find(GUIpar.Res.FlagAll==0);
                MinInd = F0(MinInd);

                GUIpar.Res.IndRM = setdiff(GUIpar.Res.IndRM, MinInd);
                N = length(GUIpar.Res.OrigX);
                GUIpar.Res.FlagAll = ones(N,1);
                GUIpar.Res.FlagAll(GUIpar.Res.IndRM) = 0;

                set(GUIpar.H2,'XData',GUIpar.Res.OrigX(GUIpar.Res.FlagAll==1),...
                              'YData',GUIpar.Res.OrigY(GUIpar.Res.FlagAll==1));
               % delete the red cross from plot
               set(GUIpar.H1,'XData',GUIpar.Res.OrigX(GUIpar.Res.FlagAll==0),...
                             'YData',GUIpar.Res.OrigY(GUIpar.Res.FlagAll==0));
              
            end

            switch lower(GUIpar.Res.FunBehav)
             case 'i'
                GUIpar = call_pi_fun(GUIpar);
            end

         case 's'
            fprintf('return all deleted points within rectangular area\n')
            if (isempty(find(GUIpar.Res.FlagAll==0))),
                fprintf(' plot does not contain deleted points (red croses) - skip\n')
            else
                Rect = getrect;
                Xmin = Rect(1); Xmax = Rect(1)+Rect(3);
                Ymin = Rect(2); Ymax = Rect(2)+Rect(4);
                MinInd   = find(GUIpar.Res.OrigX>Xmin & GUIpar.Res.OrigX<Xmax & ...
                                GUIpar.Res.OrigY>Ymin & GUIpar.Res.OrigY<Ymax & GUIpar.Res.FlagAll==0);

                GUIpar.Res.IndRM = setdiff(GUIpar.Res.IndRM, MinInd);
                N = length(GUIpar.Res.OrigX);
                GUIpar.Res.FlagAll = ones(N,1);
                GUIpar.Res.FlagAll(GUIpar.Res.IndRM) = 0;

                set(GUIpar.H2,'XData',GUIpar.Res.OrigX(GUIpar.Res.FlagAll==1),...
                              'YData',GUIpar.Res.OrigY(GUIpar.Res.FlagAll==1));
                % delete the red cross from plot
                set(GUIpar.H1,'XData',GUIpar.Res.OrigX(GUIpar.Res.FlagAll==0),...
                              'YData',GUIpar.Res.OrigY(GUIpar.Res.FlagAll==0));
            end

            switch lower(GUIpar.Res.FunBehav)
             case 'i'
                GUIpar = call_pi_fun(GUIpar);
            end

         case 'c'
            fprintf('return all deleted points\n')
            if (isempty(find(GUIpar.Res.FlagAll==0))),
                fprintf(' plot does not contain red croses (deleted points) - skip\n')
            else
                N = length(GUIpar.Res.OrigX);
                GUIpar.Res.FlagAll      = ones(N,1);   % 1 in; 0 crossed
                GUIpar.Res.IndRM    = [];
              
                set(GUIpar.H2,'XData',GUIpar.Res.OrigX(GUIpar.Res.FlagAll==1),...
                              'YData',GUIpar.Res.OrigY(GUIpar.Res.FlagAll==1));
                % delete the red cross from plot
                set(GUIpar.H1,'XData',GUIpar.Res.OrigX(GUIpar.Res.FlagAll==0),...
                              'YData',GUIpar.Res.OrigY(GUIpar.Res.FlagAll==0));
       
            end

            switch lower(GUIpar.Res.FunBehav)
             case 'i'
                GUIpar = call_pi_fun(GUIpar);
            end

         case 'a'
            fprintf(' Add artifial data point\n')
            [ArtX, ArtY] = ginput(1);
            N = length(GUIpar.Res.X);

            GUIpar.Res.X = [GUIpar.Res.X; ArtX];
            GUIpar.Res.Y = [GUIpar.Res.Y; ArtY];
            GUIpar.Res.OrigX = [GUIpar.Res.OrigX; ArtX];
            GUIpar.Res.OrigY = [GUIpar.Res.OrigY; ArtY];
            N = N + 1;
            GUIpar.Res.FlagAll(N) = 1;
            GUIpar.Res.Ind    = [1:1:N].';

            set(GUIpar.H2,'XData',GUIpar.Res.OrigX(GUIpar.Res.FlagAll==1),...
                          'YData',GUIpar.Res.OrigY(GUIpar.Res.FlagAll==1));

            switch lower(GUIpar.Res.FunBehav)
             case 'i'

                GUIpar = call_pi_fun(GUIpar);
            end
 
         case 'f'

            GUIpar=call_pi_fun(GUIpar);
            switch lower(GUIpar.Res.DispFit)
             case 'y'
                % display all residuals
                fprintf('  Table of residuals for best fit\n');
                fprintf('        X          Y      Residuals Y\n');
                fprintf('  %10.3f %10.3f   %8.3f\n',...
                        [GUIpar.Res.X, GUIpar.Res.Y, GUIpar.Res.FunOut{1}.Resid].');
                fprintf('Best fit RMS: %10.3f\n',GUIpar.Res.FunOut{1}.RMS);

             otherwise
                % do nothing
            end


         case 'g'
            % change the function free parameter
            % prompt the user for the value of the free parameter
            if (isempty(GUIpar.Res.FunParInd)),
               fprintf('g option is ignored\n');
            else
               R     = input(sprintf('Type a numeric value for the function free parameter : \n'),'s');
               Value = str2num(R);
               if (isempty(Value)),
                  fprintf('Entered value is illegal - use default value\n');
               else
                  GUIpar.Res.FunPar{GUIpar.Res.FunParInd} = Value;
               end
            end

            switch lower(GUIpar.Res.FunBehav)
             case 'i'
                GUIpar = call_pi_fun(GUIpar);
            end

         case 'm'
            fprintf('\n');
            fprintf('Display Menu:\n');
            plot_rm_menu;
 
         otherwise
            fprintf('\n'); 
            fprintf('Unknown option - options are:\n');
            plot_rm_menu;

    end
    %GUIpar.Res.FlagAll
    GUIpar.Res.Ind   = find(GUIpar.Res.FlagAll==1);
    GUIpar.Res.IndRM = find(GUIpar.Res.FlagAll==0);
    GUIpar.Res.X     = GUIpar.Res.OrigX(GUIpar.Res.Ind);
    GUIpar.Res.Y     = GUIpar.Res.OrigY(GUIpar.Res.Ind);
    set(gcf,'UserData',GUIpar);
end

GUIpar = get(gcf,'UserData');
Res = GUIpar.Res;

return



%---------------------------------------------------------
function GUIpar=call_pi_fun(GUIpar)
%---------------------------------------------------------


% deleted GraphicH (from last call of function)
if (isempty(GUIpar.GraphicH)==1),
   % do nothing
else
   set(GUIpar.GraphicH,'XData',[]);
   set(GUIpar.GraphicH,'YData',[]);
   GUIpar.GraphicH = [];
end

% call function
if (isempty(GUIpar.Res.CallFun)),
   fprintf('No function to call\n')
else
   if (length(find(GUIpar.Res.FlagAll==1))==0),
      fprintf('No data point left - cannot apply function\n');
   else
      varargout{:} = feval(GUIpar.Res.CallFun,...
                           GUIpar.Res.OrigX(GUIpar.Res.FlagAll==1),...
                           GUIpar.Res.OrigY(GUIpar.Res.FlagAll==1),...
                           GUIpar.Res.FunPar{:});

      % check if a graphic handle was added to plot
      AllH = get(gca,'Children');
      GUIpar.GraphicH = setdiff(AllH,[GUIpar.H1;GUIpar.H2]);

      GUIpar.Res.FunOut     = varargout;
   end
end



%--------------------------------------------------------
function plot_rm_menu;
%--------------------------------------------------------

disp('=== Menu ===')
disp('    q - quit')
disp('    x - put a red cross on the nearest point')
disp('    r - put a red cross on all the points inside a rectangular box')
disp('    d - put a red cross on all the points');
disp('    y - return the nearest point with a red cross')
disp('    s - return all the points with a red cross within a rectabgular region')
disp('    c - return all the red croses')
disp('    a - Add artificial point')
disp('    f - call the function')
disp('    g - change the function free parameter')
disp('    m - display this menu')
disp('----------')
fprintf('\n');
