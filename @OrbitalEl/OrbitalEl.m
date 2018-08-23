% A class for orbital elements.
% Package: @OrbitalEl
% Description: A class for orbital elements.
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Mar 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef OrbitalEl
    properties (SetAccess = public)
        T      % Time of periastron
        q      % Periastron distance
        e      % eccntricity
        Om     % Longitude of ascending node
        w      % Argument of perihelion
        i      % inclination
        M
        Epoch  % Epoch
        H
        G
        Name
        Number
        ObjType
        k    = 0.017202098950000;  % Gauss grav. constant (see gauss_grav_const.m)
        UserData
    end
    
    % consider add hidden properties:
    % isaligned
    
%     properties (Hidden = true)
%         %mean1
%         header
%     end
    methods

        %-------------------------
        %--- Class constructor ---
        %-------------------------
        
        function obj=OrbitalEl(varargin)
            % OrbitalEl class constructor
            obj(1).UserData = [];
            %obj = struct_def({'Header','UserData'},varargin{:});
            
        end
        
%         function obj=disp(Head)
%             obj = Head.Header
%         end
        
     
        %--------------------------
        %--- Structre functions ---
        %--------------------------
        function obj=isfield(OrbEl,Field)
            % isfield for OrbitalEl class
            obj = any(strcmp(fieldnames(OrbEl),Field));
        end

        function obj=isstruct(OrbEl)
            % isstruct for OrbitalEl class
            obj = true;  %isstruct(Sim) || isa(Sim,'SIM');
        end

    end % methods
    
    % Static methods
    methods (Static)
        function Ans=isOrbitalEl(Obj)
            % Return true if object is OrbitalEl
            % Description: Check if object is of OrbitalEl class.
            % Input  : - Object
            % Output : - {true|false}.
            %     By : Eran O. Ofek                    Oct 2015
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: S=1; OrbitalEl.isOrbitalEl(S);
            % Reliable: 2
            Ans = isa(Obj,'OrbitalEl');
        end
        
        function OrbEl=struct2orbitalel(St)
            % A structure array into an OrbitalEl object
            % Package: @OrbitalEl
            % Description: Convert a structure array into an OrbitalEl class object.
            % Input  : - A structure array with the appropriate fields.
            % Output : - An OrbitalEl object.
            % See also: OrbitalEl.m
            % License: GNU general public license version 3
            % Tested : Matlab R2015b
            %     By : Eran O. Ofek                    Mar 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: OrbEl=OrbotalEl.struct2orbitalel(Data);
            % Reliable: 2
            %--------------------------------------------------------------------------

            OrbEl = OrbitalEl;
            Nst   = numel(St);
            FN    = fieldnames(St);
            Nfn   = numel(FN);

            for Ist=1:1:Nst
                OrbEl(Ist) = OrbitalEl;
                for Ifn=1:1:Nfn
                    OrbEl(Ist).(FN{Ifn}) = St(Ist).(FN{Ifn});
                end
            end
        end
        
        function OrbEl=elements2orbitalel(varargin)
            % convert elements to orbital elements object
            % Input  : - Either [T, q, e, Om, w, i [k]], where k is optional,
            %            or an OrbitalEl object.
            %            Note that k is the Gaussian gravitational constant.
            %            Default is the constant for the solar system.
            % Output : - An OrbitalEl object.
            % License: GNU general public license version 3
            %     By : Eran O. Ofek                    Mar 2016
            %    URL : http://weizmann.ac.il/home/eofek/matlab/
            % Example: OrbEl = OrbitalEl.elements2orbitalel([T,q,e,Om,w,i]);
            % Reliable: 2

            if (OrbitalEl.isOrbitalEl(varargin{1}))
                OrbEl = varargin{1};
            elseif (isnumeric(varargin{1}))
                % assume: [T, q, e, Om, w, i]
                OrbEl = OrbitalEl;
                OrbEl.T  = varargin{1}(:,1);
                OrbEl.q  = varargin{1}(:,2);
                OrbEl.e  = varargin{1}(:,3);
                OrbEl.Om = varargin{1}(:,4);
                OrbEl.w  = varargin{1}(:,5);
                OrbEl.i  = varargin{1}(:,6);
                if (size(varargin{1},2)>6)
                    OrbEl.k = varargin{1}(:,7);
                end
            else
                error('Unsupported option');
            end
        end

    end
    
    methods
        %-------------------------------
        %--- Selection and searching ---
        %-------------------------------
        
        function SubOrbEl=search(OrbEl,Name)
            % Search OrbitalEl object by name
            
            if (ischar(Name))
                Name = {Name};
            end
            
            if (isnumeric(Name))
                Ind = Util.array.findmany(OrbEl.Number,Name);
            elseif (iscell(Name))
                Nn = numel(Name);
                Ind = [];
                for In=1:1:Nn
                    Ind = [Ind;find(strcmpi(Name{In},OrbEl.Name))];
                end
            else
                error('Unknown Name input');
            end
            
            SubOrbEl = OrbitalEl;
            FN = fieldnames(OrbEl);
            Nfn = numel(FN);
            for Ifn=1:1:Nfn
                Fsize = numel(OrbEl.(FN{Ifn}));
                if (Fsize==0)
                    SubOrbEl.(FN{Ifn}) = [];
                elseif (Fsize==1)
                    SubOrbEl.(FN{Ifn}) = OrbEl.(FN{Ifn});
                else
                    SubOrbEl.(FN{Ifn}) = OrbEl.(FN{Ifn})(Ind);
                end
            end
            
        end
        
        
        
        function A=a(OrbEl)
            % Calculate the semi-major axis [AU]
            
            A = OrbEl.q./(1-OrbEl.e);
        end

        function App=Q(OrbEl)
            % Calculate the apiastron distance [AU]
            
            App = (1+OrbEl.e).*OrbEl.q./(1-OrbEl.e);
        end
     
        function B=b(OrbEl)
            % Calculate the semi-minor axis
            
            B = OrbEl.q.*sqrt(1-OrbEl.e.^2)./(1-OrbEl.e);
            B(OrbEl.e>=1) = Inf;
        end
        
        function P=p(OrbEl)
            % Calculate the semi-latus rectum axis
            
            P = OrbEl.q.*(1+OrbEl.e);
        end
        
        function r=eccanom2radius(OrbEl,E)
            % Eccentric anomaly to radius vector
            
            r = (1 - OrbEl.e.*cos(E)).*OrbEl.q./(1-OrbEl.e);
        end
        
        function r=trueanom2radius(OrbEl,Nu)
            % Eccentric anomaly to radius vector
            
            Nu    = Nu.*ones(size(OrbEl.q));
            FlagP = OrbEl.e==1;
            FlagH = OrbEl.e>1;
            r        = OrbEl.q.*(1+OrbEl.e)./(1 + OrbEl.e.*cos(Nu));
            r(FlagP) = OrbEl.q(FlagP).*(1 + tan(0.5.*Nu(FlagP)).^2);
            r(FlagH) = OrbEl.q(FlagH).*(1 + OrbEl.e(FlagH))./(1 + OrbEl.e(FlagH).*cos(Nu(FlagH)));
            
        end
        
        function Period=P(OrbEl)
            % Calculate the orbital period [day]
            Period = (2.*pi./OrbEl.k).*(OrbEl.q./(1-OrbEl.e)).^1.5;
            % Set period to infinity
            Period(OrbEl.e>=1) = Inf;
        end
        
        function N=n(OrbEl)
            % Calculate the mean motion (n) [rad/day]
            N = 2.*pi./P(OrbEl);
            N(OrbEl.e>=1) = NaN;
        end
        
        function MeanAnom = meananom(OrbEl,t)
            % Calculate the mean anomaly
           
            MeanAnom = n(OrbEl).*(t-OrbEl.T);
        end
        
        function Nu=eccanom2trueanom(OrbEl,E)
            % Eccentric anomaly to True anomaly
            
            Nu = 2.*atan(sqrt( (1+OrbEl.e)./(1-OrbEl.e) ).*tan(0.5.*E));
            % is this correct for e>1
        end
        
        function E=trueanom2eccanom(OrbEl,Nu)
            % True anomaly to Eccentric anomaly
            
            E = 2.*atan(sqrt( (1-OrbEl.e)./(1+OrbEl.e) ).*tan(0.5.*Nu));
            % is this correct for e>1
        end
        
        function dNUdt=nudot(OrbEl,Nu)
            % Calculate the time derivative of the true anomaly
            % Description: Calculate the time derivative of the true anomaly.
            %              Correct only for e<1
            
            E = trueanom2eccanom(OrbEl,Nu); % Eccentric anomaly
            dNUdt = n(OrbEl).*sqrt(1-OrbEl.e.^2)./((1 - OrbEl.e.*cos(E)).^2);
        end
        
        function drdt=rdot(OrbEl,Nu)
            % Calculate the time derivative of the radius vector
            % Description: Calculate the time derivative of the radius
            %              vector.
            %              correct only for e<1
            % Input  : - OrbitalEl object.
            %          - True anomaly [rad].
            % Output : - dr/dt
            
            E = trueanom2eccanom(OrbEl,Nu); % Eccentric anomaly
            drdt = n(OrbEl).*a(OrbEl).*OrbEl.e.*sin(E)./(1 - OrbEl.e.*cos(E));
            
        end
        
        function V=r2vel(OrbEl,r)
            % Calculate orbital velocity from radius vector
            % Package: @OrbitalEl
            % Description: Calculate orbital velocity from radius vector
            %              Correct only for e<1
            % Input  : - OrbitalEl object.
            %          - Radius vector
            % Output : - Velocity
            
            
            A = a(OrbEl);
            V = (sqrt(2).*2.*pi.*A./P(OrbEl)).*sqrt( 1./r - 1./(2.*A) );
        end
        
        function V=trueanom2vel(OrbEl,Nu)
            % Calculate orbital velocity from true anomaly
            % Description: Calculate orbital velocity from true anomaly
            %              Correct only forVelocity from eccentric anomaly e<1
            % Input  : - OrbitalEl object.
            %          - True anomaly [rad]
            % Output : - Velocity
            
            r = trueanom2radius(OrbEl,Nu);
            A = a(OrbEl);
            V = (sqrt(2).*2.*pi.*A./P(OrbEl)).*sqrt( 1./r - 1./(2.*A) );
        end
        
        function V=eccanom2vel(OrbEl,E)
            % Velocity from Eccentric anomaly
            % Description: Velocity from Eccentric anomaly
            %              Correct only for e<1
            % Input  : - OrbitalEl object
            %          - Eccentric anomaly [rad].
            % Output : - Velocity
            
            r = eccanom2radius(OrbEl,E);
            A = a(OrbEl);
            V = (sqrt(2).*2.*pi.*A./P(OrbEl)).*sqrt( 1./r - 1./(2.*A) );
        end
        
        
        % Kepler Equation
        function [E,Nu,R]=kepler_elliptic(OrbEl,t,Tol)
            % Description: Solve Kepler equation
            
            if (nargin<3)
                Tol = 1e-8;
            end
            MeanAnom = meananom(OrbEl,t);
            E = kepler_elliptic_fast(MeanAnom,OrbEl.e,Tol);
            if (nargout>1)
                Nu = eccanom2trueanom(OrbEl,E);
                if (nargout>2)
                    R = eccanom2radius(OrbEl,E);
                end
            end
        end
        
        % kepler_parabolic
        % kepler_hyperbolic
        % elements2position
        % position2elements
        
        
        function BodyPos=trueanom2pos(OrbEl,Nu,R)
            % True anomaly to rectangular position
            % Package: @OrbitalEl
            % Description: True anomaly to rectangular position
            % Input  : - OrbitalEl object.
            %          - True anomaly [rad].
            %          - Optional radius vector. If not given will be
            %            calculated from the True anaomaly.
            % Output : * Either a single 3 column matrix of [X,Y,Z] or 
            %            X,Y,Z. Units the same as the radius vector units.
            
            if (nargin<3)
                % calc radius vector
                R = trueanom2radius(OrbEl,Nu);
            end
            BodyPos = trueanom2pos(R,Nu,OrbEl.Om,OrbEl.w,OrbEl.i);
        end
        
        function TI=thiele_innes(OrbEl)
            % Convert orbital elements to Thiele-Innes elements
            % Package: @OrbitalEl
            % Description: Convert orbital elements to Thiele-Innes
            %              orbital elements.
            % Input  : - OrbitalEl object
            % Output : - Structure array with Thiele-Innes elements.
            % See also: thiele_innes2el.m
            
            TI=thiele_innes(a(OrbEl),OrbEl.w,OrbEl.Om,OrbEl.i);
        
        end
    end
end

            
