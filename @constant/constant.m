%--------------------------------------------------------------------------
% constant class                                                     class
% Description: A static class for physical and astrophysical constants.
%              This class include many static function for constants.
%              Type "constant." followed by <tab> to see the full list of
%              functions.
%              constant.all return a structure array of selected constants.
% Input  : null
% Output : null
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Aug 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%--------------------------------------------------------------------------

classdef constant 
    
    
    % get all physical constants
    methods (Static)
        
        function Const=all(System)
            % Return the value of selected constants in a structure array
            % Package: @constant
            % Description: Return the value of selected constants
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - Structure array with all the selected constants
            %            and their value, units, error, and formula.
            % Package: @constant
            % Input  : - System. Default is 'cgs'.
            % Output : - A structure array with element per constant.
            
            List = {'a','a0','alpha','amu','au','c','day','e',...
                    'EarthM','EarthR','eps0','G','h','hbar','kB','ly',...
                    'me','mp','mu0','NA','pc','R','RAD','re',...
                    'Rydberg','sigma','sigmaT','SunL','SunM','SunR'};
                
            Nl = numel(List);
            Const = Util.struct.struct_def({'Name','Const','Units','Error','Form'},Nl,1);
            for Il=1:1:Nl
                Const(Il).Name = List{Il};
                [Const(Il).Const,Const(Il).Units,Const(Il).Error,Const(Il).Form] = constant.(List{Il});
            end                
                    
            
        end
        
        function [St,StUnits]=all_st(varargin)
            % Return the values of selected constants in a structure
            % Package: @constant
            % Description: Return the value of selected constants
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - Structure array with all the selected constants
            %            and their value, units, error, and formula.
            % Package: @constant
            % Input  : - System. Default is 'cgs'.
            % Output : - A structure with constants in each field.
            %          - A structure with the constant units in each field.
            % Example: C=constant.all_st
            
            Con  = constant.all(varargin{:});
            Ncon = numel(Con);
            for Icon=1:1:Ncon
                St.(Con(Icon).Name) = Con(Icon).Const;
                StUnits.(Con(Icon).Name) = Con(Icon).Units;
            end
            
        end
        
        function all_var(varargin)
            % Assign selected constants into the workspace
            % Package: @constant
            % Description: Assign into the workspace variables with the
            %              constant names and values.
            % Input  : - System. Default is 'cgs'.
            % Output : * null
            % Example: constant.all_var
            
            WorkSpace = 'base';
            
            Con  = constant.all(varargin{:});
            Ncon = numel(Con);
            for Icon=1:1:Ncon
                assignin(WorkSpace, Con(Icon).Name, Con(Icon).Const);
            end
            
        end
        
        function all_var_clear
            % Clear selected constants from the workspace
            % Package: @constant
            % Description: Use this function in order ro clear the
            %              variables assigin into the workspace by the
            %              constant.all_var function.
            % Example: constant.all_var; whos, constant.all_var_clear
            
            WorkSpace = 'base';
            
            Con  = constant.all;
            Ncon = numel(Con);
            for Icon=1:1:Ncon
                evalin(WorkSpace, sprintf('clear %s',Con(Icon).Name));
            end
            
        end
        
        
        
    end % statics
    
    % specific physical constants
    methods (Static)
        
        function [Const,Units,Error,Form]=au(System)
            % Return the value of the au
            % Package: @constant
            % Description: Return the value of the au
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 1.49597870691e13;      
                Units = 'cm';
            else
                % SI
                Const = 1.49597870691e11;
                Units = 'm';
            end
            Error = NaN;
            Form  = '';
              
        end % constant.au function
        
        function [Const,Units,Error,Form]=G(System)
            % Return the value of the G
            % Package: @constant
            % Description: Return the value of the G
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const   = 6.67259e-8;
                Units = 'cm^3 * gr^-1 * s^-2';
            else
                % SI
                Const = 6.67259e-11;
                Units  = 'm^3 * kg^-1 * s^-2';
            end
            Error = NaN;
            Form  = '';
              
        end % constant.G function
        
        function [Const,Units,Error,Form]=c(System)
            % Return the value of the c
            % Package: @constant
            % Description: Return the value of the c
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const   = 29979245800;
                Units = 'cm * s^-1';
            else
                % SI
                Const = 299792458;
                Units  = 'm * s^-1';
            end
            Error = 0;
            Form  = '';
           
        end % constant.c function
        
        function [Const,Units,Error,Form]=h(System)
            % Return the value of the h
            % Package: @constant
            % Description: Return the value of the h
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 6.6260755e-27;
                Units = 'cm^2 * gr * s^-1';
            else
                % SI
                Const = 6.6260755e-34;
                Units = 'm^2 * kg * s^-1';
            end
            Error = NaN;
            Form  = '';
              
        end % constant.h function
        
        function [Const,Units,Error,Form]=hbar(System)
            % Return the value of the hbar
            % Package: @constant
            % Description: Return the value of the hbar
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 1.054571800e-27;
                Units = 'cm^2 * gr * s^-1';
            else
                % SI
                Const = 1.054571800e-34;
                Units = 'm^2 * kg * s^-1';
            end
            Error = 1.3e-8;
            Form  = 'h/(2*pi)';
              
        end % constant.h function
        
        function [Const,Units,Error,Form]=e(System)
            % Return the value of the e
            % Package: @constant
            % Description: Return the value of the e
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 4.8032068e-10;
                Units = 'esu';
            else
                % SI
                Const = 1.60217733e-19;
                Units  = 'C';
            end
            Error = NaN;
            Form  = '';
              
        end % constant.e function
        
        function [Const,Units,Error,Form]=alpha(System)
            % Return the value of the alpha
            % Package: @constant
            % Description: Return the value of the alpha
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 7.2973525698e-3;
                Units = '';
            else
                % SI
                Const = 7.2973525698e-3;
                Units = '';
            end
            Error = NaN;
            Form  = 'e^2/(2*eps0*h*c)';
              
        end % constant.alpha function
        
        function [Const,Units,Error,Form]=Rydberg(System)
            % Return the value of the Rydberg constant
            % Package: @constant
            % Description: Return the value of the Rydberg constant
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 10.973731568539e4;
                Units = '';
            else
                % SI
                Const = 10.973731568539e6;
                Units = 'm^-1';
            end
            Error = NaN;
            Form  = 'me * e^4/(8*eps0^2 * h^3 * c)';
              
        end % constant.Rydberg function
        
        function [Const,Units,Error,Form]=eps0(System)
            % Return the value of the eps0 constant
            % Package: @constant
            % Description: Return the value of the eps0 constant
            %              Vacuum permittivity/permittivity of free space
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = NaN;
                Units = '';
            else
                % SI
                Const = 8.854187817620e-12;
                Units = 'F * m^-1';
            end
            Error = NaN;
            Form  = '1/(mu0 * c^2) = e^2/(2*alpha*h*c)';
              
        end % constant.eps0 function
        
        function [Const,Units,Error,Form]=mu0(System)
            % Vacuum permeability/magnetic constant
            % Package: @constant
            % Description: Return the value of the mu0 
            %              Vacuum permeability/magnetic constant.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = NaN;
                Units = '';
            else
                % SI
                Const = 4.*pi.*1e-7;
                Units  = 'H * m^-1';
            end
            Error = NaN;
            Form  = 'H * m^-1';
              
        end % constant.mu0 function
        
        function [Const,Units,Error,Form]=mp(System)
            % The proton mass
            % Package: @constant
            % Description: Return the value of the mp
            %              Proton rest mass
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 1.6726231e-24;
                Units = 'g';
            else
                % SI
                
                Const = 1.6726231e-27;
                Units  = 'kg';
            end
            Error = NaN;
            Form  = '';
              
        end % constant.mp function
        
        function [Const,Units,Error,Form]=me(System)
            % The mass of the electron
            % Package: @constant
            % Description: Return the value of the me
            %              Electron rest mass
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 9.10938997e-28;
                Units = 'g';
            else
                % SI
                Const = 9.10938997e-31;
                Units  = 'kg';
            end
            Error = 4.9e-10;
            Form  = '';
              
        end % constant.me function
        
        function [Const,Units,Error,Form]=amu(System)
            % The atomic mass unit constant
            % Package: @constant
            % Description: Return the value of the amu
            %              Atomic mass unit (m(12C)/12)
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 1.660538921e-24;
                Units = 'g';
            else
                % SI
                Const = 1.660538921e-27;
                Units = 'kg';
            end
            Error = 4.4e-8;
            Form  = '';
              
        end % constant.amu function
        
        function [Const,Units,Error,Form]=NA(System)
            % Avogadro constant
            % Package: @constant
            % Description: Return the value of the NA
            %              Avogadro constant/Avogadro number
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 6.0221412927e23;
                Units = 'mol^-1';
            else
                % SI
                Const = 6.0221412927e23;
                Units = 'mol^-1';
            end
            Error = 4.5e-10;
            Form  = '';
              
        end % constant.NA function
        
        function [Const,Units,Error,Form]=kB(System)
            % The Boltzman constant
            % Package: @constant
            % Description: Return the value of the KB
            %              Boltzmann constant.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 1.380648813e-16;
                Units = 'cm^2 * gr * s^-2 * K^-1';
            else
                % SI
                Const = 1.380648813e-23;
                Units = 'm^2 * kg * s^-2 * K^-1';
            end
            Error = 9.4e-9;
            Form  = 'R/NA';
              
        end % constant.kB function
        
        function [Const,Units,Error,Form]=R(System)
            % The Molar gas constant
            % Package: @constant
            % Description: Return the value of the R
            %              Molar gas constant.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 8.314598e7;
                Units = 'erg *mol^-1 *K^-1';
            else
                % SI
                Const = 8.314598;
                Units = 'J *mol^-1 *K^-1';
            end
            Error = NaN;
            Form  = '';
              
        end % constant.R function
        
        function [Const,Units,Error,Form]=sigma(System)
            % The Stefan-Boltzmann constant
            % Package: @constant
            % Description: Return the value of the sigma
            %              Stefan-Boltzmann constant.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 5.67037321e-5;
                Units = 'erg * cm^-2 * K^-4';
            else
                % SI
                Const = 5.67037321e-8;
                Units  = 'W * m^-2 * K^-4';
            end
            Error = 3.7e-8;
            Form  = '2 * pi^5 * kB^4/(15*h^3*c^2)';
              
        end % constant.sigma function
        
        function [Const,Units,Error,Form]=a(System)
            % The radiation constant
            % Package: @constant
            % Description: Return the value of the a
            %              Radiation constant.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 7.5657316369e-15;
                Units = 'erg * cm^-3 * K^-4';
            else
                % SI
                Const = 7.5657316369e-16;
                Units  = 'J * m^-3 * K^-4';
            end
            Error = 3.7e-8;
            Form  = '8 * pi^5 * kB^4/(15*h^3*c^3)';
              
        end % constant.a function
        
        function [Const,Units,Error,Form]=sigmaT(System)
            % Thompson cross-section constant
            % Package: @constant
            % Description: Return the value of the SigmaT
            %              Thompson cross-section constant.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 6.65245854533e-25;
                Units = 'cm^2';
            else
                % SI
                Const = 6.65245854533e-29;
                Units = 'm^2';
            end
            Error = NaN;
            Form  = '8 * pi * alpha^2 * h^2 * c^2/(6*pi*me)';
              
        end % constant.sigmaT function
        
        function [Const,Units,Error,Form]=re(System)
            % Classical radius of the electron
            % Package: @constant
            % Description: Return the value of the re
            %              Classical electron radius.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 2.817940289458e-13;
                Units = 'cm';
            else
                % SI
                Const = 2.817940289458e-15;
                Units = 'm';
            end
            Error = 2.1e-11;
            Form  = '';
              
        end % constant.re function
        
        function [Const,Units,Error,Form]=a0(System)
            % Bohr radius
            % Package: @constant
            % Description: Return the value of the a0
            %              Bohr radius.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 5.291772109217e-9;
                Units = 'cm';
            else
                % SI
                Const = 5.291772109217e-11;
                Units = 'm';
            end
            Error = 3.2e-12;
            Form  = 'h/(2*pi*me*c*alpha)';
              
        end % constant.a0 function
        
        function [Const,Units,Error,Form]=SunM(System)
            % The Sun mass
            % Package: @constant
            % Description: Return the value of the SunM
            %              Sun mass.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 1.98892e33;
                Units = 'g';
            else
                % SI
                Const = 1.98892e30;
                Units = 'kg';
            end
            Error = 1.3e-4;
            Form  = '';
              
        end % constant.SunM function
        
        function [Const,Units,Error,Form]=SunR(System)
            % Return the value of the Sun Radius
            % Package: @constant
            % Description: Return the value of the SunR
            %              Sun radius.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 6.96342e10;
                Units = 'cm';
            else
                % SI
                Const = 6.96342e8;
                Units = 'm';
            end
            Error = 9.3e-5;
            Form  = '';
              
        end % constant.SunR function
        
        function [Const,Units,Error,Form]=EarthM(System)
            % Return the value of the Earth Mass
            % Package: @constant
            % Description: Return the value of the EarthM
            %              Earth mass.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 5.9722e27;
                Units = 'g';
            else
                % SI
                Const = 5.9722e24;
                Units = 'kg';
            end
            Error = NaN;
            Form  = '';
              
        end % constant.EarthM function
        
        function [Const,Units,Error,Form]=EarthR(System)
            % Return the value of the Earth mean radius.
            % Package: @constant
            % Description: Return the value of the EarthR
            %              Earth mean radius.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 6371.00e5;
                Units = 'cm';
            else
                % SI
                Const = 6371.00e3;
                Units = 'm';
            end
            Error = NaN;
            Form  = '';
              
        end % constant.EarthR function
        
        function [Const,Units,Error,Form]=JupiterM(System)
            % Return the value of the Jupiter Jupiter mass
            % Package: @constant
            % Description: Return the value of the JupiterM
            %              Jupiter mass.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 1.8991741e30;
                Units = 'g';
            else
                % SI
                Const = 1.8991741e27;
                Units = 'kg';
            end
            Error = NaN;
            Form  = '';
              
        end % constant.JupiterM function
        
        function [Const,Units,Error,Form]=JupiterR(System)
            % Return the value of the Jupiter mean radius.
            % Package: @constant
            % Description: Return the value of the JupiterR
            %              Jupiter mean radius.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 69911e5;
                Units = 'cm';
            else
                % SI
                Const = 69911e3;
                Units = 'm';
            end
            Error = NaN;
            Form  = '';
              
        end % constant.JupiterR function
        
        function [Const,Units,Error,Form]=SunL(System)
            % Return the value of the Sun luminosity.
            % Package: @constant
            % Description: Return the value of the SunL
            %              Sun luminosity.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 3.839e33;
                Units = 'erg * s^-1';
            else
                % SI
                Const = 3.839e26;
                Units = 'W';
            end
            Error = NaN;
            Form  = '';
              
        end % constant.SunL function
        
        function [Const,Units,Error,Form]=pc(System)
            % Return the value of the Parsec
            % Package: @constant
            % Description: Return the value of the pc
            %              Parsec.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 3.08567758e18;
                Units = 'cm';
            else
                % SI
                Const = 3.08567758e16;
                Units = 'm';
            end
            Error = NaN;
            Form  = '';
              
        end % constant.pc function
        
        function [Const,Units,Error,Form]=ly(System)
            % Return the value of the Light year
            % Package: @constant
            % Description: Return the value of the ly
            %              Light year.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 9.4607304725808e17;
                Units = 'cm';
            else
                % SI
                Const = 9.4607304725808e15;
                Units = 'm';
            end
            Error = NaN;
            Form  = '';
              
        end % constant.ly function
        
        function [Const,Units,Error,Form]=k(System)
            % Return the value of the Gauss gravitational constant
            % Package: @constant
            % Description: Return the value of the k
            %              Gauss gravitational constant.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 0.01720209895;
                Units = '';
            else
                % SI
                Const = 0.01720209895;
                Units = '';
            end
            Error = NaN;
            Form  = '';
              
        end % constant.k function
        
        function [Const,Units,Error,Form]=PlanckMass(System)
            % Return the value of the Planck Mass
            % Package: @constant
            % Description: Return the value of the Planck Mass
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 2.176470e-5;
                Units = 'g';
            else
                % SI
                Const = 2.176470e-8;
                Units = 'kg';
            end
            Error = 2.3e-5;
            Form  = '(hbar*c/G)^0.5';
              
        end % constant.PlanckMass function
        
        function [Const,Units,Error,Form]=PlanckT(System)
            % Return the value of the Planck Temperature
            % Package: @constant
            % Description: Return the value of the Planck Temperature
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 1.416808e32;
                Units = 'K';
            else
                % SI
                Const = 1.416808e32;
                Units = 'K';
            end
            Error = 2.3e-5;
            Form  = '((hbar*c^5/G)^0.5)/kB';
              
        end % constant.PlanckT function
        
        function [Const,Units,Error,Form]=PlanckLength(System)
            % Return the value of the Planck Length
            % Package: @constant
            % Description: Return the value of the Planck Length
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 1.616229e-33;
                Units = 'cm';
            else
                % SI
                Const = 1.616229e-35;
                Units = 'm';
            end
            Error = 2.3e-5;
            Form  = '(hbar/(mp*c))';
              
        end % constant.PlanckLength function
        
        function [Const,Units,Error,Form]=PlanckTime(System)
            % Return the value of the PlanckTime
            % Package: @constant
            % Description: Return the value of the Planck Time
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative Error
            %          - Formula
            
            if (nargin==0)
                System = true;
            else
                if strcmp(System,'SI')
                    System = false;
                else
                    System = true;
                end
            end
            if (System)
                % cgs
                Const = 5.39116e-44;
                Units = 's';
            else
                % SI
                Const = 5.39116e-44;
                Units = 's';
            end
            Error = 2.3e-5;
            Form  = '(hbar*G/(c.^5)).^0.5';
              
        end % constant.PlanckTime function
        

        
    end % Static
    
    
    % specific astronomical constants
    methods (Static)
        function [Const,Units,Error,Form]=day(System)
            % Return the value of the civil day
            % Package: @constant
            % Description: Return the value of the civil day
            %              Civil day.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative error
            %          - Formula
            
            Const = 86400;
            Units = 's';
            if (nargout>2)
                Error = 0;
                Form  = '';
            end
        end
        
        function [Const,Units,Error,Form]=RAD(System)
            % Return the value of deg in radians
            % Package: @constant
            % Description: Return the value of deg in radians.
            %              Civil day.
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative error
            %          - Formula
            
            Const = 180./pi;
            Units = 'deg';
            if (nargout>2)
                Error = 0;
                Form  = '';
            end
        end
        
        function [Const,Units,Error,Form]=SunMercuryMassRatio(System)
            % Return the value of Sun/Mercury mass ratio
            % Package: @constant
            % Description: Return the value of Sun/Mercury mass ratio
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative error
            %          - Formula
            
            Const = 6023600;
            Units = '';
            if (nargout>2)
                Error = 0;
                Form  = '';
            end
        end
        
        function [Const,Units,Error,Form]=SunVenusMassRatio(System)
            % Return the value of Sun/Venus mass ratio
            % Package: @constant
            % Description: Return the value of Sun/Venus mass ratio
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative error
            %          - Formula
            
            Const = 408523.71;
            Units = '';
            if (nargout>2)
                Error = 0;
                Form  = '';
            end
        end
        
        function [Const,Units,Error,Form]=SunMarsMassRatio(System)
            % Return the value of Sun/Mars mass ratio
            % Package: @constant
            % Description: Return the value of Sun/Mars mass ratio
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative error
            %          - Formula
            
            Const = 3098708;
            Units = '';
            if (nargout>2)
                Error = 0;
                Form  = '';
            end
        end
        
        function [Const,Units,Error,Form]=SunJupiterMassRatio(System)
            % Return the value of Sun/Jupiter mass ratio
            % Package: @constant
            % Description: Return the value of Sun/Jupiter mass ratio
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative error
            %          - Formula
            
            Const = 1047.3486;
            Units = '';
            if (nargout>2)
                Error = 0;
                Form  = '';
            end
        end
        
        function [Const,Units,Error,Form]=SunSaturnMassRatio(System)
            % Return the value of Sun/Saturn mass ratio
            % Package: @constant
            % Description: Return the value of Sun/Saturn mass ratio
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative error
            %          - Formula
            
            Const = 3497.898;
            Units = '';
            if (nargout>2)
                Error = 0;
                Form  = '';
            end
        end
        
        function [Const,Units,Error,Form]=SunUranusMassRatio(System)
            % Return the value of Sun/Uranus mass ratio
            % Package: @constant
            % Description: Return the value of Sun/Uranus mass ratio
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative error
            %          - Formula
            
            Const = 22902.98;
            Units = '';
            if (nargout>2)
                Error = 0;
                Form  = '';
            end
        end
        
        function [Const,Units,Error,Form]=SunNeptuneMassRatio(System)
            % Return the value of Sun/Neptune mass ratio
            % Package: @constant
            % Description: Return the value of Sun/Neptune mass ratio
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative error
            %          - Formula
            
            Const = 19412.24;
            Units = '';
            if (nargout>2)
                Error = 0;
                Form  = '';
            end
        end
        
        function [Const,Units,Error,Form]=SunEarthMoonMassRatio(System)
            % Return the value of Sun/(Earth+Moon) mass ratio
            % Package: @constant
            % Description: Return the value of Sun/(Earth+Moon) mass ratio
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative error
            %          - Formula
            
            Const = 328900.561400;
            Units = '';
            if (nargout>2)
                Error = 0;
                Form  = '';
            end
        end
        
        function [Const,Units,Error,Form]=EarthMoonMassRatio(System)
            % Return the value of Earth/Moon mass ratio
            % Package: @constant
            % Description: Return the value of Earth/Moon mass ratio
            % Input  : - System: 'cgs'|'SI'. Default is 'cgs'.
            % Output : - The value of the constant
            %          - Units
            %          - Relative error
            %          - Formula
            
            Const = 81.30056;
            Units = '';
            if (nargout>2)
                Error = 0;
                Form  = '';
            end
        end
        
    end % static
    
end % end class
            
