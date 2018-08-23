function [JD,ExpTime,Sim]=julday(Sim,varargin)
% Use header information to calculate Julian Day of an HEAD/SIM object.
% Package: @HEAD
% Description: Get or calculate the Julian day from an Header object
%              or e.g., a SIM images object. The JD is calculated based on
%              the header information. This function is looking for
%              the time and date in a list of header keywords and convert
%              it to JD. The program can also use the exposure time to
%              calculate the mid exposure time.
% Input  : - An Header object (or e.g., SIM object).
%            See sim_julday.m for a function that can get as input
%            a FITS file name.
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Dictionary' - Dictionary of header keywords to retrieve
%                           and how to intrepret them. This is a two column
%                           cell array. The first column is the header
%                           keyword name, and the second is its type.
%                           Possible types are: 'JD', 'MJD', 'date',
%                           'hour'.
%                           Default is {'JD','JD'; 'OBSJD','JD'; 'MJD','MJD'; 'OBSMJD', 'MJD';'UTC-OBS','date'; 'DATE-OBS','date'; 'DATE','date'; 'UTC','time'; 'OBS-DATE','date'; 'TIME-OBS','time'}
%                           The program will prioretize the keywords by
%                           the order. 'date' type keywords are either
%                           'YYYY-MM-DD HH:MM:SS.frac' like objects
%                           or 'YYYY-MM-DD'. In the later case, the program
%                           will look also for 'hour' object 'HH:MM:SS'
%            'AddDic'     - An additional user dicetionary to be appended
%                           at the begining of the dictionary.
%                           Default is {}.
%            'FunParKey'  - If the date/time is not found within the
%                           dictionary keywords, then the program will
%                           read the keywords specified in this cell
%                           array and use them to calculate the JD
%                           using some formula.
%                           e.g., {'a','b'}. Default is [].
%            'Fun'        - Function to calculate the JD.
%                           Default is @(a,b) (a+b)./86400;
%            'ExpTimeKey' - A cell array of header keywords that may
%                           contain the exposure time.
%                           Default is {'AEXPTIME','EXPTIME'}.
%                           The keywords will prioretize by their order.
%            'DefExpTime' - Default exposure time to use if ExpTimeKey
%                           is not found in the image header. Default is 0.
%                           If empty, then will return an error if exposure
%                           time is not found.
%            'ExpTimeUnits'- Cell array indicating the units of the exposure
%                           time keywords. Default is {'s','s'}.
%                           See convert.units.m for options.
%            'OutTime'    - The output time corresponding to:
%                           {'mid','start','end'} of the exposure.
%                           Default is 'mid'.
%            'OutType'    - Output time type {'JD','MJD'}. Default is 'JD'.
%            'UpdateHead' - Update/write JD to SIM header {true|false}.
%                           Default is true.
%            'UpdateKey'  - Name of header keyword in which to write the
%                           JD. Default is 'JD_MID'.
% Output : - Vector of JDs per input image.
%          - Vector of ExpTime per input image.
%          - An Header object array with the updated JD keyword.
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Oct 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: JD=julday(Sim);
% Reliable: 2
%--------------------------------------------------------------------------
import Util.cell.*
import Util.string.*

%ImageField  = 'Im';
HeaderField = 'Header';
%FileField   = 'ImageFileName';
%MaskField   = 'Mask';
%BackImField = 'BackIm';
%ErrImField  = 'ErrIm';

DefV.Dictionary   = {'JD','JD'; 'OBSJD','JD'; 'MJD','MJD'; 'OBSMJD', 'MJD';'UTC-OBS','date'; 'DATE-OBS','date'; 'DATE','date'; 'UTC','time'; 'OBS-DATE','date'; 'TIME-OBS','time'};
DefV.AddDic       = {}; %{'UTC-OBS','date'};
DefV.FunParKey    = [];
DefV.Fun          = @(a,b) (a+b)./86400;
DefV.ExpTimeKey   = {'AEXPTIME','EXPTIME'};
DefV.DefExpTime   = 0;
DefV.ExpTimeUnits = {'s','s'};
DefV.OutTime      = 'mid';   % {'mid','start','end'}
DefV.OutType      = 'jd';    % {'JD','MJD'}
DefV.UpdateHead   = true;
DefV.UpdateKey    = 'JD_MID';
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


InPar.Dictionary = [InPar.AddDic; InPar.Dictionary];


if (isstruct(Sim) || SIM.issim(Sim) || isastspec(Sim))
    % no need to read images
else
    Sim = images2sim(Sim,varargin{:});
end
Nsim = numel(Sim);

%[CellTime] = sim_getkeyvals(Sim,InPar.Dictionary(:,1).','ConvNum',false);
[CellTime] = mgetkey(Sim,InPar.Dictionary(:,1).',true,false);
%[CellExpTime] = sim_getkeyvals(Sim,InPar.ExpTimeKey,'ConvNum',false);
[CellExpTime] = mgetkey(Sim,InPar.ExpTimeKey,true,false);


JD      = zeros(size(Sim));
ExpTime = zeros(size(Sim)).*NaN;
for Isim=1:1:Nsim
    % read all relevantr header keywords into a cell array
    CellTime1 = CellTime(Isim,:);
    % look for the first not NaN value
    Inn = find(~isnan_cell(CellTime1),1,'first');
    
    
    %Inn = find(cellfun(@isnan,CellTime,'UniformOutput',false)==false,1,'first');
    if (isempty(Inn))
        % no time found in header based on the time dictionary
        % attempt using Fun
        if (~isempty(InPar.FunParKey))
            %ValFK    = sim_getkeyvals(Sim,InPar.FunParKey,'ConvNum',true);
            ValFK    = mgetkey(Sim,InPar.FunParKey,true,true);
            JD(Isim) = InPar.Fun(ValFK{1},ValFK{2});
        else
            % no time is available
            JD(Isim) = NaN;
        end
    else
        % time found in header
        switch lower(InPar.Dictionary{Inn,2})
            case 'jd'
                JD(Isim) = str2double_check(CellTime1{Inn});
            case 'mjd'
                JD(Isim) = str2double_check(CellTime1{Inn}) + 2400000.5;
            case 'date'
                %DateVec = date_str2vec(CellTime1{Inn});
                DateVec = convert.str2date(CellTime1{Inn});
                if (length(DateVec)<6)
                    % need to read hour separtly
                    Itime=find(~isnan_cell(CellTime1(:)) & strcmp(InPar.Dictionary(:,2),'time'),1,'first');
                    
                    if (isempty(Itime))
                        JD(Isim) = 0; %NaN;
                    else
                        JD(Isim) = convert.date2jd(DateVec(:,[3 2 1])) + convert.hour_str2frac(CellTime1{Itime});
                    end
                    
                else
                    JD(Isim) = convert.date2jd(DateVec(:,[3 2 1 4 5 6]));
                end
                
            case 'hour'
                % do nothing - treated in 'date'
            otherwise
                error('Not in time keywords dictionary');
        end
    end
    
    
    % calculate mid exposure time
    switch lower(InPar.OutTime)
        case 'start'
            % do nothing
        otherwise
            CellExpTime1 = CellExpTime(Isim,:);
            % look for the first not NaN value
            Inn = find(~isnan_cell(CellExpTime1),1,'first');
            
            if (isempty(Inn))
                if (isempty(InPar.DefExpTime))
                    error('Exposure time not found');
                else
                    ExpTime(Isim) = InPar.DefExpTime;
                    Inn = 1;
                end
            else
                ExpTime(Isim) = str2double_check(CellExpTime1{Inn});
            end
            
            ConvFactor = convert.units(InPar.ExpTimeUnits{Inn},'day');
            switch lower(InPar.OutTime)
                case 'mid'
                    JD(Isim) = JD(Isim) + 0.5.*ConvFactor.*ExpTime(Isim);
                case 'end'
                    JD(Isim) = JD(Isim) + ConvFactor.*ExpTime(Isim);
                otherwise
                    error('Unknown OutTime option');
            end
            
    end
    
    if (InPar.UpdateHead)
        Sim(Isim) = add_key(Sim(Isim),...
                            sprintf('%s',InPar.UpdateKey), JD(Isim), 'JD for middle of exposure');
    end
    
 
end


switch lower(InPar.OutType)
    case 'jd'
        % do nothing
    case 'mjd'
        JD = JD - 2400000.5;
    otherwise
        error('Unknown OutType option');
end

        

        
        
