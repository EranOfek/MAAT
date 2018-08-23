function [Sim,Flag]=flag_bad_columns(Sim,varargin)
% SHORT DESCRIPTION HERE
% Package: @SIM
% Description: 
% Input  : - 
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jul 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

ImageField = SIM.ImageField;

DefV.ImageOp             = 'nan';
DefV.Dim                 = 1;
DefV.Thresh              = 10;
DefV.MedFiltBlock        = 100;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


Nsim = numel(Sim);
for Isim=1:1:Nsim
    
    Size = size(Sim(Isim).(ImageField));
    
    % global median
    GlobalMed   = nanmedian(Sim(Isim).(ImageField)(:));
    
    % collapse by X/Y-axis
    LineMed     = nanmedian(Sim(Isim).(ImageField),InPar.Dim);
    
    LineMedFilt = medfilt1(LineMed,InPar.MedFiltBlock);
    LineMedBS   = LineMed - LineMedFilt;
    
    LineStd     = Util.stat.rstd(LineMedBS(:));
    
    z = LineMedBS./LineStd;
    
    BadColFlag = z<(-InPar.Thresh);
    
    switch lower(InPar.ImageOp)
        case 'nan'
            % Replace bad columns with NaNs
            
            if (InPar.Dim==1)
                Sim(Isim).(ImageField)(:,BadColFlag) = NaN;
            else
                Sim(Isim).(ImageField)(BadColFlag,:) = NaN;
            end
        case 'replace'
            % Replace bad colum by median value
            
            
            if (InPar.Dim==1)
                Sim(Isim).(ImageField)(:,BadColFlag) = ones(Size(1),1)*LineMedFilt(BadColFlag);
            else
                Sim(Isim).(ImageField)(BadColFlag,:) = LineMedFilt(BadColFlag)*ones(Size(2),1);
            end
            
        case 'interp'
            % Interpolate over bad columns
            
        otherwise
            error('Unknown ImageOp option');
    end
    
    
    Flag(Isim).BadColFlag = BadColFlag;
end

    