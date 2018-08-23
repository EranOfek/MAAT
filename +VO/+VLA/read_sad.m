function [Src,Info]=read_sad(File)
% Read AIPS SAD files.
% Package: VO.VLA
% Description: Read AIPS SAD files.
% Input  : - SAD file name
% Output : - Structure array with the SAD information per source.
%            .H - Flag indicating high point in the residual.
%            .L - Flag indicating low point in the residual.
%            .S - Flag indicating large RMS.
%          - Info.RA, Info.Dec : reference position 
%            Info.RMS.
%            All fluxes are given in mJy
% Tested : Matlab 7.11
%     By : Eran O. Ofek                    Jul 2011
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reliable: 2
%------------------------------------------------------------------------------

FID = fopen(File,'r');
StopReading = 0;
Counter = 0;
while (~feof(FID) && StopReading==0)
    Counter = Counter + 1;
    Line = fgetl(FID);
    if (length(Line)<6)
        % skip
    else
       SourceIndex = Util.string.str2num_nan(Line(1:5));
       if (~isnan(SourceIndex))
           if (strcmp(Line(7),'*'))
               Factor2 = 1000;
           else
               Factor2 = 1;
           end
           
           if (strcmp(Line(7),'L'))
              Src(SourceIndex).L = 1;
           else
              Src(SourceIndex).L = 0;
           end
           if (strcmp(Line(7),'H'))
              Src(SourceIndex).H = 1;
           else
              Src(SourceIndex).H = 0;
           end   
           if (strcmp(Line(7),'S'))
              Src(SourceIndex).S = 1;
           else
              Src(SourceIndex).S = 0;
           end             
           Src(SourceIndex).Peak  = Util.string.str2num_nan(Line(8:16)).*Factor.*Factor2;
           Src(SourceIndex).Dpeak = Util.string.str2num_nan(Line(18:24)).*Factor.*Factor2;
           Src(SourceIndex).Flux  = Util.string.str2num_nan(Line(26:34)).*Factor.*Factor2;
           Src(SourceIndex).Dflux = Util.string.str2num_nan(Line(36:42)).*Factor.*Factor2;
           SrcStrRA  = Line(45:57);
           SrcStrDec = Line(59:71);
           Src(SourceIndex).RA     = celestial.coo.convertdms(SrcStrRA,'SHb','r');
           Src(SourceIndex).Dec    = celestial.coo.convertdms(SrcStrDec,'SDb','R');
           Src(SourceIndex).Dx     = Util.string.str2num_nan(Line(73:80));
           Src(SourceIndex).Dy     = Util.string.str2num_nan(Line(82:88));
           Src(SourceIndex).Maj    = Util.string.str2num_nan(Line(91:98));
           Src(SourceIndex).Min    = Util.string.str2num_nan(Line(100:106));
           Src(SourceIndex).PA     = Util.string.str2num_nan(Line(108:112));
           Src(SourceIndex).Dmaj   = Util.string.str2num_nan(Line(115:120));
           Src(SourceIndex).Dmin   = Util.string.str2num_nan(Line(122:127));
           Src(SourceIndex).Dpa    = Util.string.str2num_nan(Line(129:131));

   
        
       else
           % look for information
           if (~isempty(strfind(Line,'Component widths')))
               StopReading = 1;
           end 
           
           
           if (~isempty(strfind(Line,'Reference Center:')))
               % Reference center
               StrRA  = Line(19:31);
               StrDec = Line (33:45);
               Info.RA     = celestial.coo.convertdms(StrRA,'SHb','r');
               Info.Dec    = celestial.coo.convertdms(StrDec,'SDb','R');
           
           end
        
           if (~isempty(strfind(Line,'Errors determined by theory from RMS')))
               % Reference center
               RMS = Util.string.str2num_nan(Line(37:43));
               UnitsR = Line(45:end);
               switch lower(UnitsR)
                   case 'microjy'
                       FactorR = 1e-3;
                   otherwise
                       error('Unknwon Units option');
               end
               Info.RMS = RMS.*FactorR;
            
           end
                
            if (~isempty(strfind(Line,'Fluxes expressed in units of')))
               % Reference center
               Units  = Line(30:end);
               switch lower(Units)
                   case 'microjy/beam'
                       Factor = 1e-3;
                   otherwise
                       error('Unknwon Units option');
               end
            end
        end 
    end
end
         
         


fclose(FID);
