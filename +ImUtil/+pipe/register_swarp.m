function [FileNameResamp,R,Sim,FlagGA]=register_swarp(Images,RA,Dec,varargin)
% Register images using astrometry.m and SWarp.
% Package: ImUtil
% Description: Register images using astrometry.m and SWarp.
% Input  : - Cell arrray containing list of images to swarp.
%          - J2000 R.A. [radians]
%          - J2000 Dec. [radians]
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - 
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Jun 2018
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

ImageField  = SIM.ImageField;
HeaderField = SIM.HeaderField;

DefV.BackType             = 'MANUAL'; % ['AUTO' | 'MANUAL']
DefV.GainCorrect          = true;
DefV.BackDefault          = 0;
DefV.BackSize             = 128;
DefV.SimResamp            = true;
DefV.deleteTmp            = true;
DefV.PixScale             = 1.0;
DefV.InterpMethod         = 'LANCZOS2';
DefV.OverSampling         = [1 1];
DefV.ImageSize            = [1000 1000];
DefV.read2simPar          = {};
DefV.ExecAstrometry       = true;
DefV.AstrometryPar        = {};
DefV.Flip                 = [-1 -1; 1 -1];

InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (SIM.issim(Images))
    Sim = Images;
    clear Images;
else
    Sim = FITS.read2sim(Images,InPar.read2simPar{:});
end
Nsim = numel(Sim);

if (InPar.GainCorrect)
    Sim = gain_correct(Sim);
end

% execute astrometry.m
FlagGA = true(Nsim,1);
if (InPar.ExecAstrometry)
    
    for Isim=1:1:Nsim
        try
            [R(Isim),Sim(Isim)] = astrometry(Sim(Isim),'Flip',InPar.Flip,InPar.AstrometryPar{:},'RA',RA,'Dec',Dec);
        catch
            FlagGA(Isim) = false;
        end
    end
    R   = R(FlagGA);
    Sim = Sim(FlagGA);
else
    R = [];
end


% save SimA as fits images
Nsim     = numel(Sim);
FileName = cell(1,Nsim);
for Isim=1:1:Nsim
    FileName{Isim}       = sprintf('Tmp_%05d.fits',Isim);
    FileNameResamp{Isim} = sprintf('Tmp_%05d.resamp.fits',Isim);
    [Flag,HeaderInfo] = FITS.write_old(Sim(Isim).(ImageField),FileName{Isim},Sim(Isim).(HeaderField));
end

ImUtil.Im.swarp(FileName,'RA',RA,'Dec',Dec,...
                         'SWarpPars',{'COMBINE','N',...
                                    'BACK_TYPE',InPar.BackType,...
                                    'BACK_DEFAULT',InPar.BackDefault,...
                                    'BACK_SIZE',InPar.BackSize,...
                                    'PIXELSCALE_TYPE','MANUAL',...
                                    'PIXEL_SCALE',sprintf('%f',InPar.PixScale),...
                                    'RESAMPLING_TYPE',InPar.InterpMethod,...
                                    'OVERSAMPLING',sprintf('%d,%d',InPar.OverSampling),...
                                    'COMBINE_TYPE','WEIGHTED',...
                                    'IMAGE_SIZE',sprintf('%d,%d',InPar.ImageSize)});
                                
if (InPar.SimResamp)
    Sim = FITS.read2sim(FileNameResamp);
end

if (InPar.deleteTmp)
    delete('Tmp_*.fits');
end
