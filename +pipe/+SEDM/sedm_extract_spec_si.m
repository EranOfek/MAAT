function [Object]=sedm_extract_spec_si(SI,varargin)
%--------------------------------------------------------------------------
% sedm_extract_spec_si function                                       SEDM
% Description: 
% Input  : - 
%          * Arbitrary number of pairs or arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            --- Additional parameters
%            Any additional key,val, that are recognized by one of the
%            following programs:
% Output : - 
% Tested : Matlab R2013a
%     By : Eran O. Ofek                    Aug 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: 
% Reliable: 
%--------------------------------------------------------------------------

DefV.SpexValField = 'SpecMedianVal';

DefV.WaveRange    = [450 850];
DefV.SpecField    = 'FlatSpexSpecFit';
DefV.MeanFun      = @median;
DefV.MarkerSize   = 30;

DefV.AperRad      = 80;
DefV.AnnRad       = [110 200];
DefV.WaveInt      = (350:0.5:900).';
DefV.InterpMethod = 'linear';
DefV.MedFiltBS    = 20;
DefV.MinCorr      = 0.8; % minimum correlation in SI.Arc_BestCorr]
DefV.FlatSI       = [];
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

if (~isfield(SI,InPar.SpexValField)),
    SI = sedm_plot_seg_int(SI,'WaveRange',InPar.WaveRange,...
                              'SpecField',InPar.SpecField,...
                              'MeanFun',InPar.MeanFun,...
                              'Plot',false);
end

Val = [SI.(InPar.SpexValField)].';
Val(Val<0) = 0;

 H = scatter([SI.MeanX].',[SI.MeanY].',InPar.MarkerSize,Val,'filled');
 set(H,'Marker','p')
 colorbar
 
 % Select object to extract
 Iobj = 0;
 BT = 1;
 while (BT==1),
     fprintf('--- Select object to extract ---\n');
     fprintf('Use Left click to select object, other clicks to quit\n');
     
     [X,Y,BT] = ginput(1);
     %X=1114
     %Y=1154
     
     if (BT==1),
         BT = 0;
         % extract object
         Iobj = Iobj+1;
         % position of objects
         Object(Iobj).X = X;
         Object(Iobj).Y = Y;
         
         % Distance
         Dist = Util.Geom.plane_dist([SI.MeanX].',[SI.MeanY].',X,Y);
         IndAper = find(Dist<InPar.AperRad & [SI.Arc_BestCorr].'>InPar.MinCorr);
         IndAnn  = find(Dist>InPar.AnnRad(1) & Dist<=InPar.AnnRad(2) & [SI.Arc_BestCorr].'>InPar.MinCorr);
         
         % plot aperture
         hold on;
         H(Iobj)   = plot_ellipse([X,Y],[InPar.AperRad, InPar.AperRad],0,0,'k');
         Ha1(Iobj) = plot_ellipse([X,Y],[InPar.AnnRad(1), InPar.AnnRad(1)],0,0,'k');
         Ha2(Iobj) = plot_ellipse([X,Y],[InPar.AnnRad(2), InPar.AnnRad(2)],0,0,'k');
 
         Object(Iobj).NspexAper = numel(IndAper);
         Object(Iobj).NspexAnn  = numel(IndAnn);
         
         SpecMat   = zeros(Object(Iobj).NspexAper,numel(InPar.WaveInt));
         SpecMatX  = repmat([SI(IndAper).MeanX].',1,numel(InPar.WaveInt));
         SpecMatY  = repmat([SI(IndAper).MeanY].',1,numel(InPar.WaveInt));
         
%         Colors=generate_colors(Object(Iobj).NspexAper);
%figure(3);
         for Ia=1:1:Object(Iobj).NspexAper,
             
             NN = ~isnan(SI(IndAper(Ia)).WaveCalib);
             SpecMat(Ia,:) = interp1(SI(IndAper(Ia)).WaveCalib(NN),SI(IndAper(Ia)).(InPar.SpecField)(NN),...
                                     InPar.WaveInt,InPar.InterpMethod).';
             
%plot(SI(IndAper(Ia)).WaveCalib(NN),SI(IndAper(Ia)).(InPar.SpecField)(NN),'k-','Color',Colors(Ia,:))
%             hold on;
             
         end
         Hspec = [ones(Object(Iobj).NspexAper,1), SpecMatX(:,1), SpecMatY(:,1), SpecMatX(:,1).*SpecMatY(:,1)];
         
         
         
         
         BackMat = zeros(Object(Iobj).NspexAnn,numel(InPar.WaveInt));
         BackMatX  = repmat([SI(IndAnn).MeanX].',1,numel(InPar.WaveInt));
         BackMatY  = repmat([SI(IndAnn).MeanY].',1,numel(InPar.WaveInt));

         
         for Ia=1:1:Object(Iobj).NspexAnn,
             
             NN = ~isnan(SI(IndAnn(Ia)).WaveCalib);
             BackMat(Ia,:) = interp1(SI(IndAnn(Ia)).WaveCalib(NN),SI(IndAnn(Ia)).(InPar.SpecField)(NN),...
                                     InPar.WaveInt,InPar.InterpMethod).';
                                 
             %plot(SI(IndAnn(Ia)).WaveCalib(NN),SI(IndAnn(Ia)).(InPar.SpecField)(NN),'k-','Color',Colors(Ia,:))
             %hold on;
         end
         % fit a linear surface to background values in annulus
         % fit model: B = B0 + alpha*X + beta*Y + gamma*X*Y
         % fit per wavelength, simultanously for all wavelengths
         Hback = [ones(Object(Iobj).NspexAnn,1), BackMatX(:,1), BackMatY(:,1), BackMatX(:,1).*BackMatY(:,1)];
         ParBack = Hback\BackMat; 
         
         % subtract background from spectra in aperture
         BS_SpecMat = SpecMat - medfilt1(Hspec*ParBack,InPar.MedFiltBS,[],2);
         Object(Iobj).SpecAperBS = sum(BS_SpecMat,1);
         
         
         SpecMatNorm = bsxfun(@times,SpecMat,1./SpecMat(1,:));
         FlagW = InPar.WaveInt>InPar.WaveRange(1) & InPar.WaveInt<InPar.WaveRange(2);
         W  = nanmedian(SpecMatNorm(:,FlagW),2);   % weight per spexcell
         Object(Iobj).Wave     = InPar.WaveInt;
         Object(Iobj).SpecAperW = sum(bsxfun(@times,SpecMat,W),1)./sum(W);
         Object(Iobj).SpecAperS = sum(SpecMat,1);
         Object(Iobj).SpecAnnS  = sum(BackMat,1);
         
         
         
         
     end
 end
 
 
 
 
                                   
