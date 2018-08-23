function [SelectedMaxima,Contrast_UnitsRMS]=find_contrast_peaks(X,Y,ContrastRMS,varargin)
%--------------------------------------------------------------------------
% find_contrast_peaks function                                      ImSpec
% Description: Given a list of [X, Y], find the position of local maxima
%              in the list which are above a certain contrast above the
%              local rms.
%              The rms is defined as the 68 percentile of the extramum
%              Y-value distribution, while the contrast is defined
%              as the Y-value of the local maxima subtract from the
%              preceeding local minima.
% Input  : - X
%          - Y
%          - Threshold contrast in units of the rms (default is 4).
%          * Arbitrary number of pairs of arguments ...,keyword,value,...
%            The followinf keywords are available:
%            'PeakAlgo' - Algorithm to use:
%                         'RMSmin' - RMS is estimated from distribution
%                                    of minima, and contrast from diference
%                                    between local maxima and the next minima.
%                                    Good when the data contains many nearby
%                                    peaks. Default.
%                          'win'   - RMS is calculating using stdfilt1.m.
%                          'mf'    - Matched filter using matched_filter.m.
%            'stdpar'   - Cell array of parameters to pass to pass to
%                         stdfilt1.m. Default is {30,[],0.68}.
%            'mfpar'    - Cell array of parameters to pass to
%                         matched_filter.m. {Filter,Subtract,Thresh}.
%                         Default is {'median'}.
%            'RemNaN'   - Remove NaNs from X/Y list.
%                         {'y'|'n'}, default is 'y'.
% Output : - Selected maxima above threshold [X, Y, 2nd derivative]
%          - Contrast of local maxima in units opf the rms.
% Tested : Matlab 7.10
%     By : Eran O. Ofek                    May 2010
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [SelectedMaxima,Contrast_UnitsRMS]=find_contrast_peaks(X,Y);
% Reliable: 2
%--------------------------------------------------------------------------

Def.ContrastRMS = 4;
if (nargin==2),
   ContrastRMS = Def.ContrastRMS;
end
if (isempty(ContrastRMS)),
   ContrastRMS  = Def.ContrastRMS;
end

DefV.PeakAlgo = 'RMSmin';
DefV.stdpar   = {30,[],0.68};
DefV.mfpar    = {'median'};
DefV.RemNaN   = 'y';
InPar = set_varargin_keyval(DefV,'n','use',varargin{:});

switch lower(InPar.RemNaN)
    case 'y'
       Innan = find(~isnan(X) & ~isnan(Y));
       X = X(Innan);
       Y = Y(Innan);
    otherwise
        % do nothing
end
SelectedMaxima    = zeros(0,3);
Contrast_UnitsRMS = [];
if (~isempty(X)),
    
    if ((X(end)-X(1))<0),
        X = flipud(X);
        Y = flipud(Y);
    end

    % find extramum

    Extram = find_local_extramum(X,Y);
    % remove points with 2nd derivative = 0
    I = abs(Extram(:,3))>0;
    Extram = Extram(I,:);
    % assuming each minima is followed by a maxima
    % estimate the rms

    %Y_ErrCL = err_cl(Y);
    %Inoise = find(Y<Y_ErrCL(1,2));
    %err_cl(Y(Inoise))


    switch lower(InPar.PeakAlgo)
         case 'rmsmin'
            ErrCL = err_cl(diff(Extram(:,2)));
            RMS   = ErrCL(1,2)./2.2;   % 2.2 is the factor required to transform the rms of peaks to rms.

            % find maxima/minima by looking for negative 2nd derivative
            Imax = find(Extram(:,3)<0);

            Extram = [Extram;[Inf Inf 0]];

            Contrast_UnitsRMS = (Extram(Imax,2) - Extram(Imax+1,2))./RMS;
            Ip = find(Contrast_UnitsRMS>ContrastRMS);
            SelectedMaxima = Extram(Imax(Ip),:);
            Contrast_UnitsRMS = Contrast_UnitsRMS(Ip);
         case 'win'
            Ystd = stdfilt1(Y,InPar.stdpar{:});
            Std_At_X = interp1(X,Ystd,Extram(:,1),'linear');
            Contrast_UnitsRMS = Extram(:,2)./Std_At_X;

            Ip = find(Contrast_UnitsRMS>ContrastRMS & Extram(:,3)<0);
            SelectedMaxima = Extram(Ip,:);
            Contrast_UnitsRMS = Contrast_UnitsRMS(Ip);

         case 'mf'
            % matched filter
            error('mf option doesnt exist yet');
            [Peaks,Y,Thresh]=matched_filter(X,InPar.mfpar{:}); % Filter,Subtract,Thresh
            
     otherwise
        error('Unknown PeakAlgo option');
    end

    %plot(X,Y);
    %hold on;
    %plot(Extram(Imax(Ip),1),Extram(Imax(Ip),2),'ko')
end