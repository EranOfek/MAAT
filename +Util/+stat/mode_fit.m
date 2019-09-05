function [Mode,StD]=mode_fit(Array,varargin)
% Estimate the mode of an array by fitting a Gaussian to its histogram.
% Package: @Util.stat
% Description: Estimate the mode of an array by fitting a Gaussian to
%              the histogram of the array around its median.
%              Return also the Sigma of the Gaussian fit.
%              If best fit is out of range then set the mode to the maximum
%              value of the histogram and std to sqrt of the max value.
% Input  : - Array for which to calculate the mode and StD.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ElementPerBin' - Typical number of points per histogram bin.
%            'TrimEdge2D' - If the input is a 2D array, trim the edges
%                           by this number of pixels. Default is 5.
%            'Percent'    - The percentile of the data for which to fit
%                           the mode [lower upper]. This means that
%                           only values between the lower and upper
%                           percentiles will be used.
%                           Default is [0.025 0.9].
%            'Nbad'       - If histogram difference contains jumps
%                           which are larger than this value then these
%                           jumps will be removed from the fit.
%                           INACTIVE.
%            'BinSize'    - Bin size. If empty then will attempt to
%                           calculate bin size. Default is empty.
%            'JoinOut'    - Joining the Mode and Std into the first output
%                           argument. Default is false.
% Output : - Mode.
%          - Fitted sigma.
% License: GNU general public license version 3
% Tested : Matlab R2014a
%     By : Eran O. Ofek                    Apr 2015
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Mode,StD]=Util.stat.mode_fit(randn(1000,1000))
% Reliable: 2
%--------------------------------------------------------------------------
import Util.fit.*

SIGMA1 = 0.6827;


DefV.ElementPerBin     = 100;
DefV.TrimEdge2D        = 5;
DefV.Percent           = [0.025 0.9];
DefV.JoinOut           = false;
DefV.Rem0              = true;
DefV.Nbad              = 100;
DefV.BinSize           = [];
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (isempty(Array))
    Mode = NaN;
    StD  = NaN;
else
    
    if (min(size(Array))>1 && InPar.TrimEdge2D>0)
        % Trim edges of 2D array to remove problems near image bounderies
        Array = Array(InPar.TrimEdge2D:1:end-InPar.TrimEdge2D,InPar.TrimEdge2D:1:end-InPar.TrimEdge2D);
    end
    Array = single(Array);
    
    if numel(unique(Array))==1
        Mode = Array(1);
        StD  = 0;
    else
        

        % remove NaNs

        % remove zeros
        %Array = Array(Array~=0);

        if (isempty(InPar.BinSize))
            Array      = sort(Array(:));
            %Array      = sort(Array(~isnan(Array)));
            Nel        = numel(Array);
            Lower      = Array(ceil(Nel.*InPar.Percent(1)));
            Upper      = Array(floor(Nel.*InPar.Percent(2)));

            %Lower      = prctile(Array(:),InPar.Percent(1).*100);
            %Upper      = prctile(Array(:),InPar.Percent(2).*100);
            %Percentile = err_cl(Array(:),InPar.Percent);
            BinSize = (Upper-Lower).*InPar.ElementPerBin./Nel;

            % there is no point to have resolution higher than
            % the minmimum possible error on the mean:
            MaxRes = 0.5.*sqrt((Upper-Lower).*0.5)./sqrt(Nel);
            BinSize = max(BinSize,MaxRes);

        else
            Lower      = min(Array(:));
            Upper      = max(Array(:));
            BinSize    = 1;

        end

        Edges   = (Lower-BinSize:BinSize:Upper+BinSize).';

        if (InPar.Rem0)
            Array = Array(Array(:)~=0);
        end

        % very slow
        % tic;
        % [N]=histc(Array(:),Edges);
        % toc


        %UA = unique(Array);
        %if (numel(UA)==1),
        if (range(Array)<(2.*BinSize))
            warning('Image contains a single unique pixel value');
            UA = unique(Array);
            Mode = UA;
            StD  = 0;
        else
            [N]=histcounts(Array(:),Edges)';

            Edges = Edges(1:end-1);
            Edges = Edges + 0.5.*(Edges(2)-Edges(1));
            %N = medfilt1(N,20);

            %FlagGood = abs(N(2:end)-N(1:end-1))<InPar.Nbad;
            %bar(Edges,N);
            %plot(Edges,N)
            %Res=fit_gauss1d(Edges(FlagGood),N(FlagGood),1);

            % remove highest point
            [MaxN] = max(N);
            Flag = N<MaxN;
            Edges = Edges(Flag);
            N     = N(Flag);
            Res=Util.fit.fit_gauss1d(Edges,N,1);
            if (Res.X0<min(Edges) || Res.X0>max(Edges))
                warning('Best fit mode is out of range - set mode to max value and std to rstd');
                [~,MaxI] = max(N);
                %Mode = Edges(MaxI);
                %StD  = sqrt(Mode);

                Mode = median(Array);
                StD  = Util.stat.rstd(Array);
            else
                %plot(Res.Resid)
                %input('hi')
                Mode = Res.X0;
                StD  = Res.Sigma;
            end
        end
    end
end


if (InPar.JoinOut)
    Mode = [Mode, StD];
end