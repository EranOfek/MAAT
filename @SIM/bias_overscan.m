function [Sim,BiasSim]=bias_overscan(Sim,varargin)
% Read the overscan bias region and subtract it from image. 
% Package: @SIM
% Description: Read the overscan bias region and subtract it from image. 
% Input  : - A SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField' - A cell array (or a single string) of field
%                          names on which to execute the function.
%                          Default is {'Im'}.
%            'OverSec'   - Overscan section region.
%                          Either a string containing an header keyword
%                          name containing the bias overscan region,
%                          or [Xmin Xmax Ymin Ymax].
%                          Default is {'BIASSEC','OVERSEC'}.
%            'FinalSec'  - The final image section region.
%                          This should be used if the bias overscan region
%                          should be trimmed.
%                          This is either a string containing an header
%                          keyword name containing the bias overscan
%                          region, or [Xmin Xmax Ymin Ymax], or [].
%                          If empty then do not trim image.
%                          Default is 'CCDSEC'.
%            'OverSecDim'- The dimension over which to collapse the
%                          overscan region: 'x'|'y'|'auto'.
%                          The 'auto' option will attempt to find the
%                          right dimension by using the longest dimension.
%                          Default is 'auto'.
%            'OverSecMethod' - Function handle for the collapse of the
%                          overscan region method.
%                          This function should be of the form:
%                          Fun(Matrix,Dim,OverSePar{:}).
%                          E.g., @nanmean, @nanmedian, @rmean.
%                          Default is @nanmean.
%            'OverSecPar' - Additional arguments to pass to the
%                          OverSecMethod function.
%                          Default is {}.
%            'SubtractBias'- Subtract bias overscan from image: true|false.
%                          Default is true.
%            'GetKeyMethod' - Method by which to get the keyword value out
%                          of multiple keywords. See getkey_fromlist.m
%                          for options. Default is 'first'.
% Output : - An overscan bias subtracted SIM object.
%          - A SIM object containing the collapsed overscan bias for
%            each image.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Apr 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Sim,BiasSim]=bias_overscan(Sim);
%          % subtract bias but do not trim final image
%          [Sim,BiasSim]=bias_overscan(Sim,'FinalSec',[]);
% Reliable: 2
%--------------------------------------------------------------------------



DefV.ExecField          = {'Im'};
DefV.OverSec            = {'BIASSEC','OVERSEC'};   % keyword | [Xmin Xmax Ymin Ymax]
DefV.FinalSec           = 'CCDSEC';    % key | [Xmin...] | []
DefV.OverSecDim         = 'auto';
DefV.OverSecMethod      = @nanmean; %@nanmedian;
DefV.OverSecPar         = {};
DefV.SubtractBias       = true;
DefV.GetKeyMethod       = 'first';
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


% try to read overscan region from header
OverSec = ccdsec(Sim,InPar.OverSec,true,InPar.GetKeyMethod);
% ccdsec returns a 4 column matrix - line per image


% make sure ExecField is a cell array
if (~iscell(InPar.ExecField))
    InPar.ExecField = {InPar.ExecField};
end

Nf   = numel(InPar.ExecField);
Nsim = numel(Sim);

if (nargout>1)
    BiasSim = SIM(size(Sim));
end
for Isim=1:1:Nsim
    % for each SIM element
    for If=1:1:Nf
        % for each field
        Sec = OverSec(Isim,:);
        if (any(isnan(Sec)))
            warning('Overscan region was not found for image %d - no correction applied',Isim);
        else
        
            % overscan region
            Scan = Sim(Isim).(InPar.ExecField{If})(Sec(3):Sec(4),Sec(1):Sec(2));
            switch lower(InPar.OverSecDim)
                case 'x'
                    Dim = 2;
                case 'y'
                    Dim = 1;
                case 'auto'
                    % select overscan dimension by its longest axis
                    % the flip is needed for exchanging [I,J]->[X,Y]
                    [~,Dim] = max(fliplr(size(Scan)));
                otherwise
                    error('Unknown OverSecDim option');
            end
            
            % collapse overscan region into a line
            BiasLine = InPar.OverSecMethod(Scan,Dim,InPar.OverSecPar{:});
            
            % populate BiasSim
            if (nargout>1)
                BiasSim(Isim).(InPar.ExecField{If}) = BiasLine;
            end
            
            % subtract bias line from image
            if (InPar.SubtractBias)
                % Check that the dimensions of the image and bias line are
                % consistent
                SizeIm = size(Sim(Isim).(InPar.ExecField{If}));
                if any(SizeIm==length(BiasLine))
                    % subtract bias
                    Sim(Isim).(InPar.ExecField{If}) = bsxfun(@minus,Sim(Isim).(InPar.ExecField{If}),BiasLine);
                else
                    error('In image %d - Discrepancy between image size and bias line length',Isim);
                end
            end
            
        end
        
    end  % over all fields
    
end  % over all images

% cut the final section
if (~isempty(InPar.FinalSec))
    % trim all images in SIM
    Sim = trim_image(Sim,InPar.FinalSec); %,'ExecField',InPar.ExecField);
end
    


