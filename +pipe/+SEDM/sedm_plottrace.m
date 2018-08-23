function sedm_plottrace(SegmentsInfo,Image,varargin);
%--------------------------------------------------------------------------
% sedm_plottrace function                                             SEDM
% Description: plot a spexcell trace/s on image.
% Input  : - SegmentsInfo structure array for all the traces you like to
%            plot.
%          - Image name. If empty or not provided then will attempt to
%            plot over current image in ds9.
%          * Arbitrary number of pairs of ...,key,val,.. arguments.
%            The following keywords are available:
%            'XstartOffset' - slope start point relative to segment start.
%                      Default is 30.
%            'XendOffset' - slope end point relative to segment end.
%                      Default is 130.
%            'Step'  - plotting step in ds9. Default is 30 pix.
%            'OffsetType' - Offset type {'measured','smooth','original'}.
%                      Default is 'measured'.
%            'SlopeType' - Slope type {'measured','smooth','original'}.
%                      Default is 'measured'.
% Output : null
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Aug 2013
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: sedm_plottrace(SegmentsInfo(1130));
% Reliable: 2
%--------------------------------------------------------------------------


% Example: sedm_plottrace(SegmentsInfo);

DefV.XstartOffset       = 30;   
DefV.XendOffset         = 130;
DefV.Step               = 30;
DefV.OffsetType         = 'measured';   % {'measured','smooth','original'}
DefV.SlopeType          = 'measured';   % {'measured','smooth','original'}

InPar = set_varargin_keyval(DefV,'y','use',varargin{:});

Def.Image = [];;
if (nargin==1),
    Image  = Def.Image;
end

if (~isempty(Image)),
   ds9_disp(Image);
end

Nseg = length(SegmentsInfo);
for Iseg=1:1:Nseg,
    Xplot = (SegmentsInfo(Iseg).MinX-InPar.XstartOffset:InPar.Step:SegmentsInfo(Iseg).MaxX+InPar.XendOffset).'-SegmentsInfo(Iseg).MeanX;
    Par1 = SegmentsInfo(Iseg).Par1;
    switch lower(InPar.OffsetType)
        case 'original'
            % use Par1 as is
        case 'measured'
            Par1(2) = Par1(2) + SegmentsInfo(Iseg).MeasuredOffset;
        case 'smooth'
            Par1(2) = Par1(2) + SegmentsInfo(Iseg).PredictedOffset;
        otherwise
            error('Unknown OffsetType option');
    end
    switch lower(InPar.SlopeType)
        case 'original'
            % use Par1 as is
        case 'measured'
            Par1(1) = Par1(1);
        case 'smooth'
            Par1(1) = SegmentsInfo(Iseg).PredictedSlope;
        otherwise
            error('Unknown OffsetType option');
    end
        
        
    Yplot1 = polyval(Par1,Xplot);
    ds9_plotregion(Xplot+SegmentsInfo(Iseg).MeanX,Yplot1,'Type','circle','Size',1,'color','red');
end