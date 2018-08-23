function [Bias1,Res]=flag_bias_variable(Bias1,Bias2,varargin)
% Flag pixels with slowly variable bias levels
% Package: @SIM
% Description: Flag pixels with slowly variable bias levels
% Input  : - A SIM object containing bias images and their std images (in
%            the 'ErrIm' field).
%          - If empty, then assumes that the first input argument contains
%            multiple bias images of the same CCD and flag variable pixels.
%            Alternatively, this can be a SIM object with size identical to
%            that of the first input argument. In this case, each SIM
%            elemement corresponds to different CCD, and search for
%            variability between corresponding elemnts of the first and
%            second input arguments.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
% Output : - A SIM object containing the first input argument, in which the
%            MASK image field is updated with the variable bias pixels.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Nov 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% [Bias1,Res]=flag_bias_variable(Bias1,[],'Nbias1',16,'Nbias2',16);
% Reliable: 2
%--------------------------------------------------------------------------

ImageField   = SIM.ImageField;
ErrField     = SIM.ErrField;


DefV.VarThreshSigma       = 8;
DefV.Bit_BiasVariable     = 'Bit_BiasVariable';
DefV.Nbias1               = [];
DefV.Nbias2               = [];
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if (nargin==1)
    Bias2 = [];
end

if (isempty(Bias2))
    % search for variability among the many images of the same CCD
    if (numel(Bias1)<2)
        error('Must supply at least 2 bias images');
    end
    
    error('This option [Bias2 is not supplied] doesnt work yet');
else
    Nim1 = numel(Bias1);
    Nim2 = numel(Bias2);
    if (Nim1~=Nim2)
        error('Number of elements in Bias1 and Bias2 must be identical');
    end
    
    if (numel(InPar.Nbias1)==1)
        InPar.Nbias1 = InPar.Nbias1.*ones(Nim1);
    end
    if (numel(InPar.Nbias2)==1)
        InPar.Nbias2 = InPar.Nbias2.*ones(Nim1);
    end
    
    Res = Util.struct.struct_def({'Max','Min','Std','Rstd','Nflag'},size(Bias1));
    for Iim=1:1:Nim1
        
        SqN1 = sqrt(InPar.Nbias1(Iim)) - 1;
        SqN2 = sqrt(InPar.Nbias2(Iim)) - 1;
        
        DiffIm = Bias1(Iim).(ImageField) - Bias2(IndCCD,Iamp).(ImageField);
        ErrIm  = sqrt( (Bias1(Iim).(ErrField)./SqN1).^2 + (Bias2(Iim).(ErrField)./SqN2).^2 );
        Ratio = DiffIm./ErrIm;
        
        Flag = abs(Ratio)>InPar.VarThreshSigma;
        
        Res(Iim).Max  = max(Ratio(:));
        Res(Iim).Min  = min(Ratio(:));
        Res(Iim).Std  = std(Ratio(:));
        Res(Iim).Rstd = Util.stat.rstd(Ratio(:));
        Res(Iim).Nflag= sum(Flag(:));
        
    end
end
