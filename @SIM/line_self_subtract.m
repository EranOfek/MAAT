function [Sim,Line]=line_self_subtract(Sim,varargin)
% Subtract the mean of each line/row from itself in images in a SIM object.
% Package: @SIM
% Description: Subtract from each line or row in a SIM image its mean
%              value. This may be useful when the bias subtraction is not
%              perfect and some line to line variations are still visible.
% Input  : - A SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'ExecField' - The field on which to execute the self line
%                          subtraction. Default is 'Im'.
%            'Dim'       - The dimension of the self subtraction.
%                          1 will operate on lines (rows), while 2 will
%                          operate on columns. Default is 1.
%            'LineBackMethod' - Function handle for the method by which to
%                          calculate the line mean value.
%                          Fun(Matrix,Dim,LineBackMethodPar{:})
%                          Default is @rmean.
%            'LineBackMethodPar' - A cell array of additional parameters to
%                          pass to the 'LineBackMethod' function.
%                          Default is {[0.15 0.15]}.
%            'AddMeanLine' - A logical flag indicating if to add back to
%                          the image the global mean of the subtracted
%                          lines. Default is true.
% Output : - A SIM object with the self subtracted lines images.
%          - A SIM object of the subtracted lines.
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    Jul 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=line_self_subtrcat(Sim);
% Reliable: 2
%--------------------------------------------------------------------------

ImageField     = 'Im';

DefV.ExecField          = ImageField;
DefV.Dim                = 1;
DefV.LineBackMethod     = @rmean;
DefV.LineBackMethodPar  = {[0.15 0.15]};
DefV.AddMeanLine        = true;
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

VecOtherDim = [2 1];
OtherDim    = VecOtherDim(InPar.Dim);

Nsim = numel(Sim);
Line = SIM(size(Sim));
for Isim=1:1:Nsim
    if (isa(InPar.LineBackMethod,'function_handle'))
        Line(Isim).(InPar.ExecField) = InPar.LineBackMethod(Sim(Isim).(InPar.ExecField),...
                                                             OtherDim,...
                                                             InPar.LineBackMethodPar{:});
    else
        error('Unknown LineBackMethod option');
    end
    if (InPar.AddMeanLine)
        MeanLine = mean(Line(Isim).(InPar.ExecField));
    else
        MeanLine = 0;
    end
    
    Sim(Isim).(InPar.ExecField) = bsxfun(@minus,Sim(Isim).(InPar.ExecField),Line(Isim).(InPar.ExecField)) + MeanLine;
    
end

            