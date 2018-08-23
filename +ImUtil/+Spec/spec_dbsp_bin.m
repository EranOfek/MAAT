function Bin=spec_dbsp_bin(KeyVal,Axis)
%--------------------------------------------------------------------------
% spec_dbsp_bin function                                            ImSpec
% Description: 
% Input  : - A string containinng the DBSP binning value
%            (i.e., the value of the 'CCDSUM' in the DBSP header).
%          - Axis direction {'x'|'y'|'xy'}.
% Output : - Binning in the requested direction.
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Mar 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Bin=spec_dbsp_bin('1 1','xy');
% Reliable: 2
%--------------------------------------------------------------------------

[BinAll]=sscanf(KeyVal,'%d %d');
switch lower(Axis)
    case 'x'
        Bin = BinAll(1);
    case 'y'
        Bin = BinAll(2);
    case 'xy'
        Bin = BinAll.';
    otherwise
        error('Unknown Axis option');
end
        

