function Sim=interp2_nan(Sim,varargin)
% Interpolate over NaNs in a SIM array.
% Package: @SIM
% Description: Interpolate over NaNs in a SIM array.
%              Use the 3rd party function inpaint_nans.m.
% Input  : - A SIM object.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Method'     - Interpolation method. See inpaint_nans.m for
%                           options. Default is 0.
%            'ExecField'  - A cell array (or a single string) of field
%                           names on which to execute the function.
%                           Default is {'Im'}.
%            'CCDSEC'     - CCD section on which to execute the function.
%                           This is either a string indicating the header
%                           keyword containing the CCDSEC or a vector of
%                           [Xmin, Xmax, Ymin, Ymax].
%                           Note that this option do not trim the image.
%                           Use SIM/trim_image.m for trimming.
%            'ReplaceKey' - A three column cell array of {key,val,comment},
%                           or an HEAD object which keywords to replace
%                           or add to the SIM header.
%            'AddKey'     - Like 'ReplaceKey', but adding keywords,
%                           without replacment.
%            'MaskFun'    - Function that sets the bit mask:
%                           Fun(Sim1,MaskFunPar{:}).
%                           Default is empty.
%            'MaskFunPar' - Additional parameters to pass to MaskFun.
%                           Default is {}.
% Output : - A SIM object in which NaN values are replaced by interpolated
%            values.
% See also: inpaint_nans.m
% License: GNU general public license version 3
% Tested : Matlab R2015b
%     By : Eran O. Ofek                    May 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Sim=interp2_nan(Sim);
% Reliable: 2
%--------------------------------------------------------------------------


DefV.Method             = 0;
DefV.ExecField          = {'Im'};
DefV.CCDSEC             = [];
DefV.ReplaceKey         = {};
DefV.AddKey             = {};
DefV.MaskFun            = [];
DefV.MaskFunPar         = {};
%InPar = set_varargin_keyval(DefV,'n','use',varargin{:});
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

Sim = ufun2sim(Sim,@inpaint_nans,'FunAddPar',  {InPar.Method},...
                                 'ExecField',  InPar.ExecField,...
                                 'CCDSEC',     InPar.CCDSEC,...
                                 'ReplaceKey', InPar.ReplaceKey,...
                                 'AddKey',     InPar.AddKey,...
                                 'MaskFun',    InPar.MaskFun,...
                                 'MaskFunPar', InPar.MaskFunPar);
                       
                       