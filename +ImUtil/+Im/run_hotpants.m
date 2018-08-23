function [Stat,Res]=run_hotpants(RefIm,NewIm,DiffIm,HPars)
% Run hotpants (image subtraction) from matlab.
% Package: ImUtil.Im
% Description: Run hotpants (image subtraction) from matlab.
% Input  : - FITS reference image.
%          - FITS new image to subtract from reference image.
%          - Output FITS difference image.
%          - Cell array of optional hotpants parameters.
%            These are pairs of keywords and values.
%            All values should be given in a string format.
%            See hotpans_parameters.txt for available keywords and defaults.
%            Default is an empty cell array.
%            Alternatively, this may be a string of all parameters.
%            e.g., '-rss 30 -tl -0.1 -il -0.1 -tu 1 -iu 1'.
% Output : - Status returned by the system command after running hotpants.
%          - Result returned by the system command after running hotpants.
% Tested : Matlab 7.13
%     By : Eran O. Ofek                    Jun 2012
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: ImUtil.Im.run_hotpants('ref.fits','new.fits','diff.fits','-rss 30 -tl -0.1 -il -0.1 -tu 1 -iu 1');
%          ImUtil.Im.run_hotpants('ref.fits','new.fits','diff.fits',{'rss','30','tl','-0.1','il','-0.1','tu','1','iu','1'});
%------------------------------------------------------------------------------
MSDir = Util.files.which_dir(mfilename);
Prog = sprintf('%s/../bin/Hotpans/hotpants_v5.1.10b/hotpants',MSDir);

Def.HPars = {};
if (nargin==3),
   HPars = Def.HPars;
elseif (nargin==4),
   % do nothing
else
   error('Illegal number of input arguments');
end

% constract arguments of hotpants parameters
if (ischar(HPars)),
   ParString = HPars;
elseif (iscell(HPars)),
   N = length(HPars);
   ParString = '';
   for I=1:2:N-1,
      ParString = sprintf('%s -%s %s',ParString,HPars{I},HPars{I+1});
   end
else
   error('HPars type is unseported');
end

% Example for GALEX image subtraction:
% hotpants -inim 1.fits -tmplim 2.fits -outim diff2.fits -rss 30 -tl -0.1 -il -0.1 -tu 1 -iu 1

% run hotpants
RunString = sprintf('%s -inim %s -tmplim %s -outim %s %s',Prog,NewIm,RefIm,DiffIm,ParString);
[Stat,Res] = system(RunString);
