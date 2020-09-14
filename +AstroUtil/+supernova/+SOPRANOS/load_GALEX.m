function [jd, mjd, cps, cpserr, flux, fluxerr, prevRef, percRef ,RA, Dec, DiscMJD, z ,...
    mInd, mAB, mERR, lInd, limitmAB] = load_GALEX(sn_name)
% load GALEX data from the GALEX/PTF experiment
% Package: AstroUtil.supernove.SOPRANOS
% Description: load GALEX data from the GALEX/PTF experiment
% Input  : - sn_name
% Output : - Table with file contents
%               
% See also: AstroUtil.supernova.SOPRANOS.calcGrid
% Tested : Matlab 9.5
%     By : Noam Ganot                      Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% AstroUtil.supernova.SOPRANOS.load_GALEX('PTF12gnt');
% Reliable: 2
%--------------------------------------------------------------------------

load InfoAll.mat

ind = find(strcmp([Info.Name],sn_name(4:end)));
if isempty(ind) && strcmp(sn_name,'12sim')
    load InfoTest.mat
    ind = 1;
end
    

N = size(Info(ind(1)).LC,1);
info = Info(ind(1));
for I=ind(2:end)
    info.LC = [info.LC;Info(I).LC];
end


RA = info.RA;
Dec = info.Dec;
DiscMJD = convert.time(info.DiscJD,'jd','mjd');
z = info.z;

% http://galexgi.gsfc.nasa.gov/docs/galex/FAQ/counts_background.html
% NUV: Flux [erg sec-1 cm-2 A-1] = 2.06 x 10^-16 x CPS
% NUV: mAB = -2.5 x log10(CPS) + 20.08
GALEXF = 2.06e-16; %[erg sec-1 cm-2 A-1]/COUNT

MJD = convert.time(info.LC(:,1),'jd','mjd');
[MJD,ind] = sort(MJD);
cps = info.LC(ind,9)./info.LC(ind,4);      % CPS
cpserr = info.LC(ind,10)./info.LC(ind,4);  % CPS
jd = info.LC(ind,1);

prevRef = mean(info.LC(:,19)./info.LC(:,25))*GALEXF;
percRef = prctile(info.LC(:,9)./info.LC(:,4),25)*GALEXF;

ind = (cps~=inf);
jd = jd(ind);
mjd = MJD(ind);
flux = cps(ind)*GALEXF;  % [erg sec-1 cm-2 A-1]
fluxerr = cpserr(ind)*GALEXF;  % [erg sec-1 cm-2 A-1]
cps = cps(ind);
cpserr = cpserr(ind);

if (nargout > 12)
    mInd = (cps>0);
    mAB  = -2.5 * log10(cps(mInd)) + 20.08;
%     mERR = 2.5 * log10(1+cpserr(ind)./cps(ind));
    mERR = 2.5/log(10)*cpserr(mInd)./cps(mInd);
end

if (nargout > 15)
    lInd = (cps<=0);
    limitmAB = -2.5 * log10(cps(lInd) + 3*cpserr(lInd)) + 20.08;
end

end