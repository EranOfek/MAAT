function [CumN,N]=src_number_count(Mag,Band)
% The cumulative number of sources in the GALEX-NUV band at high Gal. lat.
% Package: VO.GALEX
% Description: A simplistic, order of magnitude estimate of the surfae area
%              of sources in the GALEX-NUV band at high Galactic latitude
%              based on interpolation of the Bianchi et al. (2007) data in
%              the 18 to 22 mag range.
% Input  : - AB mag
%          - Band. Default is 'NUV'.
% Output : - Surface number of stars brighter than the requested magnitude
%            per deg^2 at high Galactic latitude.
%          - Differential number od stars per deg^2 per mag bin.
%     By : Eran O. Ofek                    Oct 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: N=VO.GALEX.src_number_count(21)
% Reliable: 2
%--------------------------------------------------------------------------

error('Not clear to me if Bianchi et al. (2007) give diff or cum dist')

if (nargin==1)
    Band = 'NUV';
end


switch lower(Band)
    case 'nuv'
        % Bianchi et al. 2007 - NUV
        A = -5.0398;
        B = 0.36;
        N = 10.^(A + B.*Mag);
        CumN = (10.^A.*(10.^(Mag.*B) - 1))./(B.*log(10));

    otherwise
        error('Unknown Band option');
end