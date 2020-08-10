function [Coo,Lat]=xy2coo(Head,X,Y,OutType)
% Convert X/Y in WCS or HEAD object to celestial coordinates
% Package: @WorldCooSys
% Description: Given an WCS or HEAD object (or e.g., a SIM image, etc) read
%              the WCS from the header and convert X/Y to
%              Longitude/Latitude.
%              Support SIP and PV distortions.
% Input  : - An HEAD object (or SIM, AstCat).
%          - Vector of X coordinates [pixels].
%            [radian, sexagesimal string or vector] see convertdms.m
%          - Vector of Y coordinates [pixels].
%            [radian, sexagesimal string or vector] see convertdms.m
%          - Output type: 'astcat'|'struct'|'mat'.
%            Default is 'astcat'. However, if two ouput arguments are
%            requested than force ouput to 'mat'.
% Output : - Either an AstCat object or structure array with
%            the Long/Lat coordinates, or a vector of the Long
%            coordinates [radian].
%          - Vector of Lat [radian].
% Tested : Matlab R2014a
%     By : Eran O. Ofek, Tali Engel        Dec 2014
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: Coo = xy2coo(S,100,[200;100]);
%          [Long,Lat] = xy2coo(S(9),100,200);
% Reliable: 2
%--------------------------------------------------------------------------
import Util.struct.*

Def.OutType = 'astcat';    % 'mat'|'struct'|'astcat' 
if (nargin==3)
    OutType = Def.OutType;
elseif (nargin==4)
    % do nothing
else
    error('Illegal number of input arguments');
end


if (nargout>1)
    % if two output arguments
    % then force output to 'mat'
    OutType = 'mat';
end

% Get WCS from Header
if (HEAD.ishead(Head))
    % An HEAD object
    if (any(isempty_wcs(Head)))
        % The WorldCooSys is HEAD is empty
        Head = populate_wcs(Head);
    end
end
WCS = Head;

% define output
Nwcs = numel(WCS);
switch lower(OutType)
    case 'astcat'
        Coo = AstCat(size(WCS));
    case 'struct'
        Coo = struct_def({'Cat','ColCell'},Nwcs,1);
    otherwise
        % do nothing
end
    
% make sure that X/Y have the same size
if (numel(X)~=numel(Y))
    if (numel(X)==1)
        X = X.*ones(size(Y));
    elseif (numel(Y)==1)
        Y  = Y.*ones(size(X));
    else
        error('Axes have different lengths');
    end
end


% for each Header
for Iwcs=1:1:Nwcs
    
    %if (isnan(WCS(Iwcs).WCS.CTYPE1))
    if isnan(WCS(Iwcs).WCS.CTYPE{1})
        % no valid WCS in header
        Long = nan(size(X));
        Lat  = nan(size(Y));
    else    
        ProjType1 = WCS(Iwcs).WCS.CTYPE{1}(6:8);  % see also read_ctype.m
        ProjType2 = WCS(Iwcs).WCS.CTYPE{2}(6:8);
        if (~strcmp(ProjType1,ProjType2))
            error('Axes have different orojection types');
        end

        RAD = 180./pi;
        % transformation
        % relevant for both projection types:
        if (~strcmp(WCS(Iwcs).WCS.CUNIT{1},WCS(Iwcs).WCS.CUNIT{2}))
            error('CUNIT1 must be identical to CUNIT2');
        end

        switch lower(WCS(Iwcs).WCS.CUNIT{1})
            case {'deg','degree'}
                Factor = RAD;
            case {'rad','radian'}
                Factor = 1;
            otherwise
                error('Unknown CUNIT option');
        end

        U = X - WCS(Iwcs).WCS.CRPIX(1);
        V = Y - WCS(Iwcs).WCS.CRPIX(2);

        % Deal with SIP distortions
        if (isfield_notempty(WCS(Iwcs).WCS,'sip'))
            %UV = [U.' V.'];
            if (isstruct(WCS(Iwcs).WCS.sip))
                Ap     = WCS(Iwcs).WCS.sip.Ap;
                Aq     = WCS(Iwcs).WCS.sip.Aq;
                Acoeff = WCS(Iwcs).WCS.sip.Acoeff;
                Bp     = WCS(Iwcs).WCS.sip.Bp;
                Bq     = WCS(Iwcs).WCS.sip.Bq;
                Bcoeff = WCS(Iwcs).WCS.sip.Bcoeff;

                F = zeros(length(X),1);
                G = zeros(length(X),1);

                % the distortion correction is different for each point (X,Y):
                for I = 1:length(X)
                    F(I) = sum(Acoeff.*(U(I).^Ap).*(V(I).^Aq));
                    G(I) = sum(Bcoeff.*(U(I).^Bp).*(V(I).^Bq));
                end

                U = U + F;
                V = V + G;
            end
        end

        switch lower(ProjType1)
            case {'tan','tpv'}
                %[Long,Lat] = xy2sky_tan(WCS,X,Y,HDUnum);
                if (size(U,2)==1)
                    U = U.';
                end
                if (size(V,2)==1)
                    V = V.';
                end
                %D = WCS.CD*[X - WCS.CRPIX1; Y - WCS.CRPIX2]; % + [WCS.CRVAL1; WCS.CRVAL2]
                D = WCS(Iwcs).WCS.CD*[U;V];  
                D = D./Factor;
                DX = D(1,:);
                DY = D(2,:);
                [Long,Lat] = celestial.proj.pr_ignomonic(DX.',DY.',[WCS(Iwcs).WCS.CRVAL(1), WCS(Iwcs).WCS.CRVAL(2)]./Factor);

            case 'ait'
                %[Long,Lat] = xy2sky_ait(WCS,X,Y,HDUnum);
                if (WCS(Iwcs).WCS.CRVAL(1)~=0 || WCS(Iwcs).WCS.CRVAL(2)~=0)
                    error('This version does not support CRVAL ne 1');
                end
                %X     = (X - WCS.CRPIX1).*WCS.CDELT1;
                %Y     = (Y - WCS.CRPIX2).*WCS.CDELT2;
                X = U.*WCS(Iwcs).WCS.CDELT(1);
                Y = V.*WCS(Iwcs).WCS.CDELT(2);
                [Long,Lat] = celestial.proj.pr_ihammer_aitoff(X,Y,Factor);

            otherwise
                error('Unsupported projection type: %s',ProjType1);
        end

    end
    
     switch lower(OutType)
        case 'astcat'
            Coo(Iwcs).Cat      = [Long,Lat];
            Coo(Iwcs).ColCell  = {'Long','Lat'};
            Coo(Iwcs).ColUnits = {'rad','rad'};
            Coo                = colcell2col(Coo);
        case 'struct'
            Coo(Iwcs).Cat     = [Long, Lat];
            Coo(Iwcs).ColCell = {'Long','Lat'};
        case 'mat'
            if (nargout<2)
                Coo = [Long, Lat];
            else
                Coo = Long;
                %Lat = Lat;
            end
                
            if (Iwcs>1)
                warning('mat output is requested for multi element object - mat contains data from last object only');
            end
        otherwise
            error('Unknown OutType option');
    end
    
end

