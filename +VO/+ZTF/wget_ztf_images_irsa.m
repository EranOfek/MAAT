function [Files,Links,Prop,Cat]=wget_ztf_images_irsa(RA,Dec,varargin)
% Query and retrieve ZTF images from the IRSA archive
% Package: VO.ZTF
% Description: Query and retrieve ZTF images from the IRSA archive.
% Input  : - J2000.0 R.A. [radians, [H M S], or sexagesimal string], or
%            a string containing object name (e.g., 'm31').
%            If empty, then the query is done without position (only the
%            WHERE clause).
%          - J2000.0 R.A. [radians, [sign D M S], or sexagesimal string].
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'GetFiles'- Get files (true) or just query images (false).
%                        Default is true.
%            'GetN'    - Get specific image {'all','last','first'}.
%                        Default is 'all'.
%            'ImType'  - ZTF image type to query: 'sci' | 'raw' | 'cal'.
%                        Default is 'sci'.
%            'Product' - Image product type:
%                        'log'|'mask'|'image'|'scilog'|'sex'|'dao'|'daopsf'|'diff'|'diffpsf'|...
%                        Default is 'image'.
%            'Where'   - String containing WHERE clause (e.g.,
%                        'field=600 and ccdid=2'). Default is ''.
%            'Intersect'- Options are: 'COVERS' | 'ENCLOSED' | 'CENTER' |
%                        'OVERLAPS'.
%                        Default is 'OVERLAPS'.
%            'Size'    - Position search size (heiht, [width]) [deg].
%                        Default is 0.
%            'queryPar'- Cell array of additional key,val parameters to
%                        pass to VO.ZTF.irsa_query_ztf_images.
%            'constructPar'- Cell array of additional key,val parameters to
%                        pass to VO.ZTF.irsa_image_link.
%            'User'    - String containing the IRSA/IPAC user name, or a
%                        a cell array containing a file name of user/pass.
%                        Default is
%                        {'/home/eran/matlab/passwords/ztf_ipac_pass'}.
%            'Pass'    - String containing the IRSA/IPAC password,.
%                        Default is [].
%            'pwgetExtra'- Extra parameter to pass to www.pwget.
%                        Default is '--no-check-certificate --load-cookies=cookies.txt'.
%            'MaxGet'  - Max parallel wget to pass to www.pwget.
%                        Default is 15.
% Output : - Cell array of retrieved file names.
%          - Cell array of links.
%          - Structure array of image properties.
%          - AstCat object containing the queried table.
% License: GNU general public license version 3
%     By : Eran O. Ofek                    Nov 2017
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: [Files,Links,Prop,Cat]=VO.ZTF.wget_ztf_images_irsa(358./RAD,23./RAD);
%          Files = VO.ZTF.wget_ztf_images_irsa([],[],'Where','field=600 and ccdid=2','ImType','raw'); 
%          [Files,Links] = VO.ZTF.wget_ztf_images_irsa([],[],'Where','field=600 and ccdid=2','ImType','sci','Product','diff','GetFiles',false); 
% Reliable: 2
%--------------------------------------------------------------------------


DefV.GetFiles             = true;
DefV.GetN                 = 'all';
DefV.ImType               = 'sci';
DefV.Product              = 'image';
DefV.Where                = '';
DefV.Intersect            = 'OVERLAPS';
DefV.Size                 = 0;
DefV.QueryPar             = {};
DefV.constructPar         = {};
DefV.User                 = {'/home/eran/matlab/passwords/ztf_ipac_pass'}; 
DefV.Pass                 = [];
DefV.pwgetExtra           = '--no-check-certificate --load-cookies=cookies.txt';
DefV.MaxGet               = 15;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

% get user/pass from passwords file
if (iscell(InPar.User))
    [InPar.User,InPar.Pass]=Util.files.read_user_pass_file(InPar.User{1});
end

Cat=VO.ZTF.irsa_query_ztf_images(RA,Dec,'Where',InPar.Where,'Intersect',InPar.Intersect,'Size',InPar.Size,...
                                 'User',InPar.User,'Pass',InPar.Pass,...
                                 InPar.QueryPar{:});
                             
Prop  = VO.ZTF.irsa_table2prop(Cat,'ImType',InPar.ImType,'Product',InPar.Product);      

Links = VO.ZTF.irsa_image_link(Prop,InPar.constructPar{:});

if (InPar.GetFiles)
    switch lower(InPar.GetN)
        case 'all'
            LinksR = Links;
        case 'first'
            LinksR = Links{1};
        case 'last'
            LinksR = Links{end};
        otherwise
            error('Unknown GetN option');
    end
    
    Files = www.pwget(LinksR,InPar.pwgetExtra,InPar.MaxGet,[],'wget');
else
    Files = {};
end
