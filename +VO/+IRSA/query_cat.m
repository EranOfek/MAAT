function Cat=query_cat(CatName,RA,Dec,Radius,varargin)
% Query IPAC/IRSA catalog.
% Package: VO.IRSA
% Description: Query IPAC/IRSA catalog.
% Input  : - Catalog name (e.g., 'wise_allwise_p3as_psd').
%          - Object name (e.g., 'M81') or coordinates in [H M S],
%            sexagesimal string or radians.
%          - Object declination [Sign D M S], sexagesimal or radians.
%            If empty then assume that the RA contains the object name.
%          - Search radius [arcsec]. Default is 60.
%          * Arbitrary number of pairs of arguments: ...,keyword,value,...
%            where keyword are one of the followings:
%            'Type'   - Search type 'circ'|'cone'|'box'. Default is 'cone'.
%            'Columns'- Cell array of column names to query.
%                        If empty retrieve all columns. Default is empty.
% Output : - 
% Tested : Matlab R2011b
%     By : Eran O. Ofek                    Sep 2016
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Reference: http://irsa.ipac.caltech.edu/applications/Gator/GatorAid/irsa/catsearch.html
% Example: Cat=VO.IRSA.query_cat('wise_allwise_p3as_psd','M81');
% Reliable: 2
%--------------------------------------------------------------------------

Def.Dec    = [];
Def.Radius = 60;
if (nargin==2)
    Dec    = Def.Dec;
    Radius = Def.Radius;
elseif (nargin==3)
    Radius = Def.Radius;
else
    % do nothing
end


DefV.Type               = 'cone';
DefV.Units              = 'arcsec';
DefV.Columns            = [];
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


if (isempty(Dec))
    % assume RA contain object name
    ObjStr = RA;
else
    % convert to RA,Dec string
    RA  = celestial.coo.convertdms(RA,'gH','SH');
    Dec = celestial.coo.convertdms(Dec,'gD','SD');
    RA  = regexprep(RA,':','h+');
    RA(7)  = 'm';
    RA(15) = 's';
    Dec  = regexprep(Dec,':','d+');
    Dec(8)  = 'm';
    Dec(15) = 's';
    ObjStr  = sprintf('%s%s',RA,Dec);
end

%curl -o search_results.tbl "http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?catalog=wise_allwise_p3as_psd&spatial=cone&radius=300&radunits=arcsec&objstr=00h+42m+44.32s+41d+16m+08.5s&size=300&outfmt=1&selcols=w1mpro,w2mpro,w3mpro,w4mpro"

if (strcmp(InPar.Type,'circ'))
    InPar.Type = 'cone';
end

Url = VO.IRSA.irsa_db_url;
%Type = 'cone';   % 'cone'|'box'|...
%Radius = 300;
%ObjStr = '00h+42m+44.32s+41d+16m+08.5s';
Size   = Radius;  % required if 'box'
OutFormat = 1;  % 0: HTML (default)
                % 1: ASCII table 
                % 2: SVC (software handshaking structure) message
                % 3: VO Table in XML
                % 6: program interface in XML
%Columns = {'w1mpro','w2mpro','w3mpro','w4mpro'};               

if (isempty(InPar.Columns))
    % retrieve all column names
    Cols = VO.IRSA.wget_cat_columns(CatName);
    InPar.Columns = Cols.Name;
end


Ncol    = numel(InPar.Columns);

AllCol = '';
for Icol=1:1:Ncol
    if (Icol==1)
        Delim = '';
    else
        Delim = ',';
    end
    AllCol = sprintf('%s%s%s',AllCol,Delim,InPar.Columns{Icol});
end

                
Tbl = urlread(sprintf('%snph-query?catalog=%s&spatial=%s&radius=%f&radunits=%s&objstr=%s&size=%f&outfmt=1&selcols=%s',...
                        Url,CatName,InPar.Type,Radius,InPar.Units,ObjStr,Size,AllCol));
                    % 'post',{'Username','username';'Password','password'} );
            
Cat  = VO.IRSA.read_ipac_table(Tbl,'str');
                
                    
                    