function prepareSNdata(snname,redshift,RAh,decd)
% Prepare Observation data and ancillary data for a SN
% Package: AstroUtil.supernove.SOPRANOS
% Description: Prepare Observation data and ancillary data for a SN
% Input  : - The name of the supernova (The observations are expected to be
%            in a CSV file, txt file, or out_PTF48R file with the SN name).
%          - the SN redshift
%          - the SN RA in sexagesimal hours
%          - the SN dec in sexagesimal degrees
% Output : - A file named <sn_name>_data.mat
%               
% See also: AstroUtil.supernova.SOPRANOS.calcGrid
%           AstroUtil.supernova.SOPRANOS.prepare_LCs
% Tested : Matlab 9.5
%     By : Noam Ganot                      Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstroUtil.supernova.SOPRANOS.prepareSNdata('PTF12gnt')
% Reliable: 2
%--------------------------------------------------------------------------

if strcmpi(snname(1:3),'ztf')
    txtfile = sprintf('%s.txt',snname);
    csvfile = sprintf('%s.csv',snname);
    if exist(txtfile,'file')
        Table = AstroUtil.supernova.SOPRANOS.readZTFtxt(txtfile);
        bands = AstroUtil.supernova.SOPRANOS.ztfBandsTxt(Table);
    elseif exist(csvfile,'file')
        ztfTable = AstroUtil.supernova.SOPRANOS.readZTF(csvfile);
        bands = AstroUtil.supernova.SOPRANOS.ztfBands(ztfTable);
    end
    RArad  = celestial.coo.convertdms(RAh,'SH','r');
    decRad = celestial.coo.convertdms(decd,'SD','r');
elseif strcmpi(snname(1:3),'ptf')
    [Table, RA, dec, redshift] = AstroUtil.supernova.SOPRANOS.readPTFout(snname);
    bands = ztfBandsTxt(Table);
    
    RArad = RA./180*pi;
    decRad =dec./180*pi;
    
    RAh = celestial.coo.convertdms(RA,'d','SH');
    decd = celestial.coo.convertdms(dec,'d','SD');
end

save(sprintf('%s_data.mat',snname),'snname','redshift','bands', 'RAh', 'decd', 'RArad','decRad','-v7.3');