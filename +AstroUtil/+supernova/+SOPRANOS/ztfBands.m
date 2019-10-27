function bands=ztfBands(ztfTab)
% Split the table into different bands and read the filter Obj
% Package: AstroUtil.supernove.SOPRANOS
% Description: Split the table into different bands and read the filter Obj
% Input  : - The table read by readZTF
% Output : - A file named <sn_name>_<model>_LCs.mat
%               
% See also: AstroUtil.supernova.SOPRANOS.readZTF
% Tested : Matlab 9.5
%     By : Noam Ganot                      Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstroUtil.supernova.SOPRANOS.ztfBands(Tab)
% Reliable: 2
%--------------------------------------------------------------------------

filters = unique(ztfTab.filter);
instruments = unique(ztfTab.instrument);

iband = 0;
for ifilter=1:length(filters)
    filter = filters(ifilter);
    for iinstrument=1:length(instruments)
        instrument=instruments(iinstrument);
        bandTab = ztfTab(ztfTab.filter == filter & ztfTab.instrument == instrument,:);
        if ~isempty(bandTab)
            iband = iband+1;
            bands{iband}.filtername = filter;
            bands{iband}.instrumentname = instrument;
            bands{iband}.limits = bandTab(bandTab.magpsf==99,bandTab.Properties.VariableNames);
            bands{iband}.limits.filter = []; bands{iband}.limits.instrument = [];
            bands{iband}.LC    = bandTab(bandTab.magpsf~=99&bandTab.sigmamagpsf~=99,bandTab.Properties.VariableNames);
            bands{iband}.LC.filter = []; bands{iband}.LC.instrument = [];
            switch instrument
                case "P48+ZTF"
                    switch filter
                        case 'r'
                            bands{iband}.filterObj=AstFilter.get('PTF','r');
                        case 'g'
                            bands{iband}.filterObj=AstFilter.get('PTF','g');
                        otherwise
                            fprintf('unknown filter %s for P48',filter);
                    end
                case "P60+SEDM"
                      switch lower(filter)
                        case 'r'
                            bands{iband}.filterObj=AstFilter.get('SDSS','r');
                        case 'g'
                            bands{iband}.filterObj=AstFilter.get('SDSS','g');
                        case 'i'
                            bands{iband}.filterObj=AstFilter.get('SDSS','i');
                        case 'u'
                            bands{iband}.filterObj=AstFilter.get('SDSS','u');
                        otherwise
                            fprintf('unknown filter %s for P60',filter);
                      end
                case "Swift+UVOT"
                    switch filter
                        case 'UVW1'
                            bands{iband}.filterObj=AstFilter.get('Swift/UVOT','UVW1');
                        case 'UVW2'
                            bands{iband}.filterObj=AstFilter.get('Swift/UVOT','UVW2');
                        case 'UVM2'
                            bands{iband}.filterObj=AstFilter.get('Swift/UVOT','UVM2');
                        case 'u'
                            bands{iband}.filterObj=AstFilter.get('Swift/UVOT','U');
                        case 'B'
                            bands{iband}.filterObj=AstFilter.get('Swift/UVOT','B');
                        case 'V'
                            bands{iband}.filterObj=AstFilter.get('Swift/UVOT','V');
                        otherwise
                            fprintf('unknown filter %s for Swift/UVOT',filter);
                    end
            end
            if isempty(bands{iband}.filterObj)
                error('Cannot find %s %s in AstFilterCat.mat\n',instrument,filter);
            end
            FluxPnt = convert.flux(bands{iband}.LC.magpsf,'AB','cgs/A',bands{iband}.filterObj.pivot_wl_photon,'A');
            FluxErr = bands{iband}.LC.sigmamagpsf.*FluxPnt*log(10)/2.5;
            bands{iband}.LC.MJD  = convert.time(bands{iband}.LC.jdobs,'JD','MJD');
            bands{iband}.LC.flux = FluxPnt;
            bands{iband}.LC.fluxerr = FluxErr;
            bands{iband}.LC.ErrM = FluxErr;
            bands{iband}.LC.ErrP = FluxErr;
            LimitPnt = convert.flux(bands{iband}.limits.limmag,'AB', 'cgs/A',bands{iband}.filterObj.pivot_wl_photon,'A');
            bands{iband}.limits.MJD  = convert.time(bands{iband}.limits.jdobs,'JD','MJD');
            bands{iband}.limits.flux = LimitPnt;
        end
    end
end