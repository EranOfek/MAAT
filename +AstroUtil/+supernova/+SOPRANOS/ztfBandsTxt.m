function bands=ztfBandsTxt(Tab)
% Split the table into different bands and read the filter Obj
% Package: AstroUtil.supernove.SOPRANOS
% Description: Split the table into different bands and read the filter Obj
% Input  : - The table read by readZTFtxt
% Output : - A file named <sn_name>_<model>_LCs.mat
%               
% See also: AstroUtil.supernova.SOPRANOS.readZTFtxt
% Tested : Matlab 9.5
%     By : Noam Ganot                      Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstroUtil.supernova.SOPRANOS.ztfBands.Txt(Tab)
% Reliable: 2
%--------------------------------------------------------------------------

filters = unique(Tab.filter);
instruments = unique(Tab.instr);

iband = 0;
for ifilter=1:length(filters)
    filter = filters(ifilter);
    for iinstr=1:length(instruments)
        instr=instruments(iinstr);
        bandTab = Tab(Tab.filter == filter & Tab.instr == instr,:);
        if ~isempty(bandTab)
            iband = iband+1;
            bands{iband}.filtername = filter;
            bands{iband}.instrumentname = instr;
            bands{iband}.raw    = bandTab;
            vars = bandTab.Properties.VariableNames;
            bands{iband}.LC     = bandTab(bandTab.mag~=99&bandTab.magerr~=99,bandTab.Properties.VariableNames);
            bands{iband}.LC.filter = []; bands{iband}.LC.instr = [];
            bands{iband}.limits = bandTab(bandTab.mag==99,bandTab.Properties.VariableNames);
            bands{iband}.limits.filter = []; bands{iband}.limits.instr = [];
            switch instr
                case "P48+ZTF"
                    switch filter
                        case {'r_p48','r'}
                            bands{iband}.filterObj=AstFilter.get('ZTF','r');
                        case {'g_p48', 'g'}
                            bands{iband}.filterObj=AstFilter.get('ZTF','g');
                        case {'i_p48','i'}
                            bands{iband}.filterObj=AstFilter.get('ZTF','i');
                        otherwise
                            fprintf('unknown filter %s for P48',filter);
                    end
                case "P60+SEDM"
                      switch lower(filter)
                        case {'r_sdss','r'}
                            bands{iband}.filterObj=AstFilter.get('SDSS','r');
                        case {'g_sdss','g'}
                            bands{iband}.filterObj=AstFilter.get('SDSS','g');
                        case {'i_sdss','i'}
                            bands{iband}.filterObj=AstFilter.get('SDSS','i');
                        case {'u_sdss','u'}
                            bands{iband}.filterObj=AstFilter.get('SDSS','u');
                        otherwise
                            fprintf('unknown filter %s for P60',filter);
                      end
                case {"SWIFT+UVOT","Swift+UVOT"}
                    switch filter
                        case 'UVW1'
                            bands{iband}.filterObj=AstFilter.get('Swift/UVOT','UVW1');
                        case 'UVW2'
                            bands{iband}.filterObj=AstFilter.get('Swift/UVOT','UVW2');
                        case 'UVM2'
                            bands{iband}.filterObj=AstFilter.get('Swift/UVOT','UVM2');
                        case {'u_swift','U'}
                            bands{iband}.filterObj=AstFilter.get('Swift/UVOT','U');
                        case {'b_swift','B'}
                            bands{iband}.filterObj=AstFilter.get('Swift/UVOT','B');
                        case {'v_swift','V'}
                            bands{iband}.filterObj=AstFilter.get('Swift/UVOT','V');
                        otherwise
                            fprintf('unknown filter %s for Swift/UVOT',filter);
                    end
                case "GALEX"
                    switch filter
                        case 'NUV'
                            bands{iband}.filterObj=AstFilter.get('GALEX','NUV');
                        case 'FUV'
                            bands{iband}.filterObj=AstFilter.get('GALEX','FUV');
                        otherwise
                            fprintf('unknown filter %s for GALEX',filter);
                    end
                 case "P48+PTF"
                    switch filter
                        case 'r_p48'
                            bands{iband}.filterObj=AstFilter.get('PTF','R');
                        case 'g_p48'
                            bands{iband}.filterObj=AstFilter.get('PTF','g');
                        otherwise
                            fprintf('unknown filter %s for PTF',filter);
                    end
                       
                      
            end
            if isempty(bands{iband}.filterObj)
                error('Cannot find %s %s in AstFilterCat.mat\n',instr,filter);
            end
            if ~ismember('MJD', bands{iband}.LC.Properties.VariableNames) 
                bands{iband}.LC.MJD  = convert.time(bands{iband}.LC.jd,'jd','MJD');
                bands{iband}.limits.MJD  = convert.time(bands{iband}.limits.jd,'jd','MJD');
            end
            bands{iband}.LC  = sortrows(unique(bands{iband}.LC),'MJD');
            bands{iband}.limits  = sortrows(unique(bands{iband}.limits),'MJD');
        end
    end
end

