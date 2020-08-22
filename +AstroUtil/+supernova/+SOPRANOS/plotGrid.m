function plotGrid(filename,mintpoints,peakNumber,filter,suffix)
% Plot the likelihood grid calculated by calcGrid
% Package: AstroUtil.supernove.SOPRANOS
% Description: Plot the likelihood grid calculated by calcGrid.
% Input  : - grid file name
%          - minimal transient points constaint for each band (0 means no
%                 constraint).
%          - peakNumber for plotting the marginal distribution in case of
%                 multiple peak solution (default is 1)
%          - External filter to plot subgrid or to mainpulate the valid
%            region. The filter is a string contains matlab code which is
%            executed using eval to change the Vectors which span the grid
%            or the valid variable which defines the valid grid points.
%          - suffix for the fMax results file which describes the filter
% Output : - Figures and text plotted to the command line.
%               
% Tested : Matlab 9.5
%     By : Noam Ganot                      Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example: AstroUtil.supernova.SOPRANOS.plotGrid('PTF12ffs_msw_grid.mat', [0 1], 1, 'VecEbv=VecEbv(VecEbv<=0.09);','Ebv.lte.0.09');
% Reliable: 2
%--------------------------------------------------------------------------

if (nargin==2) || (nargin==4 && isempty(peakNumber))
    peakNumber = 1;
end

gridfile = matfile(filename);
nparams = length(gridfile.BGbands)+6;
VecRs   = gridfile.VecRs;
VecVs   = gridfile.VecVs;
Vect0   = gridfile.Vect0;
VecMs   = gridfile.VecMs;
VecFrho = gridfile.VecFrho;
VecEbv  = gridfile.VecEbv;
valid = true;

if (nargin>=4)
    eval(filter)
    is_Rs = find(gridfile.VecRs==VecRs(1));
    ie_Rs = find(gridfile.VecRs==VecRs(end));
    is_Vs = find(gridfile.VecVs==VecVs(1));
    ie_Vs = find(gridfile.VecVs==VecVs(end));
    is_t0 = find(gridfile.Vect0==Vect0(1));
    ie_t0 = find(gridfile.Vect0==Vect0(end));
    is_Ms = find(gridfile.VecMs==VecMs(1));
    ie_Ms = find(gridfile.VecMs==VecMs(end));
    is_Fr = find(gridfile.VecFrho==VecFrho(1));
    ie_Fr = find(gridfile.VecFrho==VecFrho(end));
    is_Eb = find(gridfile.VecEbv==VecEbv(1));
    ie_Eb = find(gridfile.VecEbv==VecEbv(end));
    if length(valid)>1
        valid=valid(is_t0:ie_t0,is_Rs:ie_Rs,is_Vs:ie_Vs,is_Ms:ie_Ms,is_Fr:ie_Fr,is_Eb:ie_Eb);
    end
else
    is_Rs = 1; ie_Rs = length(VecRs);
    is_Vs = 1; ie_Vs = length(VecVs);
    is_t0 = 1; ie_t0 = length(Vect0);
    is_Ms = 1; ie_Ms = length(VecMs);
    is_Fr = 1; ie_Fr = length(VecFrho);
    is_Eb = 1; ie_Eb = length(VecEbv);
end

PDFmap = zeros(length(Vect0),length(VecRs),length(VecVs),length(VecMs),length(VecFrho),length(VecEbv));
valid = valid&true(size(PDFmap));
for it=is_t0:ie_t0
    chi2 = double(gridfile.chi2(it,is_Rs:ie_Rs,is_Vs:ie_Vs,is_Ms:ie_Ms,is_Fr:ie_Fr,is_Eb:ie_Eb));
    dof  = double(gridfile.points(it,is_Rs:ie_Rs,is_Vs:ie_Vs,is_Ms:ie_Ms,is_Fr:ie_Fr,is_Eb:ie_Eb))-nparams.*ones(size(chi2));
    PDFmap(it-is_t0+1,:,:,:,:,:)= chi2pdf(chi2,dof);
    valid(it-is_t0+1,:,:,:,:,:)=valid(it-is_t0+1,:,:,:,:,:)&(dof>0);
end
clear chi2 dof

if nargin==1
    mintpoints = zeros(size(gridfile.NTransient));
end
NTrans = cell2mat(gridfile.NTransient(1,1));
valid_minpoints = NTrans(is_t0:ie_t0,is_Rs:ie_Rs,is_Vs:ie_Vs,is_Ms:ie_Ms,is_Fr:ie_Fr)>=mintpoints(1);
for iband = 2:length(gridfile.NTransient)
    NTrans = cell2mat(gridfile.NTransient(1,iband));
    valid_minpoints=valid_minpoints&NTrans(is_t0:ie_t0,is_Rs:ie_Rs,is_Vs:ie_Vs,is_Ms:ie_Ms,is_Fr:ie_Fr)>=mintpoints(iband);
end
valid = valid & valid_minpoints;
clear NTrans valid_minpoints;

PDFmap(~valid)=0;

if length(gridfile.VecFrho)>1
    PDF_Rv = squeeze(trapz(Vect0,trapz(VecMs,trapz(log10(VecFrho),trapz(VecEbv,PDFmap,6),5),4),1));
else
    PDF_Rv = squeeze(trapz(Vect0,trapz(VecMs,trapz(VecEbv,PDFmap,6),4),1));
end

valid_Rv = squeeze(any(any(any(any(valid,6),5),4),1));

[hPDF,Vs,Rs,iVs,iRs]=plotPDF(PDF_Rv, valid_Rv, VecVs, VecRs);
[hCDF,~,~,SigmaVs,SigmaRs,VsPDF,RsPDF] = plotCDF(PDF_Rv, valid_Rv, VecVs, VecRs);
marginalizedPDF=squeeze(PDFmap(:,iRs,iVs,:,:,:));
[val, imax] = max(marginalizedPDF(:));
if length(gridfile.VecFrho)>1
    [it0,iMs,ifrho,iEbv] = ind2sub(size(marginalizedPDF),imax);
else
    [it0,iMs,iEbv] = ind2sub(size(marginalizedPDF),imax);
    ifrho = 1;
end

t0 = Vect0(it0);
Ms = VecMs(iMs);
frho = VecFrho(ifrho);
Ebv = VecEbv(iEbv);

% when filter is applied the grid plotted may be a subset of the one in
% the grid file. Find the orignial values in order to get the chi2 and
% dof values of the peak out of the grid file.
it0g   = find(gridfile.Vect0   == t0);
iRsg   = find(gridfile.VecRs   == Rs);
iVsg   = find(gridfile.VecVs   == Vs);
iMsg   = find(gridfile.VecMs   == Ms);
ifrhog = find(gridfile.VecFrho == frho);
iEbvg  = find(gridfile.VecEbv  == Ebv);
chi2_val = gridfile.chi2(it0g,iRsg,iVsg,iMsg,ifrhog,iEbvg);
dof_val = gridfile.points(it0g,iRsg,iVsg,iMsg,ifrhog,iEbvg)-nparams;


% find all secondary peaks higher the 1% then the absoulte peak.
PDF_Rv(PDF_Rv<0.01*max(PDF_Rv(:)))=0;

[iRs_peaks,iVs_peaks]=ind2sub(size(PDF_Rv),find(imregionalmax(PDF_Rv,8)));
Rs_peaks = VecRs(iRs_peaks,1);
Vs_peaks = VecVs(iVs_peaks,1);
for ipeak = 1:length(Rs_peaks)
    marginalizedPDF=squeeze(PDFmap(:,iRs_peaks(ipeak),iVs_peaks(ipeak),:,:,:));
    [~, imax] = max(marginalizedPDF(:));
    if length(gridfile.VecFrho)>1
        [it0,iMs,ifrho,iEbv] = ind2sub(size(marginalizedPDF),imax);
    else
        [it0,iMs,iEbv] = ind2sub(size(marginalizedPDF),imax);
        ifrho = 1;
    end
    
    t0_peaks(ipeak) = Vect0(it0);
    Ms_peaks(ipeak) = VecMs(iMs);
    frho_peaks(ipeak) = VecFrho(ifrho);
    Ebv_peaks(ipeak) = VecEbv(iEbv);

    % when filter is applied the grid plotted may be a subset of the one in
    % the grid file. Find the orignial values in order to get the chi2 and
    % dof values of the peak out of the grid file.
    it0g   = find(gridfile.Vect0   == Vect0(it0));
    iRsg   = find(gridfile.VecRs   == Rs_peaks(ipeak));
    iVsg   = find(gridfile.VecVs   == Vs_peaks(ipeak));
    iMsg   = find(gridfile.VecMs   == VecMs(iMs));
    ifrhog = find(gridfile.VecFrho == VecFrho(ifrho));
    iEbvg  = find(gridfile.VecEbv  == VecEbv(iEbv));
    chi2_peaks(ipeak) = gridfile.chi2(it0g,iRsg,iVsg,iMsg,ifrhog,iEbvg);
    dof_peaks(ipeak) = gridfile.points(it0g,iRsg,iVsg,iMsg,ifrhog,iEbvg)-nparams;
end

marginalizedPDF = squeeze(trapz(VecVs/10^8.5,trapz(VecRs,PDFmap,2),3));
[ht0, t0PDF, Sigmat0] = plotMargProbt0(marginalizedPDF, VecMs, VecEbv, log10(VecFrho), Vect0, [], t0);
[hMs, MsPDF, SigmaMs] = plotMargProbMs(marginalizedPDF, VecMs, VecEbv, log10(VecFrho), Vect0, [],Ms);
[hEbv, EbvPDF, SigmaEbv] = plotMargProbEbv(marginalizedPDF, VecMs, VecEbv, log10(VecFrho), Vect0, [], Ebv);
if length(gridfile.VecFrho)>1
    [hfrho, frhoPDF, SigmaFrho] = plotMargProbFrho(marginalizedPDF, VecMs, VecEbv, VecFrho, Vect0, [], frho);
else
    SigmaFrho = [0 0];
end

for ipeak = 1:length(iRs_peaks)
    [RssigmaM(ipeak),~,RssigmaP(ipeak)]=oneSigmaMove(VecRs,RsPDF,Rs_peaks(ipeak));
    RssigmaM(ipeak) = Rs_peaks(ipeak) - RssigmaM(ipeak);
    RssigmaP(ipeak) = RssigmaP(ipeak) - Rs_peaks(ipeak);

    [VssigmaM(ipeak),~,VssigmaP(ipeak)]=oneSigmaMove(VecVs/10^8.5,VsPDF,Vs_peaks(ipeak)/10^8.5);
    VssigmaM(ipeak) = Vs_peaks(ipeak)/10^8.5 - VssigmaM(ipeak);
    VssigmaP(ipeak) = VssigmaP(ipeak) - Vs_peaks(ipeak)/10^8.5;
    
    [t0sigmaM(ipeak),~,t0sigmaP(ipeak)]=oneSigmaMove(Vect0,t0PDF,t0_peaks(ipeak));
    t0sigmaM(ipeak) = t0_peaks(ipeak) - t0sigmaM(ipeak);
    t0sigmaP(ipeak) = t0sigmaP(ipeak) - t0_peaks(ipeak);

    [MssigmaM(ipeak),~,MssigmaP(ipeak)]=oneSigmaMove(VecMs,MsPDF,Ms_peaks(ipeak));
    MssigmaM(ipeak) = Ms_peaks(ipeak) - MssigmaM(ipeak);
    MssigmaP(ipeak) = MssigmaP(ipeak) - Ms_peaks(ipeak);

    [EbvsigmaM(ipeak),~,EbvsigmaP(ipeak)]=oneSigmaMove(VecEbv,EbvPDF,Ebv_peaks(ipeak));
    EbvsigmaM(ipeak) = Ebv_peaks(ipeak) - EbvsigmaM(ipeak);
    EbvsigmaP(ipeak) = EbvsigmaP(ipeak) - Ebv_peaks(ipeak);

    [frhosigmaM(ipeak),~,frhosigmaP(ipeak)]=oneSigmaMove(log10(VecFrho),frhoPDF,log10(frho_peaks(ipeak)));    
    frhosigmaM(ipeak) = frho_peaks(ipeak)- 10^frhosigmaM(ipeak);
    frhosigmaP(ipeak) = 10^frhosigmaP(ipeak) - frho_peaks(ipeak);
    
    P(ipeak) = chi2pdf(double(chi2_peaks(ipeak)), double(dof_peaks(ipeak)));
end

[P,peakInd] = sort(P,'descend');

hPDFpeaks =plot(hPDF.Children(3),Vs_peaks/10^8.5,Rs_peaks,'+');
figure(hCDF); subplot(4,4,[1:3 5:7 9:11]); hold on; hCDFpeaks = plot(Vs_peaks/10^8.5,Rs_peaks,'k+');



fprintf('Maximal proablility %4.3f, %4.2f/%d (chi2/dof)\nRs=%6.2f_{-%5.2f}^{+%5.2f}\nv_{s*,8.5}=%6.4f_{-%5.4f}^{+%5.4f}\nt_0=%6.2f_{-%4.2f}^{+%4.2f}\nMs=%3.1f_{-%3.1f}^{+%3.1f}\nf_rho=%6.5f_{-%6.5f}^{%6.5f}\nEbv=%7.5f_{-%7.5f}^{%7.5f}\n',...
    val, chi2_val, dof_val, ...
    Rs, Rs-SigmaRs(1), SigmaRs(2)-Rs, ...
    Vs/10^8.5, (Vs)/10^8.5-SigmaVs(1), SigmaVs(2)-(Vs)/10^8.5, ...
    t0, t0-Sigmat0(1), Sigmat0(2)- t0, ...
    Ms, Ms-SigmaMs(1), SigmaMs(2)-Ms, ...
    frho, frho-SigmaFrho(1), SigmaFrho(2)-frho, ...
    Ebv, Ebv-SigmaEbv(1), SigmaEbv(2)-Ebv          );

nPeaks = length(Rs_peaks);
fprintf('\n');
line = 'Parameter          ';
for iPeak=1:nPeaks,line = [line sprintf('& Peak \\\\#%d                    ',iPeak)];end
line = [line '\\\\\n'];
fprintf(line);
fprintf('\\hline\n');

line='$R_{*}$ [$R_\\sun$] ';
for iPeak=1:nPeaks 
    line = [line sprintf('& $%4.0f_{-%4.0f}^{+%4.0f}$      ', Rs_peaks(peakInd(iPeak)), RssigmaM(peakInd(iPeak)), RssigmaP(peakInd(iPeak)))];
end
line = [line '\\\\\n'];
fprintf(line);

line ='$v_{s*,8.5}$       ';
for iPeak=1:nPeaks
    line = [line sprintf('& $%4.2f_{-%4.2f}^{+%4.2f}$      ', Vs_peaks(peakInd(iPeak))/10^8.5, VssigmaM(peakInd(iPeak)), VssigmaP(peakInd(iPeak)))];
end
line = [line '\\\\\n'];
fprintf(line);

line ='$M_{ej}$ [$M_\\sun$]';
for iPeak=1:nPeaks
    line = [line sprintf('& $%4.1f_{-%4.1f}^{+%4.1f}$      ', Ms_peaks(peakInd(iPeak)), MssigmaM(peakInd(iPeak)), MssigmaP(peakInd(iPeak)))];
end
line = [line '\\\\\n'];
fprintf(line);

line ='$t_{\\rm ref}$ [MJD]';
for iPeak=1:nPeaks
    line = [line sprintf('& $%6.2f_{-%4.2f}^{+%4.2f}$  ', t0_peaks(peakInd(iPeak)), t0sigmaM(peakInd(iPeak)), t0sigmaP(peakInd(iPeak)))];
end
line = [line '\\\\\n'];
fprintf(line);

line ='$f_{\\rho}$         ';
for iPeak=1:nPeaks
    line = [line sprintf('& $%5.3f_{-%5.3f}^{+%5.3f}$   ', frho_peaks(peakInd(iPeak)), frhosigmaM(peakInd(iPeak)), frhosigmaP(peakInd(iPeak)))];
end
line = [line '\\\\\n'];
fprintf(line);

line ='$E_{B-V}$          ';
for iPeak=1:nPeaks
    line = [line sprintf('& $%5.3f_{-%5.3f}^{+%5.3f}$   ', Ebv_peaks(peakInd(iPeak)), EbvsigmaM(peakInd(iPeak)), EbvsigmaP(peakInd(iPeak)))];
end
line = [line '\\\\\n'];
fprintf(line);

line ='$\\chi^2$/dof       ';
for iPeak=1:nPeaks
    line = [line sprintf('& %6.2f/%3u                  ',chi2_peaks(peakInd(iPeak)), dof_peaks(peakInd(iPeak)))];
end
line = [line '\\\\\n'];
fprintf(line);

fprintf('\\hline\n\n')

save grid_max t0 Rs Vs Ms frho Ebv Rs_peaks Vs_peaks Ms_peaks frho_peaks t0_peaks Ebv_peaks P
bands = gridfile.bands;
for iBand = 1:length(bands)
    BG = cell2mat(gridfile.bg(1,iBand));
    if length(BG)>1
        fprintf('%s %s background: %6.4f\n',bands{iBand}.instrumentname, bands{iBand}.filtername, BG(it0, iRs, iVs, iMs, iFrho, iEbv));
        bg(iBand) = BG(it0, iRs, iVs, iMs, iFrho, iEbv);
    else
        fprintf('%s %s background: %6.4f\n',bands{iBand}.instrumentname, bands{iBand}.filtername, BG);
        bg(iBand) = BG;
    end
end
clear BG

if nargin<4
    filter_suffix = '';
elseif nargin==4
    filter_suffix = '_filter';
else
    filter_suffix = sprintf('_%s',suffix);
end

if strcmp(filename(end-3:end),'.mat')
    filename = filename(1:end-4);
end

if all(mintpoints==0)
    results_fname=sprintf('fMax_results_%s%s.mat',filename,filter_suffix);
    figs_name = sprintf('%s%s',filename,filter_suffix);
else
    results_fname=sprintf('fMax_results_%s_%s%s.mat',filename,sprintf('%d',mintpoints),filter_suffix);
    figs_name = sprintf('%s_%s%s',filename,sprintf('%d',mintpoints),filter_suffix);
end

if ~exist(results_fname,'file')
    redshift = gridfile.redshift;
    ProgType = gridfile.ProgType;
    bands   = gridfile.bands; 
    model = gridfile.model;
    transientStart = gridfile.transientStart;
    transientEnd   = gridfile.transientEnd;
    save fmax_input filename Vs_peaks Rs_peaks bg results_fname mintpoints VecRs VecVs Vect0 VecMs VecFrho VecEbv redshift ProgType bands model transientStart transientEnd
    if isunix
        !screen -dmS findMaximum nice -n 19 matlab -nodisplay -nosplash -nodesktop -r "load('fmax_input.mat','Rs_peaks','Vs_peaks','bg','results_fname', 'mintpoints');[Vs, Rs, Ms, Ebv, Rv,  Vect0, PDFt0, t0, chi2values] = AstroUtil.supernova.SOPRANOS.findMaximum('fmax_input.mat',Vs_peaks,Rs_peaks,bg,mintpoints);save(results_fname, 'Vs', 'Rs', 'Ms', 'Ebv', 'Rv', 'Vect0', 'PDFt0', 't0', 'chi2values');exit" -logfile findMaximum.log &
    elseif ispc
        !start /BELOWNORMAL matlab -nosplash -nodesktop -minimize -r "load('fmax_input.mat','Rs_peaks','Vs_peaks','bg','results_fname', 'mintpoints');[Vs, Rs, Ms, Ebv, Rv,  Vect0, PDFt0, t0, chi2values] = AstroUtil.supernova.SOPRANOS.findMaximum('fmax_input.mat',Vs_peaks,Rs_peaks,bg,mintpoints);save(results_fname, 'Vs', 'Rs', 'Ms', 'Ebv', 'Rv', 'Vect0', 'PDFt0', 't0', 'chi2values');exit" -logfile "findMaximum.log"
    end
    return
else
    load(results_fname);
end

% find confidence level for each peak and find their probablity order
nPeaks = length(chi2values);
P = zeros(nPeaks,1);
for iPeak = 1:nPeaks
    [chi2values(iPeak).all.RssigmaM,~,chi2values(iPeak).all.RssigmaP]=oneSigmaMove(VecRs,RsPDF,chi2values(iPeak).all.Rs);
    chi2values(iPeak).all.RssigmaM = chi2values(iPeak).all.Rs - chi2values(iPeak).all.RssigmaM;
    chi2values(iPeak).all.RssigmaP = chi2values(iPeak).all.RssigmaP - chi2values(iPeak).all.Rs;

    [chi2values(iPeak).all.VssigmaM,~,chi2values(iPeak).all.VssigmaP]=oneSigmaMove(VecVs/10^8.5,VsPDF,chi2values(iPeak).all.Vs/10^8.5);
    chi2values(iPeak).all.VssigmaM = chi2values(iPeak).all.Vs/10^8.5 - chi2values(iPeak).all.VssigmaM;
    chi2values(iPeak).all.VssigmaP = chi2values(iPeak).all.VssigmaP - chi2values(iPeak).all.Vs/10^8.5;
    
    [chi2values(iPeak).all.t0sigmaM,~,chi2values(iPeak).all.t0sigmaP]=oneSigmaMove(Vect0,t0PDF,chi2values(iPeak).all.t0);
    chi2values(iPeak).all.t0sigmaM = chi2values(iPeak).all.t0 - chi2values(iPeak).all.t0sigmaM;
    chi2values(iPeak).all.t0sigmaP = chi2values(iPeak).all.t0sigmaP - chi2values(iPeak).all.t0;

    [chi2values(iPeak).all.MssigmaM,~,chi2values(iPeak).all.MssigmaP]=oneSigmaMove(VecMs,MsPDF,chi2values(iPeak).all.Ms);
    chi2values(iPeak).all.MssigmaM = chi2values(iPeak).all.Ms - chi2values(iPeak).all.MssigmaM;
    chi2values(iPeak).all.MssigmaP = chi2values(iPeak).all.MssigmaP - chi2values(iPeak).all.Ms;

    [chi2values(iPeak).all.EbvsigmaM,~,chi2values(iPeak).all.EbvsigmaP]=oneSigmaMove(VecEbv,EbvPDF,chi2values(iPeak).all.Ebv);
    chi2values(iPeak).all.EbvsigmaM = chi2values(iPeak).all.Ebv - chi2values(iPeak).all.EbvsigmaM;
    chi2values(iPeak).all.EbvsigmaP = chi2values(iPeak).all.EbvsigmaP - chi2values(iPeak).all.Ebv;

    [chi2values(iPeak).all.frhosigmaM,~,chi2values(iPeak).all.frhosigmaP]=oneSigmaMove(log10(VecFrho),frhoPDF,log10(chi2values(iPeak).all.frho));    
    chi2values(iPeak).all.frhosigmaM = chi2values(iPeak).all.frho - 10^chi2values(iPeak).all.frhosigmaM;
    chi2values(iPeak).all.frhosigmaP = 10^chi2values(iPeak).all.frhosigmaP - chi2values(iPeak).all.frho;
    
    P(iPeak) = chi2pdf(chi2values(iPeak).all.chi2, chi2values(iPeak).all.dof);
end
[~,peakInd] = sort(P,'descend');

delete(hPDFpeaks);
delete(hCDFpeaks);
clear Rs_peaks Vs_peaks
for iPeak = 1:nPeaks
    Rs_peaks(iPeak) = chi2values(iPeak).all.Rs;
    Vs_peaks(iPeak) = chi2values(iPeak).all.Vs;
end
figure(hPDF);hPDFpeaks =plot(hPDF.Children(3),Vs_peaks/10^8.5,Rs_peaks,'+');
savefig(hPDF,sprintf('%s_PDF.fig',figs_name));
figure(hCDF); subplot(4,4,[1:3 5:7 9:11]); hold on; hCDFpeaks = plot(Vs_peaks/10^8.5,Rs_peaks,'k+');
savefig(hCDF,sprintf('%s_CDF.fig',figs_name));

ht0 = plotMargProbt0(marginalizedPDF, VecMs, VecEbv, log10(VecFrho), Vect0, ht0, chi2values(peakInd(peakNumber)).all.t0);
savefig(ht0,sprintf('%s_t_ref.fig',figs_name));
hMs = plotMargProbMs(marginalizedPDF, VecMs, VecEbv, log10(VecFrho), Vect0, hMs, chi2values(peakInd(peakNumber)).all.Ms);
savefig(hMs,sprintf('%s_Ms.fig',figs_name));
hEbv = plotMargProbEbv(marginalizedPDF, VecMs, VecEbv, log10(VecFrho), Vect0, hEbv, chi2values(peakInd(peakNumber)).all.Ebv);
savefig(hEbv,sprintf('%s_Ebv.fig',figs_name));
hfrho = plotMargProbFrho(marginalizedPDF, VecMs, VecEbv, VecFrho, Vect0, hfrho, chi2values(peakInd(peakNumber)).all.frho);    
savefig(hfrho,sprintf('%s_frho.fig',figs_name));


fprintf('\n');
line = 'Parameter          ';
for iPeak=1:nPeaks,line = [line sprintf('& Peak \\\\#%d                    ',iPeak)];end
line = [line '\\\\\n'];
fprintf(line);
fprintf('\\hline\n');

line='$R_{*}$ [$R_\\sun$] ';
for iPeak=1:nPeaks 
    line = [line sprintf('& $%4.0f_{-%4.0f}^{+%4.0f}$      ', chi2values(peakInd(iPeak)).all.Rs, chi2values(peakInd(iPeak)).all.RssigmaM, chi2values(peakInd(iPeak)).all.RssigmaP)];
end
line = [line '\\\\\n'];
fprintf(line);

line ='$v_{s*,8.5}$       ';
for iPeak=1:nPeaks
    line = [line sprintf('& $%4.2f_{-%4.2f}^{+%4.2f}$      ', chi2values(peakInd(iPeak)).all.Vs/10^8.5, chi2values(peakInd(iPeak)).all.VssigmaM, chi2values(peakInd(iPeak)).all.VssigmaP)];
end
line = [line '\\\\\n'];
fprintf(line);

line ='$M_{ej}$ [$M_\\sun$]';
for iPeak=1:nPeaks
    line = [line sprintf('& $%4.1f_{-%4.1f}^{+%4.1f}$      ', chi2values(peakInd(iPeak)).all.Ms, chi2values(peakInd(iPeak)).all.MssigmaM, chi2values(peakInd(iPeak)).all.MssigmaP)];
end
line = [line '\\\\\n'];
fprintf(line);

line ='$t_{\\rm ref}$ [MJD]';
for iPeak=1:nPeaks
    line = [line sprintf('& $%6.2f_{-%4.2f}^{+%4.2f}$  ', chi2values(peakInd(iPeak)).all.t0, chi2values(peakInd(iPeak)).all.t0sigmaM, chi2values(peakInd(iPeak)).all.t0sigmaP)];
end
line = [line '\\\\\n'];
fprintf(line);

line ='$f_{\\rho}$         ';
for iPeak=1:nPeaks
    line = [line sprintf('& $%5.3f_{-%5.3f}^{+%5.3f}$   ', chi2values(peakInd(iPeak)).all.frho, chi2values(peakInd(iPeak)).all.frhosigmaM, chi2values(peakInd(iPeak)).all.frhosigmaP)];
end
line = [line '\\\\\n'];
fprintf(line);

line ='$E_{B-V}$          ';
for iPeak=1:nPeaks
    line = [line sprintf('& $%5.3f_{-%5.3f}^{+%5.3f}$   ', chi2values(peakInd(iPeak)).all.Ebv, chi2values(peakInd(iPeak)).all.EbvsigmaM, chi2values(peakInd(iPeak)).all.EbvsigmaP)];
end
line = [line '\\\\\n'];
fprintf(line);

line ='$\\chi^2$/dof       ';
for iPeak=1:nPeaks
    line = [line sprintf('& %6.2f/%3u                  ',chi2values(peakInd(iPeak)).all.chi2, chi2values(peakInd(iPeak)).all.dof)];
end
line = [line '\\\\\n'];
fprintf(line);

fprintf('\\hline\n')
end

function [f, Vs, Rs, iVs, iRs]=plotPDF(PDFmap, validmap, VecVs, VecRs)
%===================================================
% plot raw PDF
%===================================================
PDFmap = PDFmap./(trapz(VecRs,trapz(VecVs/10^8.5,PDFmap,2)));
[~, imax] = max(PDFmap(:));
[iRs, iVs] = ind2sub(size(PDFmap),imax);
Vs = VecVs(iVs);
Rs = VecRs(iRs);
f=figure;
f.Name ='PDF';
minPDF = min(PDFmap(validmap));
maxPDF = max(PDFmap(validmap));

PDFmap(~validmap) = -1;
hac = subplot(4,4,[1:3 5:7 9:11]);
hac2 = axes('Position',hac.Position);
hac2.Visible = 'off';
cmap = colormap('gray');cmap=cmap(end:-1:1,:);
[~,hcf]=contourf(hac,VecVs/10^8.5, VecRs, PDFmap,linspace(minPDF,maxPDF,64),'LineStyle','none'); colormap(hac,cmap);
hc=contour(hac2,VecVs/10^8.5,VecRs,PDFmap,linspace(minPDF,maxPDF,20));
hac2.Visible = 'off';
colormap(hac2,'default');
linkaxes([hac hac2]);

set(hac,'XAxisLocation','top');
grid minor
set(hac,'XScale','log','YScale','log');
set(hac2,'XScale','log','YScale','log');
xlabel(hac,'$v_{s*,8.5}$','Interpreter','Latex');
ylabel(hac,'$R/R_{\odot}$','Interpreter','Latex');
set(hac2,'CLimMode','Manual');
hold(hac2,'on');
contourf(hac2,VecVs/10^8.5,VecRs,1e-300.*~validmap,[1e-300 1e-300], 'k');
haR=subplot(4,4,[4 8 12]);
hR=semilogy(trapz(VecVs/10^8.5,(validmap.*PDFmap),2),VecRs);
grid minor
linkaxes([hac hac2 haR],'y')
xlabel('PDF [$R_\odot/R$]','Interpreter','Latex');
haV=subplot(4,4,[13:15]);
hV=semilogx(VecVs/10^8.5,trapz(VecRs,(validmap.*PDFmap),1));
grid minor
linkaxes([hac hac2 haV],'x')
ylabel('PDF $[v_{s*}^{-1}]$','Interpreter','Latex');
xlabel('$v_{s*,8.5}$','Interpreter','Latex');
ylim(hac,[VecRs(1) VecRs(end)]);
hac2.Visible = 'off';
end

function [f, MaxVs, MaxRs, SigmaVs, SigmaRs, Vprob, Rprob]=plotCDF(PDFmap, validmap, VecVs, VecRs)
%=============================================
% plot PDF contour which donate to 1,2,3 sigma
% using cell area
%=============================================
f=figure;
f.Name = 'CDF';
dVs = diff(VecVs(:)).';
dRs  = diff(VecRs(:)).'; 
[MatVs, MatRs] = meshgrid(([0, dVs] + [dVs, 0])/2,([0, dRs] + [dRs, 0])/2);
area = MatVs.* MatRs;
PDFmap(~validmap) = 0;
PDFmap = PDFmap./(trapz(VecRs,trapz(VecVs,PDFmap,2)));

[sorted_pdf, idx] = sort(PDFmap(validmap),'descend');
valid_area = area(validmap);
sorted_dCDF = sorted_pdf.*valid_area(idx);
cdf = cumsum(sorted_dCDF);
tmpPDF = PDFmap(validmap);
for sigma = 1:3
    prob = normcdf(sigma,0,1)-normcdf(-sigma,0,1);
    cdfid = find(cdf>=prob,1,'first');   
    v(sigma) = tmpPDF(idx(cdfid));
end
hac = subplot(4,4,[1:3 5:7 9:11]);
contour(VecVs/10^8.5,VecRs,PDFmap,v,'-k');
set(hac,'XAxisLocation','top');
set(gca,'XScale','log','YScale','log');
set(gca,'CLimMode','Manual','CLim',[-v(2) v(1)]);
ylabel('$R/R_{\odot}$','Interpreter','Latex');
if (hac.XAxis.TickValues(1)~=hac.XAxis.Limits(1))
    hac.XAxis.TickValues = [hac.XAxis.Limits(1) hac.XAxis.TickValues];
end
if (hac.XAxis.TickValues(end) ~= hac.XAxis.Limits(end))
    hac.XAxis.TickValues = [hac.XAxis.TickValues hac.XAxis.Limits(end)];
end
Ticks = [hac.XAxis.TickValues(1)];
for itick = 1:length(hac.XAxis.TickValues)-1
    if (hac.XAxis.TickValues(itick+1)/hac.XAxis.TickValues(itick))>5
        mag = floor(log10(hac.XAxis.TickValues(itick)));
        Ticks = [Ticks 5*10^mag hac.XAxis.TickValues(itick+1)];
    else
        Ticks = [Ticks hac.XAxis.TickValues(itick+1)];
    end
end
hac.XAxis.TickValues = Ticks;
if (hac.YAxis.TickValues(1)~=hac.YAxis.Limits(1))
    hac.YAxis.TickValues = [hac.YAxis.Limits(1) hac.YAxis.TickValues];
end
if (hac.YAxis.TickValues(end) ~= hac.YAxis.Limits(end))
    hac.YAxis.TickValues = [hac.YAxis.TickValues hac.YAxis.Limits(end)];
end
Ticks = [hac.YAxis.TickValues(1)];
for itick = 1:length(hac.YAxis.TickValues)-1
    if (hac.YAxis.TickValues(itick+1)/hac.YAxis.TickValues(itick))>4
        mag = floor(log10(hac.YAxis.TickValues(itick)));
        Ticks = [Ticks 5*10^mag hac.YAxis.TickValues(itick+1)];
    else
        Ticks = [Ticks hac.YAxis.TickValues(itick+1)];
    end
end
hac.YAxis.TickValues = Ticks;
hold on
PDFmap(~validmap) = -1;
c=contourc(VecVs/10^8.5,VecRs,PDFmap,[-1 -1]);
idx = 1;
csize = size(c,2);
while (idx<csize)
    vertices = c(2,idx);
    xpoints = c(1,idx+[1:vertices]);
    ypoints = c(2,idx+[1:vertices]);
    %close the patches
    if (xpoints(1) == min(VecVs/10^8.5))
        if (ypoints(end) == max(VecRs))
            xpoints = [xpoints min(VecVs/10^8.5)];
            ypoints = [ypoints max(VecRs)];
        elseif (ypoints(end) == min(VecRs))
            xpoints = [xpoints min(VecVs/10^8.5)];
            ypoints = [ypoints min(VecRs)];
        elseif (xpoints(end) == max(VecVs/10^8.5))
            xpoints = [xpoints max(VecVs/10^8.5) min(VecVs/10^8.5)];
            ypoints = [ypoints max(VecRs)        max(VecRs)       ];
        elseif (xpoints(end) == min(VecVs/10^8.5))
            xpoints = [xpoints min(VecVs/10^8.5) max(VecVs/10^8.5) max(VecVs/10^8.5) min(VecVs/10^8.5)];
            ypoints = [ypoints max(VecRs)        max(VecRs)        min(VecRs)        min(VecRs)       ];
        end
    elseif (ypoints(1) == max(VecRs))
        if (xpoints(end) == min(VecVs/10^8.5))
            xpoints = [xpoints min(VecVs/10^8.5)];
            ypoints = [ypoints max(VecRs)];
        elseif (xpoints(end) == max(VecVs/10^8.5))
            xpoints = [xpoints max(VecVs/10^8.5) min(VecVs/10^8.5) min(VecVs/10^8.5)];
            ypoints = [ypoints min(VecRs)        min(VecRs)        max(VecRs)];
        end
    elseif (ypoints(1) == min(VecRs))
        if (xpoints(end) == min(VecVs/10^8.5))
%             xpoints = [xpoints min(VecVs/10^8.5) max(VecVs/10^8.5)   max(VecVs/10^8.5)];
%             ypoints = [ypoints max(VecRs)        max(VecRs)          min(VecRs)        ];
            xpoints = [xpoints min(VecVs/10^8.5)];
            ypoints = [ypoints min(VecRs)       ];
        end
    end
    patch(xpoints,ypoints,[0.7 0.7 0.7]);
    idx = idx+1+vertices;
end
grid minor
hac.YLim = [VecRs(1) VecRs(end)];
haR=subplot(4,4,[4 8 12]);
Rprob = trapz(VecVs,(validmap.*PDFmap),2);
Rprob = Rprob./trapz(VecRs,Rprob);
hR=semilogy(Rprob,VecRs,'-k');
linkaxes([hac,haR],'y');
hac.YAxis.Limits = [VecRs(1) VecRs(end)];
haR.YAxis.TickValues = hac.YAxis.TickValues;
haR.YAxis.TickLabels = hac.YAxis.TickLabels;
xlabel('PDF [$R_\odot/R$]','Interpreter','Latex');
hold on
grid minor

[downR,downRval,upR,upRval]=oneSigmaMove(VecRs,Rprob);
line([0 downRval], downR.*[1 1],'Color','k');
line([0 upRval], upR.*[1 1],'Color','k');

[~, imax] = max(Rprob);
MaxRs   = VecRs(imax);
SigmaRs = [downR upR];

set(haR,'YAxisLocation','right','XTickLabelRotation',-90)
haR.XAxis.Exponent = 0;
haV=subplot(4,4,[13:15]);
Vprob = trapz(VecRs,(validmap.*PDFmap),1);
Vprob = Vprob./trapz(VecVs/10^8.5,Vprob);
hV=semilogx(VecVs/10^8.5,Vprob,'-k');
linkaxes([hac,haV],'x');
haV.XAxis.TickValues = hac.XAxis.TickValues;
haV.XAxis.TickLabels = hac.XAxis.TickLabels;
ylabel('PDF $[v_{s*}^{-1}]$','Interpreter','Latex');
xlabel('$v_{s*,8.5}$','Interpreter','Latex');
hold on
grid minor

[downV,downVval,upV,upVval]=oneSigmaMove(VecVs/10^8.5,Vprob);
line(downV.*[1 1], [0 downVval], 'Color','k');
line(upV.*[1 1], [0 upVval], 'Color','k');

[~, imax] = max(Vprob);
MaxVs   = VecVs(imax);
SigmaVs = [downV upV];
end

function f=plotPvalue(map, VecVs, VecRs)
%===================================================
% plot P value
%===================================================
validmap = map.valid_VtoR;
[~, idx] = min(map.pvalue,[],1);

Pmap = zeros(size(squeeze(map.pvalue(1,:,:))));
for iVs=1:length(VecVs),for iRs=1:length(VecRs),Pmap(iRs,iVs)=map.pvalue(idx(1,iRs,iVs),iRs,iVs);end,end

f=figure;
f.Name ='P-value';

hc=contour(VecVs/10^8.5,VecRs,Pmap,[0.01 0.05]);% colorbar;
set(gca,'XScale','log','YScale','log');
xlabel('$v_{s*,8.5}$','Interpreter','Latex');
ylabel('$R/R_{\odot}$','Interpreter','Latex');
hold on
contourf(VecVs/10^8.5,VecRs,1e-9.*~validmap,[1e-9 1e-9], 'k');

end

function [f, t0PDF, Sigmat0]=plotMargProbt0(PDFmap, VecMs, VecEbv, Vecfrho, Vect0, fig, maxt0)
%===============================================
% plot t0 PDF
%===============================================
if (nargin>=6)&&(~isempty(fig))
    f=fig;
    LineStyle = '-';
    figure(f);
else
    f=figure;
    f.Name = 'P(t0)';
    LineStyle = '--';
end

if length(Vecfrho)>1
    t0PDF = trapz(Vecfrho,trapz(VecMs,trapz(VecEbv,PDFmap,4),2),3);
else
    t0PDF = trapz(VecMs,trapz(VecEbv,PDFmap,3),2);
end
t0PDF  = t0PDF./trapz(Vect0,t0PDF);
plot(Vect0,t0PDF,'LineStyle',LineStyle,'Color','k');
hold on

ax = gca;
ax.XAxis.Exponent = 0;

if (nargin==7)
    t0 = maxt0;
    [downt0,downt0val,upt0,upt0val]=oneSigmaMove(Vect0,t0PDF,maxt0);
else
    [~,imax] = max(t0PDF);
    t0 = Vect0(imax);
    [downt0,downt0val,upt0,upt0val]=oneSigmaMove(Vect0,t0PDF);
end
line((downt0).*[1 1], [0 downt0val], 'Color','k','LineStyle',LineStyle);
line((upt0).*[1 1], [0 upt0val], 'Color','k','LineStyle',LineStyle);

Sigmat0=[downt0 upt0];

xlabel('$t_{\rm ref}\, [MJD]$','Interpreter','Latex');
ylabel('$\rm PDF\, [dy^{-1}]$','Interpreter','Latex');
drawnow
end

function [f, MsPDF, SigmaMs]=plotMargProbMs(PDFmap, VecMs, VecEbv, Vecfrho, Vect0, fig, maxMs)
%===============================================
% plot Ms PDF
%===============================================
if (nargin>=6)&&(~isempty(fig))
    f=fig;
    LineStyle = '-';
    figure(f);
else
    f=figure;
    f.Name = 'P(Ms)';
    LineStyle = '--';
end

if length(Vecfrho)>1
    MsPDF = trapz(Vecfrho,trapz(Vect0,trapz(VecEbv,PDFmap,4),1),3);
else
    MsPDF = trapz(Vect0,trapz(VecEbv,PDFmap,3),1);
end
MsPDF  = squeeze(MsPDF./trapz(VecMs,MsPDF));
plot(VecMs,MsPDF,'LineStyle',LineStyle,'Color','k');
hold on

if (nargin==7)
    [downMs,downMsval,upMs,upMsval]=oneSigmaMove(VecMs,MsPDF,maxMs);
    Ms = maxMs;
else
    [downMs,downMsval,upMs,upMsval]=oneSigmaMove(VecMs,MsPDF);
    [~,imax] = max(MsPDF);
    Ms = VecMs(imax);
end
line(downMs.*[1 1], [0 downMsval], 'Color','k','LineStyle',LineStyle);
line(upMs.*[1 1], [0 upMsval], 'Color','k','LineStyle',LineStyle);

SigmaMs = [downMs upMs];

xlabel('$\rm M_{ej}\, [M_{\odot}]$','Interpreter','Latex');
ylabel('$\rm PDF\, [M_{\odot}^{-1}]$','Interpreter','Latex');
drawnow

end

function [f, EbvPDF, SigmaEbv]=plotMargProbEbv(PDFmap, VecMs, VecEbv, Vecfrho, Vect0, fig, maxEbv)
%===============================================
% plot Ebv PDF
%===============================================
if (nargin>=6)&&(~isempty(fig))
    f=fig;
    LineStyle = '-';
    figure(f);
else
    f=figure;
    f.Name = 'P(Ebv)';
    LineStyle = '--';
end

if length(Vecfrho)>1
    EbvPDF = trapz(Vecfrho,trapz(VecMs,trapz(Vect0,PDFmap,1),2),3);
else
    EbvPDF = trapz(VecMs,trapz(Vect0,PDFmap,1),2);
end
EbvPDF  = squeeze(EbvPDF./trapz(VecEbv,EbvPDF));
plot(VecEbv,EbvPDF,'LineStyle',LineStyle,'Color','k');
hold on

if (nargin==7)
    [downEbv,downEbvval,upEbv,upEbvval]=oneSigmaMove(VecEbv,EbvPDF,maxEbv);
    Ebv = maxEbv;
else
    [downEbv,downEbvval,upEbv,upEbvval]=oneSigmaMove(VecEbv,EbvPDF);
    [~,imax] = max(EbvPDF);
    Ebv = VecEbv(imax);
end
line(downEbv.*[1 1], [0 downEbvval], 'Color','k','LineStyle',LineStyle);
line(upEbv.*[1 1], [0 upEbvval], 'Color','k','LineStyle',LineStyle);

SigmaEbv = [downEbv upEbv];

xlabel('$\rm E_{B-V}$','Interpreter','Latex');
ylabel('$\rm PDF\, [E_{B-V}^{-1}]$','Interpreter','Latex');
drawnow

end


function [f, frhoPDF, SigmaFrho]=plotMargProbFrho(PDFmap, VecMs, VecEbv, Vecfrho, Vect0, fig, maxfrho)
%===============================================
% plot frho PDF
%===============================================
if (nargin>=6)&&(~isempty(fig))
    f=fig;
    LineStyle = '-';
    figure(f);
else
    f=figure;
    f.Name = 'P(f_rho)';
    LineStyle = '--';
end

frhoPDF = trapz(VecEbv,trapz(VecMs,trapz(Vect0,PDFmap,1),2),4);
frhoPDF  = squeeze(frhoPDF./trapz(log10(Vecfrho),frhoPDF));
semilogx(Vecfrho,frhoPDF,'LineStyle',LineStyle,'Color','k');
hold on

if (nargin==7)
    [downFrho,downFrhoval,upFrho,upFrhoval]=oneSigmaMove(log10(Vecfrho),frhoPDF,log10(maxfrho));
    frho = maxfrho;
else
    [downFrho,downFrhoval,upFrho,upFrhoval]=oneSigmaMove(log10(Vecfrho),frhoPDF);
    [~,imax] = max(frhoPDF);
    frho = Vecfrho(imax);
    
end
line(10.^downFrho.*[1 1], [0 downFrhoval], 'Color','k','LineStyle',LineStyle);
line(10.^upFrho.*[1 1], [0 upFrhoval], 'Color','k','LineStyle',LineStyle);

SigmaFrho = 10.^[downFrho upFrho];

xlabel('$\rm f_{\rho}$','Interpreter','Latex');
ylabel('$\rm PDF\,[1/log_{10}(f_{\rho})]$','Interpreter','Latex');

drawnow

end

function [f, f1]=plotMargQuantity(map, validmap, QuanName, VecQuan, VecVs, VecRs)
%==================================================
% plot marginalized quantities
%===================================================
[~,maxidx] = max(map,[],3);

maxQuan = zeros(size(maxidx));
for iRs=1:length(VecRs),for iVs=1:length(VecVs),maxQuan(iRs,iVs) = VecQuan(maxidx(iRs,iVs));end,end

f=figure;
f.Name = sprintf('%s for max probability',QuanName);
hc=contourf(VecVs/10^8.5,VecRs,maxQuan); colorbar;
grid minor
set(gca,'XScale','log','YScale','log');
xlabel('$v_{s*,8.5}$','Interpreter','Latex');
ylabel(sprintf('$%s$',QuanName),'Interpreter','Latex');
set(gca,'CLimMode','Manual');
hold on
contourf(VecVs/10^8.5,VecRs,1e-9.*~validmap,[1e-9 1e-9], 'k');

f1=figure;
f1.Name = sprintf('%s for max probability as function of v_s*,8.5',QuanName);
margMap=squeeze(trapz(VecRs,map,1));
[~,maxidx] = max(margMap,[],2);
plot(VecVs/10^8.5,VecQuan(maxidx));
xlabel('$v_{s*,8.5}$','Interpreter','Latex');

end

% this function calculates the one sigma limits for a given variable and
% its PDF.
function [downvar,downval,upvar,upval]=oneSigma(var,PDF)
CDF  = cumtrapz(var,PDF);
[~, maxidx] = max(PDF);
if (CDF(maxidx)-(0.5-normcdf(-1,0,1))>=0)
    if (CDF(maxidx)+(normcdf(1,0,1)-0.5)<=1)
        ind = find(CDF>=CDF(maxidx)-0.5+normcdf(-1,0,1),1,'first')-1;
        downvar = var(ind)+((CDF(maxidx)-0.5+normcdf(-1,0,1))-CDF(ind)) * (var(ind+1)-var(ind))/(CDF(ind+1)-CDF(ind));
        downval = PDF(ind)+(PDF(ind+1)-PDF(ind))/(var(ind+1)-var(ind))*(downvar-var(ind));
        ind   = find(CDF>=CDF(maxidx)-0.5+normcdf(1,0,1),1,'first')-1;
        upvar = var(ind) + ((CDF(maxidx)-0.5+normcdf(1,0,1))-CDF(ind)) * (var(ind+1)-var(ind))/(CDF(ind+1)-CDF(ind));
        upval = PDF(ind) + (PDF(ind+1)-PDF(ind))/(var(ind+1)-var(ind))*(upvar-var(ind));        
    else
        ind = find(CDF>=2*normcdf(-1,0,1),1,'first')-1;
        downvar = var(ind) + ((2*normcdf(-1,0,1))-CDF(ind)) * (var(ind+1)-var(ind))/(CDF(ind+1)-CDF(ind));
        downval = PDF(ind)+(PDF(ind+1)-PDF(ind))/(var(ind+1)-var(ind))*(downvar-var(ind));       
        upvar   = var(end);
        upval   = PDF(end);
    end
else
    downvar = var(1);
    downval = PDF(1);
    ind = find(CDF>=(normcdf(1,0,1)-normcdf(-1,0,1)),1,'first')-1;
    upvar = var(ind) + ((normcdf(1,0,1)-normcdf(-1,0,1))-CDF(ind)) * (var(ind+1)-var(ind))/(CDF(ind+1)-CDF(ind));
    upval = PDF(ind) + (PDF(ind+1)-PDF(ind))/(var(ind+1)-var(ind))*(upvar-var(ind));        
end
end

function [downvar,downval,upvar,upval]=oneSigmaMove(var,PDF,max_var)
CDF  = cumtrapz(var,PDF);
%minlen = inf;
%index  = 1;
%dir    = 'down';

% find the shortest interval with CDF of normcdf(1,0,1)-normcdf(-1,0,1) on
% the discerte grid of var:
downlen = inf*ones(length(var),1);
uplen   = inf*ones(length(var),1);
for i=1:length(CDF)
    if (CDF(i)<2*normcdf(-1,0,1))
        ind = find(CDF>=CDF(i)+normcdf(1,0,1)-normcdf(-1,0,1),1,'first')-1;
        up = var(ind) + (CDF(i)+(normcdf(1,0,1)-normcdf(-1,0,1))-CDF(ind)) * (var(ind+1)-var(ind))/(CDF(ind+1)-CDF(ind));
        if (nargin==2)||((max_var>=var(i))&&(max_var<=up))
            downlen(i)=up-var(i);             
%         elseif (nargin>2)&&((max_var>up)&&(max_var<=var(ind+1)))
%             up = max_var;
%             CDF_max = CDF(ind) + (CDF(ind+1)-CDF(ind))/(var(ind+1)-var(ind))*(max_var-var(ind));
%             down = var(i) + (CDF_max-(normcdf(1,0,1)-normcdf(-1,0,1))-CDF(i)) * (var(i+1)-var(i))/(CDF(i+1)-CDF(i));
%             downlen(i) = up-down;
        end
    else
        % we reached the upper bound of the probability CDF(i)+1\sigma>1;
        break;
    end
end
for i=length(CDF):-1:1
    if ((1-CDF(i))<2*normcdf(-1,0,1))
        ind = find(CDF>=CDF(i)-(normcdf(1,0,1)-normcdf(-1,0,1)),1,'first')-1;
        down = var(ind) + (CDF(i)-(normcdf(1,0,1)-normcdf(-1,0,1))-CDF(ind)) * (var(ind+1)-var(ind))/(CDF(ind+1)-CDF(ind));
        if (nargin==2)||((max_var>=down)&&(max_var<=var(i)))
            uplen(i)=var(i)-down;             
%         elseif (nargin>2)&&((max_var>=var(ind))&&(max_var<down))
%             down = max_var;
%             CDF_max = CDF(ind) + (CDF(ind+1)-CDF(ind))/(var(ind+1)-var(ind))*(max_var-var(ind));
%             up = var(i-1) + (CDF_max+(normcdf(1,0,1)-normcdf(-1,0,1))-CDF(i-1)) * (var(i)-var(i-1))/(CDF(i)-CDF(i-1));
%             uplen(i) = up-down;
        end
    else
        % we reached the lower bound of the probability 1-CDF(i)<1\sigma;
        break;
    end
end

[min_down_len,idown]=min(downlen);
[min_up_len,iup]=min(uplen);

if (min_down_len<min_up_len)
    dir = 'down';
    index = idown;
    findind = find(CDF>=CDF(idown)+normcdf(1,0,1)-normcdf(-1,0,1),1,'first')-1;
    up = var(findind) + (CDF(idown)+(normcdf(1,0,1)-normcdf(-1,0,1))-CDF(findind)) * (var(findind+1)-var(findind))/(CDF(findind+1)-CDF(findind));
    if (nargin==2)||((max_var>=var(idown))&&(max_var<=up))
        minlen=up-var(idown);
        upvar=up;
    else
        pause;
%     elseif (nargin>2)&&((max_var>up)&&(max_var<=var(findind+1)))
%         dir = 'up';
%         index = findind+1;
%         findind = idown;
%         down = var(findind) + (CDF(index)-(normcdf(1,0,1)-normcdf(-1,0,1))-CDF(findind)) * (var(findind+1)-var(findind))/(CDF(findind+1)-CDF(findind));
%         minlen=var(index)-down;             
    end    
else
    dir = 'up';
    index = iup;
    findind = find(CDF<CDF(iup)-(normcdf(1,0,1)-normcdf(-1,0,1)),1,'last');
    down = var(findind) + (CDF(iup)-(normcdf(1,0,1)-normcdf(-1,0,1))-CDF(findind)) * (var(findind+1)-var(findind))/(CDF(findind+1)-CDF(findind));
    if (nargin==2)||((max_var>=down)&&(max_var<=var(iup)))
        minlen=var(iup)-down;
        downvar=down;
    else
        pause;
%     elseif (nargin>2)&&((max_var<down)&&(max_var>=var(findind)))
%         dir = 'down';
%         index = findind;
%         findind = iup-1;
%         up = var(findind) + (CDF(index)+(normcdf(1,0,1)-normcdf(-1,0,1))-CDF(findind)) * (var(findind+1)-var(findind))/(CDF(findind+1)-CDF(findind));
%         upvar = up;
%         minlen = up-var(index);            
    end
end

% fine tune the result between grid points
if strcmp(dir,'down')
    % the lower end of the interval is a grid point
    PDF_L = PDF(index);
    PDF_H = PDF(findind) + (PDF(findind+1)-PDF(findind))/(var(findind+1)-var(findind))*(upvar-var(findind));
    if (PDF_L<PDF_H) % dxL,dxH>0
        dPDF_L = (PDF(index+1)-PDF(index))/(var(index+1)-var(index));
        dPDF_H = (PDF(findind+1)-PDF(findind))/(var(findind+1)-var(findind));
        if (dPDF_L==dPDF_H)
            % if the derivatives are identical the solution is degenerated.
            % Do not move on this case;
            dxL    = 0;
            dxH    = 0;
        else
            dxL   = PDF_L./dPDF_L.*(-1+sqrt(1-dPDF_L./(dPDF_L-dPDF_H).*(1-PDF_H.^2./PDF_L.^2)));
            dxH    = (-PDF_H+sqrt(PDF_H^2+2.*dPDF_H.*(PDF_L.*dxL+0.5.*dPDF_L.*dxL.^2)))./dPDF_H;
        end
        if (dxL>=0)
            ind     = findind;
            downvar = var(index);            
            while (dxH > (var(ind+1)-upvar))
                % as long as the required movement of the upper end is
                % beyond var(ind+1) progress to var(ind+1) and calculate
                % the residual required movement:
                dxH    = (var(ind+1)-upvar); 
                dxL    = (-PDF_L+sqrt(PDF_L^2+2.*dPDF_L.*(PDF_H.*dxH+0.5.*dPDF_H.*dxH.^2)))./dPDF_L;
                if (nargin==3)&&(downvar+dxL>=max_var)
                    % if the movement will skip the constraint variable
                    % value limit the movement so the lower end of the
                    % movement will be the constraint value.
                    dxL = max_var-downvar;
                    dxH = (-PDF_H+sqrt(PDF_H^2+2.*dPDF_H.*(PDF_L.*dxL+0.5.*dPDF_L.*dxL.^2)))./dPDF_H;
                else
                    % move to edge of current cell
                    upvar   = var(ind+1);
                    PDF_H   = PDF(ind+1);
                    downvar = downvar+dxL;
                    PDF_L = PDF_L + dPDF_L * dxL;
                    ind = ind + 1;
                    if (length(var)>ind)
                        % if we did not reach the end of the array
                        % calculate the residual required movement.
                        dPDF_H = (PDF(ind+1)-PDF(ind))/(var(ind+1)-var(findind));
                        dxL    = PDF_L./dPDF_L.*(-1+sqrt(1-dPDF_L./(dPDF_L-dPDF_H).*(1-PDF_H.^2./PDF_L.^2)));
                        dxH    = (-PDF_H+sqrt(PDF_H^2+2.*dPDF_H.*(PDF_L.*dxL+0.5.*dPDF_L.*dxL.^2)))./dPDF_H;
                        if (dxL<0)
                            % the function is not monotonic, the point
                            % between the two cell is a local minimum.
                            dxL=0;
                            dxH=0;
                        end                                                       
                    else
                        % we reached the end of the array, no further
                        % movement is possible.
                        dxL=0;
                        dxH=0;
                    end
                end
            end
            % move the required dxL and dxH, which is either the original
            % required value or the residual one from the last cell edge.
            if (nargin==3)&&(downvar+dxL>=max_var)
                % if the movement will skip the constraint variable
                % value limit the movement so the lower end of the
                % movement will be the constraint value.
                dxL = max_var-downvar;
                dxH = (-PDF_H+sqrt(PDF_H^2+2.*dPDF_H.*(PDF_L.*dxL+0.5.*dPDF_L.*dxL.^2)))./dPDF_H;
            end
            upvar   = upvar+dxH;
            upval   = PDF_H + dPDF_H * dxH;
            downvar = downvar+dxL;
            downval = PDF_L + dPDF_L * dxL;
        else % dxL<0 , PDF_L<PDF_H
            % while PDF_L<PDF_H and the movement should be rightward, the
            % calculated dxL resulted negative, means leftward movement. We
            % should not arrive to this point, if we arrive do not move.
            downvar = var(index);
            downval = PDF(index);
            upval   = PDF(findind) + (PDF(findind+1)-PDF(findind))/(var(findind+1)-var(findind))*(upvar-var(findind));                            
        end
    else            % dxL, dxH<0
        if (index>1)
            % recalculate the slope of the lower end and the result
            % required movement, since it should move to the previous
            % cell.
            % the dx signs are not opposite, movement leftward is
            % defined positive.
            dPDF_L = (PDF(index)-PDF(index-1))/(var(index)-var(index-1));
            dPDF_H = (PDF(findind+1)-PDF(findind))/(var(findind+1)-var(findind));
            if (dPDF_L==dPDF_H)
                % if the derivatives are identical the solution is degenerated.
                % Do not move on this case;
                dxL    = 0;
                dxH    = 0;
            else
                dxL = PDF_L./dPDF_L .* (1 - sqrt(1-dPDF_L.*((1-PDF_H.^2./PDF_L.^2)./(dPDF_L-dPDF_H))));
                dxH = (PDF_H-sqrt(PDF_H.^2-2.*dPDF_H.*(PDF_L.*dxL-0.5.*dPDF_L.*dxL.^2)))./dPDF_H;
            end
            if (dxL>=0)
                ind = findind;
                downvar = var(index);            
                while ((upvar-dxH) < var(ind))
                    % as long as the required movement is larger than
                    % one cell, move to the cell edge and calculate
                    % the residual required movement
                    dxH = upvar-var(ind);
                    dxL = (PDF_L-sqrt(PDF_L.^2-2.*dPDF_L.*(PDF_H.*dxH-0.5.*dPDF_H.*dxH.^2)))./dPDF_L;
                    if (nargin==3)&&((upvar-dxH)<=max_var)
                        % if the residual required movement takes the
                        % constraint value outside the interval, limit
                        % the movement such that the upper end of the
                        % interval will be the constaint value.
                        dxH = upvar-max_var;
                        dxL = (PDF_L-sqrt(PDF_L.^2-2.*dPDF_L.*(PDF_H.*dxH-0.5.*dPDF_H.*dxH.^2)))./dPDF_L;
                    else
                        % move to edge of current cell
                        upvar   = var(ind);
                        PDF_H   = PDF(ind);
                        downvar = downvar-dxL;
                        PDF_L = PDF_L - dPDF_L * dxL;
                        ind = ind - 1;

                        dPDF_H = (PDF(ind)-PDF(ind-1))/(var(ind)-var(ind-1));
                        dxL = PDF_L./dPDF_L .* (1 - sqrt(1-dPDF_L.*((1-PDF_H.^2./PDF_L.^2)./(dPDF_L-dPDF_H))));
                        dxH = (PDF_H-sqrt(PDF_H.^2-2.*dPDF_H.*(PDF_L.*dxL2-0.5.*dPDF_L.*dxL.^2)))./dPDF_H;
                        if dxL<0
                            % the function is not monotonic at the
                            % upper end. This point is the local
                            % minimum.
                            dxL=0;
                            dxH=0;
                        end
                    end
                end
                % move the required dxL and dxH, which is either the original
                % required value or the residual one from the last cell edge.
                if (nargin==3)&&(upvar-dxH<=max_var)
                    % if the movement will skip the constraint variable
                    % value limit the movement so the upper end of the
                    % movement will be the constraint value.
                        dxH = upvar-max_var;
                        dxL = (PDF_L-sqrt(PDF_L.^2-2.*dPDF_L.*(PDF_H.*dxH-0.5.*dPDF_H.*dxH.^2)))./dPDF_L;
                end
                upvar   = upvar - dxH;
                upval   = PDF_H - dPDF_H * dxH;
                downvar = downvar - dxL;
                downval = PDF_L - dPDF_L * dxL;

            else
                % the function is not monotonic in its lower end. This
                % the local minimum of the function, no movement is
                % required.
                downvar = var(index);
                downval = PDF(index);
                upval   = PDF(findind) + (PDF(findind+1)-PDF(findind))/(var(findind+1)-var(findind))*(upvar-var(findind));                
            end            
        else %(index==1)
            % we are at the edge of the array, no movement is possible.
            downvar = var(index);
            downval = PDF(index);
            upval   = PDF(findind) + (PDF(findind+1)-PDF(findind))/(var(findind+1)-var(findind))*(upvar-var(findind));
        end
    end
    
else %strcmp(dir,'up')
    % the upper end of the interval is a grid point
    PDF_L = PDF(findind) + (PDF(findind+1)-PDF(findind))/(var(findind+1)-var(findind))*(downvar-var(findind));
    PDF_H = PDF(index);
    if (PDF_L>PDF_H) 
        % leftward movement
        dPDF_L = (PDF(findind+1)-PDF(findind))/(var(findind+1)-var(findind));
        dPDF_H = (PDF(index)-PDF(index-1))/(var(index)-var(index-1));
        if (dPDF_L==dPDF_H)
            % if the derivatives are identical the solution is degenerated.
            % Do not move on this case;
            dxL    = 0;
            dxH    = 0;
        else
            dxL = PDF_L./dPDF_L .* (1 - sqrt(1-dPDF_L.*((1-PDF_H.^2./PDF_L.^2)./(dPDF_L-dPDF_H))));
            dxH = (PDF_H-sqrt(PDF_H.^2-2.*dPDF_H.*(PDF_L.*dxL-0.5.*dPDF_L.*dxL.^2)))./dPDF_H;
        end
        if (dxL>=0)
            ind   = findind;
            upvar = var(index);
            while (dxL > (downvar-var(ind)))
                % as long as the required movement of the lower end is
                % beyond var(ind) progress to var(ind) and calculate the
                % residual required movement:
                dxL=downvar-var(ind);
                dxH = (PDF_H-sqrt(PDF_H.^2-2.*dPDF_H.*(PDF_L.*dxL-0.5.*dPDF_L.*dxL.^2)))./dPDF_H;
                if (nargin==3)&&((upvar-dxH)<=max_var)
                    % if the movement will skip the constraint variable
                    % value limit the movement so the upper end of the
                    % movement will be the constraint value.
                    dxH = upvar-max_var;
                    dxL = (PDF_L-sqrt(PDF_L.^2-2.*dPDF_L.*(PDF_H.*dxH-0.5.*dPDF_H.*dxH.^2)))./dPDF_L;
                else
                    % move to edge of current cell
                    upvar = upvar - dxH;
                    PDF_H = PDF_H - dPDF_H * dxH;
                    downvar = downvar - dxL;
                    PDF_L = PDF_L - dPDF_L * dxL;
                    ind = ind - 1;                    
                    if (ind>1)
                        % if we did not reach the end of the array
                        % calculate the residual required movement.
                        dPDF_L = (PDF(ind+1)-PDF(ind))/(var(ind+1)-var(ind));
                        dxL = PDF_L./dPDF_L .* (1 - sqrt(1-dPDF_L.*((1-PDF_H.^2./PDF_L.^2)./(dPDF_L-dPDF_H))));
                        dxH = (PDF_H-sqrt(PDF_H.^2-2.*dPDF_H.*(PDF_L.*dxL-0.5.*dPDF_L.*dxL.^2)))./dPDF_H;
                        if (dxL<0)
                            % the function is not monotonic, the point
                            % between the two cell is a local minimum.
                            dxL=0;
                            dxH=0;
                        end                                                       
                    else
                        % we reached the end of the array, no further
                        % movement is possible.
                        dxL=0;
                        dxH=0;
                    end
                end
            end
            % move the required dxL and dxH, which is either the original
            % required value or the residual one from the last cell edge.
            if (nargin==3)&&((upvar-dxH)<=max_var)
                % if the movement will skip the constraint variable
                % value limit the movement so the lower end of the
                % movement will be the constraint value.
                dxH = upvar-max_var;
                dxL = (PDF_L-sqrt(PDF_L.^2-2.*dPDF_L.*(PDF_H.*dxH-0.5.*dPDF_H.*dxH.^2)))./dPDF_L;
            end
            upvar   = upvar - dxH;
            upval   = PDF_H - dPDF_H * dxH;
            downvar = downvar - dxL;
            downval = PDF_L - dPDF_L * dxL;
        else % dxL<0 , PDF_L>PDF_H
            % while PDF_L>PDF_H and the movement should be leftward, the
            % calculated dxL resulted negative, means rightward movement. We
            % should not arrive to this point, if we arrive do not move.
            downval = PDF_L;
            upvar   = var(index);
            upval   = PDF_H;                            
        end
        
    else %(PDF_L<=PDF_H) 
        % rightward movement
        dPDF_L = (PDF(findind+1)-PDF(findind))/(var(findind+1)-var(findind));
        dPDF_H = (PDF(index+1)-PDF(index))/(var(index+1)-var(index));
        if (dPDF_L==dPDF_H)
            % if the derivatives are identical the solution is degenerated.
            % Do not move on this case;
            dxL    = 0;
            dxH    = 0;
        else
            dxL = PDF_L./dPDF_L .* (-1 + sqrt(1-dPDF_L.*((1-PDF_H.^2./PDF_L.^2)./(dPDF_L-dPDF_H))));
            dxH = (-PDF_H+sqrt(PDF_H.^2+2.*dPDF_H.*(PDF_L.*dxL+0.5.*dPDF_L.*dxL.^2)))./dPDF_H;
        end
        if (dxL>0)
            ind   = findind;
            upvar = var(index);
            while (dxL > (var(ind+1)-downvar))
                % as long as the required movement of the lower end is
                % beyond var(ind+1) progress to var(ind+1) and calculate the
                % residual required movement:
                dxL=var(ind+1)-downvar;
                dxH = (-PDF_H+sqrt(PDF_H.^2+2.*dPDF_H.*(PDF_L.*dxL+0.5.*dPDF_L.*dxL.^2)))./dPDF_H;
                if (nargin==3)&&((downvar+dxL)>=max_var)
                    % if the movement will skip the constraint variable
                    % value limit the movement so the lower end of the
                    % movement will be the constraint value.
                    dxL = max_var-downvar;
                    dxH = (-PDF_H+sqrt(PDF_H.^2+2.*dPDF_H.*(PDF_L.*dxL+0.5.*dPDF_L.*dxL.^2)))./dPDF_H;
                else
                    % move to edge of current cell and calculate residual
                    % dx
                    upvar = upvar + dxH;
                    PDF_H = PDF_H + dPDF_H * dxH;
                    downvar = downvar + dxL;
                    PDF_L = PDF_L + dPDF_L * dxL;
                    ind = ind + 1; 
                    
                    dPDF_L = (PDF(ind+1)-PDF(ind))/(var(ind+1)-var(ind));
                    dxL = PDF_L./dPDF_L .* (-1 + sqrt(1-dPDF_L.*((1-PDF_H.^2./PDF_L.^2)./(dPDF_L-dPDF_H))));
                    dxH = (-PDF_H+sqrt(PDF_H.^2+2.*dPDF_H.*(PDF_L.*dxL+0.5.*dPDF_L.*dxL.^2)))./dPDF_H;
                end
                if (dxL<0)
                    % the function is not monotonic, the point
                    % between the two cell is a local minimum.
                    dxL=0;
                    dxH=0;
                end                                                       
            end
            % move the required dxL and dxH, which is either the original
            % required value or the residual one from the last cell edge.
            if (nargin==3)&&((downvar+dxL)>=max_var)
                % if the movement will skip the constraint variable
                % value limit the movement so the lower end of the
                % movement will be the constraint value.
                dxL = max_var-downvar;
                dxH = (-PDF_H+sqrt(PDF_H.^2+2.*dPDF_H.*(PDF_L.*dxL+0.5.*dPDF_L.*dxL.^2)))./dPDF_H;
            end
            upvar   = upvar + dxH;
            upval   = PDF_H + dPDF_H * dxH;
            downvar = downvar + dxL;
            downval = PDF_L + dPDF_L * dxL;
        else % dxL<0 , PDF_L>PDF_H
            % while PDF_L>PDF_H and the movement should be leftward, the
            % calculated dxL resulted negative, means rightward movement. We
            % should not arrive to this point, if we arrive do not move.
            downval = PDF_L;
            upvar   = var(index);
            upval   = PDF_H;                            
        end

    end
            
end

end

