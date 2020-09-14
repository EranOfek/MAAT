function [Vs, Rs, Ms, Ebv, Rv, Vect0, PDFval, t0, chi2valuesArray] = findMaximum(filename,Vs_in,Rs_in,bg,mintpoints)
% Find numerically the maximum likelihood
% Package: AstroUtil.supernove.SOPRANOS
% Description: Find numerically the maximum likelihood
% Input  : - grid file name
%          - Vector of initial conditions of Shock velocity parameter [cm s^-1]
%          - Vector of initial conditions or Progenitor radius [R_\sun]
%            (Vs_in and Rs_in must be of the same length)
%          - background value for each band
%          - minimal transient points constraint for each band
% Output : - Table with file contents
%               
% See also: AstroUtil.supernova.SOPRANOS.calcGrid
% Tested : Matlab 9.5
%     By : Noam Ganot                      Oct 2019
%    URL : http://weizmann.ac.il/home/eofek/matlab/
% Example:
% AstroUtil.supernova.SOPRANOS.readZTFtxt('ZTF18abokyfk_unbinned.txt');
% Reliable: 2
%--------------------------------------------------------------------------

in=matfile(filename);
Vect0  = in.Vect0;
VecMs  = in.VecMs;
VecFrho= in.VecFrho;
VecEbv = in.VecEbv;
VecRs  = in.VecRs;
VecVs  = in.VecVs;
redshift = in.redshift;
ProgType = in.ProgType;
data   = in.bands; 
model = in.model;
transient = [in.transientStart in.transientEnd];
clear in

priors.t0 = Vect0;
priors.Ms = VecMs;
priors.frho = VecFrho;
priors.Ebv = VecEbv;
priors.Rs  = VecRs;
priors.Vs = VecVs;

% postfix = filename(6:end-4);
% if exist(sprintf('fminsearch_%s.mat',postfix),'file')
%     load(sprintf('fminsearch_%s.mat',postfix));
% else
for ipeak = 1:length(Vs_in)
    fprintf('Finding the maximum in the vicinity of peak #%1d, Rs=%4.0f, Vs=%4.2f\n', ipeak, Rs_in(ipeak), Vs_in(ipeak)/10^8.5);
    Rv = 3.08;

    c = parallel.pool.Constant(data);
    PDFfun = @(x) (AstroUtil.supernova.SOPRANOS.PDF(data,redshift,Vect0,x(1),x(2),VecMs,VecFrho,VecEbv,model,6,ProgType,priors,bg,transient,mintpoints,c));
    fun = @(x) (1-PDFfun(x));
    x0 = [Rs_in(ipeak), Vs_in(ipeak)];
    PDF0 = PDFfun(x0);
    [x,fval,exitflag,output] = fminsearch(fun,x0);
    if (exitflag == 0) || (1-fval>PDF0)
        converge = 1;
        Rs = x(1);
        Vs = x(2);
        [PDFvalue, chi_2, dof, chi2bands, pntsbands] = PDFfun(x);
        PDFmat = chi2pdf(chi_2,dof);
        PDFmat(dof<=6)=0;
        [~,idx] = max(PDFmat(:));
        [idxt0,idxMs,idxFrho,idxEbv] = ind2sub(size(PDFmat),idx);
        t0 = Vect0(idxt0);
        Ms = VecMs(idxMs);
        Ebv = VecEbv(idxEbv);
        frho= VecFrho(idxFrho);       
    else
        disp('fminsearch did not converge')
        converge = 0;
        [PDFvalue, chi_2, dof, chi2bands, pntsbands] = PDFfun(x0);
        PDFmat = chi2pdf(chi_2,dof);       
        PDFmat(dof<=6)=0;
        [~,idx] = max(PDFmat(:));
        [idxt0,idxMs,idxFrho,idxEbv] = ind2sub(size(PDFmat),idx);
        t0 = Vect0(idxt0);
        Ms = VecMs(idxMs);
        Ebv = VecEbv(idxEbv);
        frho= VecFrho(idxFrho);               
    end
    margt0 = trapz(log10(VecFrho),trapz(VecMs,trapz(VecEbv,PDFmat,4),2),3);
    margMs = trapz(log10(VecFrho),trapz(Vect0,trapz(VecEbv,PDFmat,4),1),3);
    margEbv=trapz(log(VecFrho),trapz(Vect0,trapz(VecMs,PDFmat,2),1),3);
    margFrho=trapz(VecEbv,trapz(Vect0,trapz(VecMs,PDFmat,2),1),4);

    chi2values.Rs       = Rs;
    chi2values.Vs       = Vs;
    chi2values.PDFvalue = PDFvalue;
    chi2values.chi2     = chi_2;
    chi2values.dof      = dof;
    chi2values.chi2bands= chi2bands;
    chi2values.pntsbands = pntsbands;
    chi2values.t0val    = t0;
    chi2values.idxt0    = idxt0;
    chi2values.Msval    = Ms;
    chi2values.idxMs    = idxMs;
    chi2values.Ebvval   = Ebv;
    chi2values.idxEbv   = idxEbv;
    chi2values.frho     = frho;
    chi2values.idxFrho  = idxFrho;
    chi2values.margt0   = margt0;
    chi2values.margMs   = margMs;
    chi2values.margEbv  = margEbv;
    chi2values.margFrho = margFrho;
    chi2values.PDFmat   = PDFmat;

    fprintf('End of step 1, Rs=%4.0f, Vs=%4.2f, Ms=%4.1f, f_rho=%5.3f, t_ref=%8.2f, Ebv=%6.4f chi2/dof=%6.2f/%d\n', Rs, Vs/10^8.5, Ms, frho, t0, Ebv, chi_2(idxt0, idxMs, idxFrho, idxEbv), dof(idxt0, idxMs, idxFrho, idxEbv));

    fun = @(x) 1-(AstroUtil.supernova.SOPRANOS.PDF(data,redshift,x(1),Rs,Vs,x(2),x(3),x(4),model,6,ProgType,priors,bg,transient,mintpoints,c));
    x0 = [t0 Ms frho Ebv];
    f0 = fun(x0);
    [x,fval,exitflag,output] = fminsearch(fun, x0);
    if (exitflag == 1) && (fval>f0)
        t0 = Vect0(idxt0);
        Ms = VecMs(idxMs);
        Ebv= VecEbv(idxEbv);
        frho = VecFrho(idxFrho);
    else
        t0 = x(1);
        Ms = x(2);
        frho = x(3);
        Ebv = x(4);
        [PDFvalue, chi_2, dof, chi2bands, pntsbands] = AstroUtil.supernova.SOPRANOS.PDF(data,redshift,x(1),Rs,Vs,x(2),x(3),x(4),model,6,ProgType,priors,bg,transient,mintpoints,c);

        chi2values.t0MsEbv.PDF  = PDFvalue;
        chi2values.t0MsEbv.chi2 = chi_2;
        chi2values.t0MsEbv.dof = dof;
        chi2values.t0MsEbv.chi2bands = chi2bands;
        chi2values.t0MsEbv.pntsbands = pntsbands;

        chi2values.t0MsEbv.t0       = t0;
        chi2values.t0MsEbv.Ms       = Ms;
        chi2values.t0MsEbv.Ebv      = Ebv;
        chi2values.t0MsEbv.frho     = frho;
    end

    fprintf('End of step 2, Rs=%4.0f, Vs=%4.2f, Ms=%4.1f, f_rho=%5.3f, t_ref=%8.2f, Ebv=%6.4f chi2/dof=%6.2f/%d\n', Rs, Vs/10^8.5, Ms, frho, t0, Ebv, chi_2, dof);

    fun = @(x) 1-(AstroUtil.supernova.SOPRANOS.PDF(data,redshift,x(1),x(2),x(3),x(4),x(5),x(6),model,6,ProgType,priors,bg,transient,mintpoints,c));
    x0 = [t0 Rs Vs Ms frho Ebv];
    f0 = fun(x0);
    [x,fval,exitflag,output] = fminsearch(fun, x0);
    if (exitflag == 1) && (fval>f0)
        t0  = x0(1);
        Rs  = x0(2);
        Vs  = x0(3);
        Ms  = x0(4);
        frho= x0(5);        
        Ebv = x0(6);
    else
        t0  = x(1);
        Rs  = x(2);
        Vs  = x(3);
        Ms  = x(4);
        frho= x(5);        
        Ebv = x(6);
        [PDFval, chi_2, dof, chi2bands, pntsbands] = AstroUtil.supernova.SOPRANOS.PDF(data,redshift,x(1),x(2),x(3),x(4),x(5),x(6),model,6,ProgType,priors,bg,transient,mintpoints,c);

        chi2values.all.PDF  = PDFval;
        chi2values.all.chi2 = chi_2;
        chi2values.all.dof = dof;
        chi2values.all.chi2bands = chi2bands;
        chi2values.all.pntsbands = pntsbands;
        chi2values.all.t0       = x(1);
        chi2values.all.Rs       = x(2);
        chi2values.all.Vs       = x(3);
        chi2values.all.Ms       = x(4);
        chi2values.all.frho     = x(5);
        chi2values.all.Ebv      = x(6);
        chi2values.all.f0       = f0;
        chi2values.all.fval     = fval;
        chi2values.all.exitflag = exitflag;
    end
    chi2valuesArray(ipeak) = chi2values;
    fprintf('End of step 3, Rs=%4.0f, Vs=%4.2f, Ms=%4.1f, f_rho=%5.3f, t_ref=%8.2f, Ebv=%6.4f chi2/dof=%6.2f/%d\n\n', Rs, Vs/10^8.5, Ms, frho, t0, Ebv, chi_2, dof);
end
% end

end