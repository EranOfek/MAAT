function [x1, x2, x3, x4, epsilon, ValueAtRoot]=QuarticSolverVec(a,b,c,d,e)
% Quartic integral solver (vectorized): [x1, x2, x3, x4]=QuarticSolverVec(a,b,c,d,e)
% Package: Util.integral
% Description: Quartic integral solver (vectorized)
%   v.0.1 - Nearly identical to QuarticSolver v. 0.4, the first successful vectorized implimentation 
%           Changed logic of ChosenSet to accomudate simultaneous convergence of sets 1 & 2
%         - Note the periodicity in nearly-convergent solutions can other
%           than four (related to text on step 4 after table 3). examples:
%           period of 5: [a,b,c,d,e]=[0.111964240308252 -0.88497524334712 -0.197876116344933 -1.07336408259262 -0.373248675102065];
%           period of 6: [a,b,c,d,e]=[-1.380904438798326 0.904866918945240 -0.280749330818231 0.990034312758900 1.413106456228119];
%           period of 22: [a,b,c,d,e]=[0.903755513939902 0.490545114637739 -1.389679906455410 -0.875910689438623 -0.290630547104907];
%           Therefore condition was changed from epsilon1(iter)==0 to epsilon1(iter)<8*eps (and similarl for epsilon2)
%         - Special case criterion of the analytical formula was changed to
%           ind=abs(4*Delta0.^3./Delta1.^2)<2*eps;  (instead of exact zero)
%         - vectorized
%   ============================================
%   - Solves for the x1-x4 roots of the quartic equation y(x)=ax^4+bx^3+cx^2+dx+e.
%     Multiple eqations can be soved simultaneously by entering same-sized column vectors on all inputs.
%   - Note the code immediatly tanslates the input parameters ["a","b","c","d","e"] to the reference paper parameters [1,a,b,c,d] for consistency,
%     and the code probably performes best when "a"=1.
% Input  : - A
%          - B
%          - C
%          - D
%          - E
% Output : * The four free output parameters are the polynomial roots.
%            The function always returns four (possibly complex) values. Multiple roots, if exist, are given multiple times. An error will result in four NaN values.
%            No convergence will result in four inf values (still?)
%          - epsilon    : the sum of the magnitudes of the correction steps on of the last itereation (should be a small number)
%          - ValueAtRoot: the value of the polynimial evaluated at the roots (AFTER division by "a" - note translation of parameters above).
%                Should be a very small number in magnitude (ideally zero).
%
% reference: 
%    Peter Strobach (2010), Journal of Computational and Applied Mathematics 234
%    http://www.sciencedirect.com/science/article/pii/S0377042710002128
% By    : Aviv Ofir                  Oct 2017
% Reliable: 1


MaxIter=16;

% INPUT CONTROL
% all-column vectors only
if size(a,1)~=size(b,1) || size(a,1)~=size(c,1) || size(a,1)~=size(d,1) || size(a,1)~=size(e,1) || ...
   size(a,2)~=1 || size(b,2)~=1 || size(c,2)~=1 || size(d,2)~=1 || size(e,2)~=1
    fprintf('ERROR: illegal input parameter sizes.\n');
    x1=inf; x2=inf; x3=inf; x4=inf;    
    return;
end

% translate input variables to the paper's
if any(a==0)
   fprintf('ERROR: a==0. Not a quartic equation.\n')
    x1=NaN; x2=NaN; x3=NaN; x4=NaN;    
    return;
else
    input_a=a;
    input_b=b;
    input_c=c;
    input_d=d;
    input_e=e;
    a=input_b./input_a;
    b=input_c./input_a;
    c=input_d./input_a;
    d=input_e./input_a;
    clear e input_a input_b input_c input_d;
end

% PRE-ALLOCATE MEMORY
% ChosenSet is used to track which input set already has a solution (=non-zero value)
ChosenSet=zeros(size(a));
[x1, x2, x3, x4, x11, x12, x21, x22, alpha01, alpha02, beta01, beta02, gamma01, gamma02, delta01, delta02, e11, e12, e13, e14, e21, e22, e23, e24, alpha1, alpha2, beta1, beta2, gamma1, gamma2, delta1, delta2, alpha, beta, gamma, delta]=deal(NaN(size(a)));

% check multiple roots -cases 2 & 3. indexed by ChosenSet=-2
test_alpha=0.5*a;
test_beta=0.5*(b-test_alpha.^2);
test_epsilon=[c-2*test_alpha.*test_beta d-test_beta.^2];
ind=all(test_epsilon==0,2);
if any(ind)
    [x1(ind), x2(ind)]=SolveQuadratic(1,test_alpha(ind),test_beta(ind));
    x3(ind)=x1(ind); x4(ind)=x2(ind);
    ChosenSet(ind)=-2;
end;

% check multiple roots -case 4. indexed by ChosenSet=-4
i=ChosenSet==0;
[x11(i), x12(i)]=SolveQuadratic(ones(sum(i),1),a(i)/2,b(i)/6);
x21(i)=-a(i)-3*x11(i);    
test_epsilon(i,1:2)=[c+x11(i).^2.*(x11(i)+3*x21(i)) d(i)-x11(i).^3.*x21(i)];
ind(i)=all(test_epsilon(i)==0,2);
if any(ind(i))
    x1(ind(i))=x11(ind(i)); x2(ind(i))=x11(ind(i)); x3(ind(i))=x11(ind(i)); x4(ind(i))=x12(ind(i));
    ChosenSet(ind(i))=-4;
end;    
x22(i)=-a(i)-3*x12(i);
test_epsilon(i,1:2)=[c+x12(i).^2.*(x12(i)+3*x22(i)) d(i)-x12(i).^3.*x22(i)];
ind(i)=all(test_epsilon(i)==0,2);
if any(ind(i))
    x1(ind(i))=x21(ind(i)); x2(ind(i))=x21(ind(i)); x3(ind(i))=x21(ind(i)); x4(ind(i))=x22(ind(i));
    ChosenSet(ind(i))=-4;
end;    
    
% General solution
% initilize
clear x11 x12 x21 x22 test_alpha test_beta test_epsilon
epsilon1=inf(size(a,1),MaxIter,1);
epsilon2=epsilon1;

i=ChosenSet==0;
fi=find(i);
[x(i,1), x(i,2), x(i,3), x(i,4)]=AnalyticalSolution(ones(sum(i),1),a(i),b(i),c(i),d(i));
[~, ind]=sort(abs(x(i,:)),2,'descend');
x1(i)=x(sub2ind(size(x),fi,ind(:,1)));
x2(i)=x(sub2ind(size(x),fi,ind(:,2)));
x3(i)=x(sub2ind(size(x),fi,ind(:,3)));
x4(i)=x(sub2ind(size(x),fi,ind(:,4)));
clear x ind

alpha01(i)=-real(x1(i)+x2(i));
beta01(i)=real(x1(i).*x2(i));
alpha02(i)=-real(x2(i)+x3(i));
beta02(i)=real(x2(i).*x3(i));
[gamma01(i), delta01(i)]=FastGammaDelta(alpha01(i),beta01(i),a(i),b(i),c(i),d(i));
[gamma02(i), delta02(i)]=FastGammaDelta(alpha02(i),beta02(i),a(i),b(i),c(i),d(i));

alpha1(i)=alpha01(i); beta1(i)=beta01(i); gamma1(i)=gamma01(i); delta1(i)=delta01(i);
alpha2(i)=alpha02(i); beta2(i)=beta02(i); gamma2(i)=gamma02(i); delta2(i)=delta02(i);

%Backward Optimizer Outer Loop
e11(i)=a(i)-alpha1(i)-gamma1(i);
e12(i)=b(i)-beta1(i)-alpha1(i).*gamma1(i)-delta1(i);
e13(i)=c(i)-beta1(i).*gamma1(i)-alpha1(i).*delta1(i);
e14(i)=d(i)-beta1(i).*delta1(i);

e21(i)=a(i)-alpha2(i)-gamma2(i);
e22(i)=b(i)-beta2(i)-alpha2(i).*gamma2(i)-delta2(i);
e23(i)=c(i)-beta2(i).*gamma2(i)-alpha2(i).*delta2(i);
e24(i)=d(i)-beta2(i).*delta2(i);

iter=0;
while iter<MaxIter && any(ChosenSet(i)==0)
    iter=iter+1;
    i=find(ChosenSet==0);
    
    [alpha1(i), beta1(i), gamma1(i), delta1(i), e11(i), e12(i), e13(i), e14(i), epsilon1(i,iter)]=BackwardOptimizer_InnerLoop(a(i),b(i),c(i),d(i),alpha1(i),beta1(i),gamma1(i),delta1(i),e11(i),e12(i),e13(i),e14(i));
    [alpha2(i), beta2(i), gamma2(i), delta2(i), e21(i), e22(i), e23(i), e24(i), epsilon2(i,iter)]=BackwardOptimizer_InnerLoop(a(i),b(i),c(i),d(i),alpha2(i),beta2(i),gamma2(i),delta2(i),e21(i),e22(i),e23(i),e24(i));

    [BestEps, j]=min([epsilon1(i,iter) epsilon2(i,iter)],[],2);
    ind=BestEps<8*eps;
    ChosenSet(i(ind))=j(ind);
    ind=~ind;
    if iter>1 && any(ind)
        LimitCycleReached=[any(bsxfun(@eq,epsilon1(i(ind),max(1,iter-4):max(1,iter-1)),epsilon1(i(ind),iter)),2) any(bsxfun(@eq,epsilon2(i(ind),max(1,iter-4):max(1,iter-1)),epsilon2(i(ind),iter)),2)];
        ChosenSet(LimitCycleReached(:,1) & ~LimitCycleReached(:,2))=1;
        ChosenSet(~LimitCycleReached(:,1) & LimitCycleReached(:,2))=2;
        ChosenSet(LimitCycleReached(:,1) & LimitCycleReached(:,2))=j(ind(LimitCycleReached(:,1)));
    end
end

i=find(ChosenSet==0);
ind=epsilon1(i,end)<epsilon2(i,end);
ChosenSet(i(ind))=1;
ChosenSet(i(~ind))=2;

% Output
i=ChosenSet==1;
alpha(i)=alpha1(i);
beta(i)=beta1(i);
gamma(i)=gamma1(i);
delta(i)=delta1(i);
epsilon1=epsilon1(i,:);
epsilon(i,1)=epsilon1(sub2ind(size(epsilon1),(1:size(epsilon1,1)).',sum(isfinite(epsilon1),2)));

i=ChosenSet==2;
alpha(i)=alpha2(i);
beta(i)=beta2(i);
gamma(i)=gamma2(i);
delta(i)=delta2(i);
epsilon2=epsilon2(i,:);
epsilon(i,1)=epsilon2(sub2ind(size(epsilon2),(1:size(epsilon2,1)).',sum(isfinite(epsilon2),2)));

i=ChosenSet>0;
[x1(i), x2(i)]=SolveQuadratic(ones(sum(i),1),alpha(i),beta(i));
[x3(i), x4(i)]=SolveQuadratic(ones(sum(i),1),gamma(i),delta(i));

if nargout>5
    ValueAtRoot=[x1.^4+a.*x1.^3+b.*x1.^2+c.*x1+d ...
                 x2.^4+a.*x2.^3+b.*x2.^2+c.*x2+d ...
                 x3.^4+a.*x3.^3+b.*x3.^2+c.*x3+d ...
                 x4.^4+a.*x4.^3+b.*x4.^2+c.*x4+d];
end


function [x1, x2, x3, x4]=AnalyticalSolution(a,b,c,d,e)
    % reference: https://en.wikipedia.org/wiki/Quartic_function#General_formula_for_roots
    p=(8*a.*c-3*b.^2)./(8*a.^2);
    q=(b.^3-4*a.*b.*c+8*a.^2.*d)./(8*a.^3);
    
    Delta0=c.^2-3*b.*d+12*a.*e;
    Delta1=2*c.^3 -9*b.*c.*d +27*b.^2.*e +27*a.*d.^2 -72*a.*c.*e;

    Q=(0.5*(Delta1+sqrt(Delta1.^2-4*Delta0.^3))).^(1/3);
    %ind=Delta0==0;
    ind=abs(4*Delta0.^3./Delta1.^2)<2*eps;
    if any(ind)
        ind2=false(size(a));
        Delta(ind)=256*a(ind).^3.*e(ind).^3-192*a(ind).^2.*b(ind).*d(ind).*e(ind).^2-128*a(ind).^2.*c(ind).^2.*e(ind).^2+144*a(ind).^2.*c(ind).*d(ind).^2.*e(ind)-27*a(ind).^2.*d(ind).^4 ...
                  +144*a(ind).*b(ind).^2.*c(ind).*e(ind).^2-6*a(ind).*b(ind).^2.*d(ind).^2.*e(ind)-80*a(ind).*b(ind).*c(ind).^2.*d(ind).*e(ind)+18*a(ind).*b(ind).*c(ind).*d(ind).^3+16*a(ind).*c(ind).^4.*e(ind) ...
                 -4*a(ind).*c(ind).^3.*d(ind).^2-27*b(ind).^4.*e(ind).^2+18*b(ind).^3.*c(ind).*d(ind).*e(ind)-4*b(ind).^3.*d(ind).^3-4*b(ind).^2.*c(ind).^3.*e(ind)+b(ind).^2.*c(ind).^2.*d(ind).^2;
        ind2(ind)=Delta(ind)~=0;
        if any(Delta(ind & ind2))
            Q(ind & ind2)=(Delta1(ind & ind2)).^(1/3);
        end
    end
    S=0.5*sqrt(-2/3*p+1./(3*a).*(Q+Delta0./Q));

    x1=-b./(4*a)-S+0.5*sqrt(-4*S.^2-2*p+q./S);
    x2=-b./(4*a)-S-0.5*sqrt(-4*S.^2-2*p+q./S);
    x3=-b./(4*a)+S+0.5*sqrt(-4*S.^2-2*p-q./S);
    x4=-b./(4*a)+S-0.5*sqrt(-4*S.^2-2*p-q./S);
   
function [gamma0, delta0]=FastGammaDelta(alpha0,beta0,a,b,c,d)
    % Table 3
    phi1=a+alpha0.^2+beta0.^2;
    phi2=alpha0.*(1+beta0);
    c1=a-alpha0+alpha0.*(b-beta0)+beta0.*c;
    c2=b-beta0+alpha0.*c+beta0.*d;
    L1=sqrt(phi1);
    L3=phi2./L1;
    L2=sqrt(phi1-phi2.^2./phi1);
    y1=c1./L1;
    y2=(c2-y1.*L3)./L2;
    delta0=y2./L2;
    gamma0=(y1-delta0.*L3)./L1;

function [alpha, beta, gamma, delta, e1, e2, e3, e4, epsilon]=BackwardOptimizer_InnerLoop(a,b,c,d,alpha,beta,gamma,delta,e1,e2,e3,e4)
    U23=alpha-gamma;
    U33=beta-delta-gamma.*U23;
    L43=-delta.*U23./U33;
    U44=beta-delta-L43.*U23;
   
    x1=e1;
    x2=e2-gamma.*x1;
    x3=e3-delta.*x1-gamma.*x2;
    x4=e4-delta.*x2-L43.*x3;
    y4=x4./U44;
    y3=(x3-U23.*y4)./U33;
    y2=x2-U23.*y3-y4;
    y1=x1-y3;
       
    alpha=alpha+y1;
    beta=beta+y2;
    gamma=gamma+y3;
    delta=delta+y4;
   
    e1=a-alpha-gamma;
    e2=b-beta-alpha.*gamma-delta;
    e3=c-beta.*gamma-alpha.*delta;
    e4=d-beta.*delta;   
    epsilon=abs(e1)+abs(e2)+abs(e3)+abs(e4);

function [x1, x2]=SolveQuadratic(a,b,c)
    % Chapter 5.6 from Numerical Recepies
    i=all(imag([a b c])==0,2);
    q(i,1)=-0.5*(b(i)+sign(b(i)).*sqrt(b(i).^2-4*a(i).*c(i)));
    i=~i;
    if any(i)
        s(i)=real(conj(b(i)).*sqrt(b(i).^2-4*a(i).*c(i)));
        s(i(s(i)<0))=-s(i(s(i)<0));
        q(i)=-0.5*(b(i)+s(i));
    end
    
    x1=q./a;
    x2=c./q;