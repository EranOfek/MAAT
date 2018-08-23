function [x1, x2, x3, x4, epsilon, ValueAtRoot]=QuarticSolver(a,b,c,d,e)
% Quartic integral solver: [x1,x2,x3,x4]=QuarticSolver(a,b,c,d,e)
% Package: Util.integral
% Description: Quartic integral solver.
%   v.0.31 - Added some correction detected from QuarticSolverVec
%          Changed logic of ChosenSet to accomudate simultaneous convergence of sets 1 & 2
%        - quickly tanslates the input parameters [a,b,c,d,e] to the reference paper parameters [1,a,b,c,d] for consistency
%        - corrected typo in alpha01/beta01 definitions
%        - corrected typos in AnalyticalSolution, added special case
%        - Note the periodicity in nearly-convergent solutions can other
%          than four (related to text on step 4 after table 3). examples:
%          period of 5: [a,b,c,d,e]=[0.111964240308252 -0.88497524334712 -0.197876116344933 -1.07336408259262 -0.373248675102065];
%          period of 6: [a,b,c,d,e]=[-1.380904438798326 0.904866918945240 -0.280749330818231 0.990034312758900 1.413106456228119];
%          period of 22: [a,b,c,d,e]=[0.903755513939902 0.490545114637739 -1.389679906455410 -0.875910689438623 -0.290630547104907];
%          Therefore condition was changed from epsilon1(iter)==0 to epsilon1(iter)<8*eps (and similarl for epsilon2)
% 	===========================================
%   - Solves for the x1-x4 roots of the quartic equation y(x)=ax^4+bx^3+cx^2+dx+e.
%   - The function always returns four values. Multiple roots, if exist, are given multiple times. No roots will result in four NaN values
%     No convergence will result in four inf values.
% Input : - A
%         - B
%         - C
%         - D
%         - E
% Output: - X1
%         - X2
%         - X3
%         - X4 
%         - X5
% 
% Reference: 
%    PeterStrobach (2010), Journal of Computational and Applied Mathematics 234
%    http://www.sciencedirect.com/science/article/pii/S0377042710002128
% By    : Aviv Ofir                  Oct 2017
% Reliable: 1


MaxIter=16;

% translate input variables to the paper's
if a==0
   fprintf('ERROR: a==0. Not a quartic equation.\n')
    x1=NaN; x2=NaN; x3=NaN; x4=NaN;    
    return;
else
    input_a=a;
    input_b=b;
    input_c=c;
    input_d=d;
    input_e=e;
    a=input_b/input_a;
    b=input_c/input_a;
    c=input_d/input_a;
    d=input_e/input_a;
    clear e input_a input_b input_c input_d;
end

% check multiple roots -cases 2 & 3
test_alpha=0.5*a;
test_beta=0.5*(b-test_alpha^2);
test_epsilon=[c-2*test_alpha*test_beta d-test_beta^2];
if all(test_epsilon==0)
    [x1, x2]=SolveQuadratic(1,test_alpha,test_beta);
    x3=x1; x4=x2;
    return;
end;

% check multiple roots -case 4
[x11, x12]=SolveQuadratic(1,a/2,b/6);
x21=-a-3*x11;    
test_epsilon=[c+x11.^2*(x11+3*x21) d-x11.^3*x21];
if all(test_epsilon==0)
    x1=x11; x2=x11; x3=x11; x4=x12;
    return;
end;    
x22=-a-3*x12;
test_epsilon=[c+x12.^2*(x12+3*x22) d-x12.^3*x22];
if all(test_epsilon==0)
    x1=x21; x2=x21; x3=x21; x4=x22;
    return;
end;    
    
% general solution    
% initilize
clear x11 x12 x21 x22 test_alpha test_beta test_epsilon
epsilon1=inf(MaxIter,1);
epsilon2=epsilon1;

[x(1), x(2), x(3), x(4)]=AnalyticalSolution(1,a,b,c,d);
[~, ind]=sort(abs(x),'descend');
x1=x(ind(1));
x2=x(ind(2));
x3=x(ind(3));
x4=x(ind(4));
clear x ind

alpha01=-real(x1+x2);
beta01=real(x1*x2);
alpha02=-real(x2+x3);
beta02=real(x2*x3);
[gamma01, delta01]=FastGammaDelta(alpha01,beta01,a,b,c,d);
[gamma02, delta02]=FastGammaDelta(alpha02,beta02,a,b,c,d);

alpha1=alpha01; beta1=beta01; gamma1=gamma01; delta1=delta01;
alpha2=alpha02; beta2=beta02; gamma2=gamma02; delta2=delta02;

%Backward Optimizer Outer Loop
e11=a-alpha1-gamma1;
e12=b-beta1-alpha1*gamma1-delta1;
e13=c-beta1*gamma1-alpha1*delta1;
e14=d-beta1*delta1;

e21=a-alpha2-gamma2;
e22=b-beta2-alpha2*gamma2-delta2;
e23=c-beta2*gamma2-alpha2*delta2;
e24=d-beta2*delta2;

iter=0;
ChosenSet=0;
while iter<MaxIter && ChosenSet==0
    iter=iter+1;
    [alpha1, beta1, gamma1, delta1, e11, e12, e13, e14, epsilon1(iter)]=BackwardOptimizer_InnerLoop(a,b,c,d,alpha1,beta1,gamma1,delta1,e11,e12,e13,e14);
    [alpha2, beta2, gamma2, delta2, e21, e22, e23, e24, epsilon2(iter)]=BackwardOptimizer_InnerLoop(a,b,c,d,alpha2,beta2,gamma2,delta2,e21,e22,e23,e24);

    [BestEps, j]=min([epsilon1(iter) epsilon2(iter)]);
    if BestEps<8*eps
        ChosenSet=j;
    elseif iter>1
        LimitCycleReached=[any(epsilon1(iter)==epsilon1(max(1,iter-4):max(1,iter-1)))  any(epsilon2(iter)==epsilon2(max(1,iter-4):max(1,iter-1)))];
        ChosenSet(LimitCycleReached(1) && ~LimitCycleReached(2))=1;
        ChosenSet(~LimitCycleReached(1) && LimitCycleReached(2))=2;
        ChosenSet(LimitCycleReached(1) && LimitCycleReached(2))=j;
    end
end

if ChosenSet==0
   if epsilon1(iter)<epsilon2(iter)
       ChosenSet=1;
   else
       ChosenSet=2;
   end
end
if ChosenSet==1
    alpha=alpha1;
    beta=beta1;
    gamma=gamma1;
    delta=delta1;
    epsilon=epsilon1(iter);
elseif ChosenSet==2
    alpha=alpha2;
    beta=beta2;
    gamma=gamma2;
    delta=delta2;
    epsilon=epsilon2(iter);
end

[x1, x2]=SolveQuadratic(1,alpha,beta);
[x3, x4]=SolveQuadratic(1,gamma,delta);

if nargout>5
    ValueAtRoot=[x1^4+a*x1^3+b*x1^2+c*x1+d ...
                 x2^4+a*x2^3+b*x2^2+c*x2+d ...
                 x3^4+a*x3^3+b*x3^2+c*x3+d ...
                 x4^4+a*x4^3+b*x4^2+c*x4+d];
end



function [x1, x2, x3, x4]=AnalyticalSolution(a,b,c,d,e)
    % reference: https://en.wikipedia.org/wiki/Quartic_function#General_formula_for_roots
    p=(8*a*c-3*b^2)/(8*a^2);
    q=(b^3-4*a*b*c+8*a^2*d)/(8*a^3);
    
    Delta0=c^2-3*b*d+12*a*e;
    Delta1=2*c^3 -9*b*c*d +27*b^2*e +27*a*d^2 -72*a*c*e;

    Q=(0.5*(Delta1+sqrt(Delta1^2-4*Delta0^3)))^(1/3);
    %if Delta0==0
    if abs(4*Delta0.^3./Delta1.^2)<2*eps;
        Delta=256*a^3*e^3-192*a^2*b*d*e^2-128*a^2*c^2*e^2+144*a^2*c*d^2*e-27*a^2*d^4 ...
             +144*a*b^2*c*e^2-6*a*b^2*d^2*e-80*a*b*c^2*d*e+18*a*b*c*d^3+16*a*c^4*e ...
             -4*a*c^3*d^2-27*b^4*e^2+18*b^3*c*d*e-4*b^3*d^3-4*b^2*c^3*e+b^2*c^2*d^2;
        if Delta~=0
            Q=(Delta1)^(1/3);
        end
    end
    S=0.5*sqrt(-2/3*p+1/(3*a)*(Q+Delta0/Q));

    x1=-b/(4*a)-S+0.5*sqrt(-4*S^2-2*p+q/S);
    x2=-b/(4*a)-S-0.5*sqrt(-4*S^2-2*p+q/S);
    x3=-b/(4*a)+S+0.5*sqrt(-4*S^2-2*p-q/S);
    x4=-b/(4*a)+S-0.5*sqrt(-4*S^2-2*p-q/S);
   
function [gamma0, delta0]=FastGammaDelta(alpha0,beta0,a,b,c,d)
    % Table 3
    phi1=a+alpha0^2+beta0^2;
    phi2=alpha0*(1+beta0);
    c1=a-alpha0+alpha0*(b-beta0)+beta0*c;
    c2=b-beta0+alpha0*c+beta0*d;
    L1=sqrt(phi1);
    L3=phi2/L1;
    L2=sqrt(phi1-phi2^2/phi1);
    y1=c1/L1;
    y2=(c2-y1*L3)/L2;
    delta0=y2/L2;
    gamma0=(y1-delta0*L3)/L1;

function [alpha, beta, gamma, delta, e1, e2, e3, e4, epsilon]=BackwardOptimizer_InnerLoop(a,b,c,d,alpha,beta,gamma,delta,e1,e2,e3,e4)
    U23=alpha-gamma;
    U33=beta-delta-gamma*U23;
    L43=-delta*U23/U33;
    U44=beta-delta-L43*U23;
   
    x1=e1;
    x2=e2-gamma*x1;
    x3=e3-delta*x1-gamma*x2;
    x4=e4-delta*x2-L43*x3;
    y4=x4/U44;
    y3=(x3-U23*y4)/U33;
    y2=x2-U23*y3-y4;
    y1=x1-y3;
       
    alpha=alpha+y1;
    beta=beta+y2;
    gamma=gamma+y3;
    delta=delta+y4;
   
    e1=a-alpha-gamma;
    e2=b-beta-alpha*gamma-delta;
    e3=c-beta*gamma-alpha*delta;
    e4=d-beta*delta;   
    epsilon=abs(e1)+abs(e2)+abs(e3)+abs(e4);

function [x1, x2]=SolveQuadratic(a,b,c)
    % Chapter 5.6 from NumericalRecepies

    if all(isreal([a b c]))
        q=-0.5*(b+sign(b)*sqrt(b^2-4*a*c));
    else
        s=real(conj(b)*sqrt(b^2-4*a*c));
        if s<0
            s=-s;
        end
        q=-0.5*(b+s);        
    end
    x1=q/a;
    x2=c/q;