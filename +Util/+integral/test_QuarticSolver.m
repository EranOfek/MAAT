% generate random coefficients for the vectorized tests
%v=3*(rand(1e6,5)-0.5);
v=3*(rand(1e3,4)-0.5); v=[ones(size(v,1),1) v];
%v=tmp_v(428622,:);

% r=[1 1e6 1e12 1e18];
% M1=[1    0   0   0; ...
%    r(1)  1   0   0; ...
%     0   r(1) 1   0; ...
%     0    0  r(1) 1; ...
%     0    0   0  r(1)];
% M2=[1    0   0  ; ...
%    r(2)  1   0  ; ...
%     0   r(2) 1  ; ...
%     0    0  r(2)];
% M3=[1    0  ; ...
%    r(3)  1  ; ...
%     0   r(3)];
% M4=[ 1  ; ...
%     r(4)];
% v=(M1*M2*M3*M4).';

import Util.integral.*

tic;
[x1, x2, x3, x4, epsilon, ValueAtRoot]=QuarticSolverVec(v(:,1),v(:,2),v(:,3),v(:,4),v(:,5));
toc
%%
PlotArraySize=min(3,ceil(sqrt(size(v,1))));
disp_x=-40:0.1:40;
v_ind=0;
for i=1:PlotArraySize
    for j=1:PlotArraySize
        v_ind=v_ind+1;
        if v_ind<=size(v,1)
            subplot(PlotArraySize,PlotArraySize,v_ind);
            cla;
            disp_y=v(v_ind,1)*disp_x.^4+v(v_ind,2)*disp_x.^3+v(v_ind,3)*disp_x.^2+v(v_ind,4)*disp_x+v(v_ind,5);
            plot(disp_x,disp_y)
            title(sprintf('%g: y=%gx^4%+gx^3%+gx^2%+gx%-g',v_ind,v(v_ind,1),v(v_ind,2),v(v_ind,3),v(v_ind,4),v(v_ind,5)))
            line([min(disp_x) max(disp_x)],[0 0],'color','k');

            Curr_x=[x1(v_ind), x2(v_ind), x3(v_ind), x4(v_ind)];
            Curr_x(abs(imag(Curr_x))>0)=[];
            hold on;
            plot(Curr_x,zeros(1,numel(Curr_x)),'*r');
            hold off;
            %axis([min(disp_x) max(disp_x) -.5 .5]);
            ax=axis;
            text(0.8*ax(1),0.8*ax(4),sprintf('%g unique real roots',numel(unique(Curr_x))))
            if numel(unique(Curr_x))>0
                text(0.2*ax(1),0.8*ax(4),sprintf([': ' repmat('%g ',1,numel(unique(Curr_x)))],unique(Curr_x)))
                text(0.8*ax(1),0.7*ax(4),sprintf('The ValueAtRoot are: %g, %g, %g, %g',ValueAtRoot(v_ind,:)));
            end
        end;
    end
end

%%

% generate random coefficients for the vectorized tests
v=3*(rand(1e5,4)-0.5); v=[ones(size(v,1),1) v];
% cases where vectorized code is different from ITSELF when the v below is part of a matrix or alone (but by one eps(1) )
% v=[1 -0.606821751702061 -1.17991954474538 -0.393589311631564 1.48805682757755]
% v=[1 -0.913233077249547 -0.903690280017304 -1.24269716935743 0.852707574286925]
%v=tmp_v;

tic;
[v1, v2, v3, v4, epsilon, ValueAtRoot]=QuarticSolverVec(v(:,1),v(:,2),v(:,3),v(:,4),v(:,5));
toc
[x1, x2, x3, x4]=deal(zeros(size(v1)));
tic;
for i=1:size(v,1)
    [x1(i,1), x2(i,1), x3(i,1), x4(i,1), epsilon(i,1), ValueAtRoot(i,:)]=QuarticSolver(v(i,1),v(i,2),v(i,3),v(i,4),v(i,5));
end
toc
sum([x1 x2 x3 x4]-[v1 v2 v3 v4])
