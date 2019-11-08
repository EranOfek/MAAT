function [OutX,OutY] = apply_tranclass_INV(TranC,X,Y,varargin)
% Apply inverse transformation to X,Y coordinates
% Package: @TranClass
% Description: Apply inverse transformation to X,Y coordinates
% Input  : - A TranClass object
%          - A vector of X coordinates
%          - A vector of Y coordinates
%          - Accuracy for convergence. Default is 1e-3.
% Output : - A vector of inverse transformed X coordinates
%          - A vector of inverse transformed Y coordinates
DefV.Thresh = 1e-8;
DefV.SolverType = 'fsolve';   % 'fsolve' | 'olditer'
InPar = InArg.populate_keyval(DefV,varargin,mfilename);
Thresh=InPar.Thresh;
IC = 1;
%
GuessX = X;
GuessY = Y;
% X par
Idim = 1;
%ParX = cell2mat(TranC(IC).Par{Idim}).';
ParX = par2vector(TranC(IC),Idim).';
% Y par
Idim = 2;
%ParY = cell2mat(TranC(IC).Par{Idim}).';
ParY = par2vector(TranC(IC),Idim).';
switch lower(InPar.SolverType)
    case 'fsolve'
        %       use matlab fsolve to get root, now support vector use
        options=optimoptions(@fsolve,'Display','off','FunctionTolerance',Thresh);%
        X=reshape(X,[length(X),1]);
        Y=reshape(Y,[length(Y),1]);
        [tmp,~]=fsolve(@(x) inv_root(TranC(IC),[X,Y],x),[X,Y],options);
        GuessX=reshape(tmp(:,1),size(GuessX));
        GuessY=reshape(tmp(:,2),size(GuessY));
      %  for k=1:length(X)
      %    
      %      [tmp,~]=fsolve(@(x) inv_root(TranC(IC),[X(k),Y(k)],x),[X(k),Y(k)],options);
      %      GuessX(k)=tmp(1);GuessY(k)=tmp(2);
      %  end
    case 'olditer'
        NotConverge = true;
        Counter = 0;
        while NotConverge
            Counter = Counter + 1;
            
            H = design_matrix(TranC(IC),'X',GuessX,'Y',GuessY);
            
            Idim = 1;
            X1 = H{Idim}*ParX;
            Idim = 2;
            Y1 = H{Idim}*ParY;
            
            DX = X1 - X;
            DY = Y1 - Y;
            
            %[DX, X, X1, GuessX]
            
            if (max(abs(DX))<Thresh && max(abs(DY))<Thresh)
                % converged
                NotConverge = false;
            end
            
            GuessX = GuessX - DX;
            GuessY = GuessY - DY;
            
        end
end

OutX = GuessX;
OutY = GuessY;
    function F = inv_root(TranC,x0,x)
        [y1,y2]=apply_tranclass(TranC,x(:,1),x(:,2));
        F=[x0(:,1)-y1,x0(:,2)-y2];
      %  [y1,y2]=apply_tranclass(TranC,x(1),x(2));
      %  F(1)=x0(1)-y1;
      %  F(2)=x0(2)-y2;
     
    end
end