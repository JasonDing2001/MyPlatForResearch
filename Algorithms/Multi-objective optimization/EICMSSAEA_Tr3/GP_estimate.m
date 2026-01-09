function [PopObj,PopObj_b,MSE,MSE_b] = GP_estimate(X,Model,numOBJ,beta)
            [N,~]  = size(X);
            PopObj = zeros(N,numOBJ);
            MSE    = zeros(N,numOBJ);
            for i = 1: N
                for j = 1 : numOBJ
                    [PopObj(i,j),~,MSE(i,j)] = predictor(X(i,:),Model{j});
                end
            end
            MSE = max(MSE,0);
            S_ = sqrt(MSE);
            
            %MSE = S_.*(MSE<=1)+MSE.*(MSE>1);% AUCB
            MSE = S_;% AUCB
            PopObj_b = PopObj;
            MSE_b = MSE;
            PopObj = PopObj+beta*MSE;% AUCB
end

