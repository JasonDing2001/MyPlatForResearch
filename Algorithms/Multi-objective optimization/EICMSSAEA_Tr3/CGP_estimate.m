function [PopCon,PopCon_b,CMSE,CMSE_b] = CGP_estimate(X,CModel,numOBJ,cbeta)
            [N,~]  = size(X);
            PopCon = zeros(N,numOBJ);
            CMSE    = zeros(N,numOBJ);
            for i = 1: N
                for j = 1 : numOBJ
                    [PopCon(i,j),~,CMSE(i,j)] = predictor(X(i,:),CModel{j});
                end
            end
            CMSE = max(CMSE,0);
            CS_ = sqrt(CMSE);
            CMSE = CS_;
            
            PopCon_b = PopCon;
            CMSE_b = CMSE;

            PopCon = PopCon+cbeta*CS_;
end

