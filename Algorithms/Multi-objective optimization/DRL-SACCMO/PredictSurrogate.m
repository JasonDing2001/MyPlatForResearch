function [PopObj,MSE,PopCon,CMSE] = PredictSurrogate(PopDec,Model,CModel)
% 使用 Kriging 预测均值与方差（方差取平方根作为标准差）

    N   = size(PopDec,1);
    M   = numel(Model);
    PopObj = zeros(N,M);
    MSE    = zeros(N,M);

    for i = 1 : N
        for j = 1 : M
            [PopObj(i,j),~,MSE(i,j)] = predictor(PopDec(i,:),Model{j});
        end
    end
    MSE = sqrt(max(MSE,0));

    numC  = numel(CModel);
    PopCon = zeros(N,numC);
    CMSE   = zeros(N,numC);
    if numC > 0
        for i = 1 : N
            for j = 1 : numC
                [PopCon(i,j),~,CMSE(i,j)] = predictor(PopDec(i,:),CModel{j});
            end
        end
        CMSE = sqrt(max(CMSE,0));
    end
end
