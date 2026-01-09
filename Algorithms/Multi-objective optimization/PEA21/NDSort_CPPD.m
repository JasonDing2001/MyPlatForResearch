function [FrontNo,MaxFNo] = NDSort_CPPD(PopObj,ObjMSE,PopCon,ConMSE,nSort)
% Constrained Probabilistic Dominance (CPD)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    N         = size(PopObj,1);
    LFP       = Feasible_Probability(PopCon,ConMSE);
    sigma     = sqrt(ObjMSE(reshape(ones(N,1)*(1:N),N*N,1),:) + repmat(ObjMSE,N,1));
    mean      = PopObj(reshape(ones(N,1)*(1:N),N*N,1),:) - repmat(PopObj,N,1);
    x_PD      = normcdf((0-mean)./sigma);
    y_PD      = 1 - x_PD;
    x_PD      = - x_PD.*LFP(reshape(ones(N,1)*(1:N),N*N,1),:);
    y_PD      = - y_PD.*repmat(LFP,N,1);
    
    [x_CPD, y_CPD]  =  PoCB(ConMSE,PopCon);
       
    dominate  = false(N);
    
    
    for i = 1 : N - 1
        for j = i + 1 : N
            % 检查第 i 个体是否支配第 j 个体
            if all(x_CPD(N * (i - 1) + j, :) <= y_CPD(N * (i - 1) + j, :)) && ~all(x_CPD(N * (i - 1) + j, :) == y_CPD(N * (i - 1) + j, :))
                dominate(i, j) = true;  % i 支配 j
            % 检查第 j 个体是否支配第 i 个体
            elseif all(x_CPD(N * (i - 1) + j, :) >= y_CPD(N * (i - 1) + j, :)) && ~all(x_CPD(N * (i - 1) + j, :) == y_CPD(N * (i - 1) + j, :))
                dominate(j, i) = true;  % j 支配 i
            end
        end
    end
    
    
    for i = 1 : N-1
        for j = i+1 : N
            if all(x_PD(N*(i-1)+j,:) <= y_PD(N*(i-1)+j,:)) && ~all(x_PD(N*(i-1)+j,:) == y_PD(N*(i-1)+j,:))
                dominate(i,j) = true;
            elseif all(x_PD(N*(i-1)+j,:) >= y_PD(N*(i-1)+j,:)) && ~all(x_PD(N*(i-1)+j,:) == y_PD(N*(i-1)+j,:))
                dominate(j,i) = true;
            end
        end
    end

    FrontNo = inf(1,N);
    MaxFNo  = 0;
    while sum(FrontNo~=inf) < min(nSort,N)
        MaxFNo                     = MaxFNo + 1;
        current                    = find(FrontNo==inf);
        dominate_                  = sum(dominate(current,current),1);
        index                      = find(dominate_==min(dominate_));
        FrontNo(current(index))    = MaxFNo;
        dominate(current(index),:) = false;
    end
end

function LFP = Feasible_Probability(PopCon,ConMSE)
    [N,M] = size(PopCon);
    LFP   = ones(N,1);
    for i = 1 : N
        for j = 1 : M
             LFP(i) = min([LFP(i),normcdf((0-(PopCon(i,j)+0*sqrt(ConMSE(i,j))))/sqrt(ConMSE(i,j)))]);
        end
    end
end

function [x_CPD, y_CPD] = PoCB(ConMSE,PopCon)
    N = size(PopCon, 1);
    CV = sum(max(0,PopCon),2);
    MSE = sum(ConMSE.^2,2);
    mean  = CV(reshape(ones(N, 1) * (1:N), N * N, 1), :) - repmat(CV, N, 1);  % 约束之间的差值
    sigma = sqrt(MSE(reshape(ones(N, 1) * (1:N), N * N, 1), :) + repmat(MSE, N, 1));
    x_CPD  = normcdf((0-mean) ./ sigma);  % 支配概率 x_PD
    y_CPD  = 1 - x_CPD;  % 非支配概率 y_PD
    x_CPD = - x_CPD;
    y_CPD = - y_CPD;
end    