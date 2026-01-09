function xbest = EPDM_U(TSDec,TSObj,PopDec,PopObj,Popreal,W,Model,Problem)
% 不考虑约束的 EPDM 指标

    [PopObj, MSE, ~, ~] = PredictSurrogate(PopDec,Model,{});

    if isempty(Popreal) || size(Popreal,2) ~= size(PopDec,2)
        % Popreal 为空或维度不一致时，视为无已选样本
        Lia = false(size(PopDec,1),1);
    else
        Lia = ismember(PopDec,Popreal,'rows');
    end
    TSDec = [TSDec;PopDec(Lia==1,:)];
    TSObj = [TSObj;PopObj(Lia==1,:)];

    PopDec = PopDec(Lia==0,:);
    PopObj = PopObj(Lia==0,:);
    MSE    = MSE(Lia==0,:);

    if isempty(PopDec)
        xbest = [];
        return;
    end

    l  = randperm(size(W,1),1);
    aW = W(l,:);

    if isempty(TSObj)
        xbest = PopDec(randi(size(PopDec,1)),:);
        return;
    end

    N = size(PopDec,1);
    theta = 5;
    CosineAngle = 1-pdist2(TSObj,aW,'cosine');
    SineAngle   = sqrt(1 - CosineAngle .^ 2);
    PDM = theta*sqrt(sum(TSObj.^2,2)).*SineAngle+mean(TSObj,2);
    minPDM = min(PDM,[],1);

    num_sample = 200; % TODO: 可调整采样数平衡精度与开销
    inScore = zeros(N,1);
    for i = 1:N
        Sigma = diag(max(MSE(i,:).^2,1e-12));
        rand_samples = mvnrnd(PopObj(i,:),Sigma,num_sample);
        Score = zeros(num_sample,1);
        for num = 1:num_sample
            CosineAngle = 1-pdist2(rand_samples(num,:),aW,'cosine');
            SineAngle   = sqrt(1 - CosineAngle .^ 2);
            Score(num,1) = minPDM-(theta*sqrt(sum(rand_samples(num,:).^2,2))*SineAngle+mean(rand_samples(num,:),2));
        end
        inScore(i,1) = mean(Score);
    end
    [~,best] = max(inScore);
    xbest = PopDec(best,:);
end
