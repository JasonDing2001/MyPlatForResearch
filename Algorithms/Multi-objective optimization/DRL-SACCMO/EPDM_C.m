function xbest = EPDM_C(TSDec,TSObj,TSCon,PopDec,PopObj,PopCon,Popreal,W,Model,CModel,Problem)
% 约束场景的 EPDM 指标（PoF + CV 修正）

    [PopObj,MSE,PopCon,CMSE] = PredictSurrogate(PopDec,Model,CModel);

    if isempty(Popreal) || size(Popreal,2) ~= size(PopDec,2)
        % Popreal 为空或维度不一致时，视为无已选样本
        Lia = false(size(PopDec,1),1);
    else
        Lia = ismember(PopDec,Popreal,'rows');
    end
    TSDec = [TSDec;PopDec(Lia==1,:)];
    TSObj = [TSObj;PopObj(Lia==1,:)];
    TSCon = [TSCon;PopCon(Lia==1,:)];

    PopDec = PopDec(Lia==0,:);
    PopObj = PopObj(Lia==0,:);
    PopCon = PopCon(Lia==0,:);
    MSE    = MSE(Lia==0,:);
    CMSE   = CMSE(Lia==0,:);

    if isempty(PopDec)
        xbest = [];
        return;
    end

    if isempty(TSObj)
        xbest = PopDec(randi(size(PopDec,1)),:);
        return;
    end

    l  = randperm(size(W,1),1);
    aW = W(l,:);

    CV_TS = sum(max(0,TSCon),2);
    feasible_TS = find(CV_TS==0);
    CV_pred = sum(max(0,PopCon),2);
    feasible_pop = find(CV_pred==0);
    N = size(PopDec,1);
    theta = 5;

    CMSE = max(CMSE,1e-12);

    if ~isempty(feasible_TS)
        refObj = TSObj(feasible_TS,:);
    else
        refObj = TSObj;
    end
    CosineAngle = 1-pdist2(refObj,aW,'cosine');
    SineAngle   = sqrt(1 - CosineAngle .^ 2);
    PDM = theta*sqrt(sum(refObj.^2,2)).*SineAngle+mean(refObj,2);
    minPDM = min(PDM,[],1);

    num_sample = 200; % TODO: 可调整采样数
    inScore = zeros(N,1);
    for i = 1:N
        Sigma = diag(max(MSE(i,:).^2,1e-12));
        rand_samples = mvnrnd(PopObj(i,:),Sigma,num_sample);
        Score = zeros(num_sample,1);
        for num = 1:num_sample
            sample = rand_samples(num,:);
            CosineAngle = 1-pdist2(sample,aW,'cosine');
            SineAngle   = sqrt(1 - CosineAngle .^ 2);
            Score(num,1) = minPDM-(theta*sqrt(sum(sample.^2,2))*SineAngle+mean(sample,2));
        end
        inScore(i,1) = mean(Score);
    end

    Pof = prod(normcdf((0-PopCon)./CMSE),2);
    if ~isempty(feasible_pop)
        cand = feasible_pop;
        Score = Pof.*inScore;
        [~,best_idx] = max(Score(cand));
        best = cand(best_idx);
    else
        minCV = min(CV_pred);
        cand = find(CV_pred == minCV);
        if numel(cand) == 1
            best = cand;
        else
            Score = Pof./(1+CV_pred);
            [~,best_idx] = max(Score(cand));
            best = cand(best_idx);
        end
    end

    xbest = PopDec(best,:);
end
