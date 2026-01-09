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

    % 将 CV 归一化并拼接到目标，用于约束子任务的角度/PDM计算
    [TSObjC,PopObjC] = AppendNormalizedCV(TSObj,TSCon,PopObj,PopCon);
    if isempty(TSObjC)
        xbest = PopDec(randi(size(PopDec,1)),:);
        return;
    end

    l  = randperm(size(W,1),1);
    aW = W(l,:);

    CV = sum(max(0,TSCon),2);
    feasible = find(CV==0);
    N = size(PopDec,1);
    theta = 5;

    CMSE = max(CMSE,1e-12);

    if ~isempty(feasible)
        CosineAngle = 1-pdist2(TSObjC,aW,'cosine');
        SineAngle   = sqrt(1 - CosineAngle .^ 2);
        PDM = theta*sqrt(sum(TSObjC.^2,2)).*SineAngle+mean(TSObjC,2);
        f_PDM = PDM(feasible);
        minPDM = min(f_PDM,[],1);

        num_sample = 200; % TODO: 可调整采样数
        inScore = zeros(N,1);
        CV_norm = PopObjC(:,end);
        for i = 1:N
            Sigma = diag(max(MSE(i,:).^2,1e-12));
            rand_samples = mvnrnd(PopObj(i,:),Sigma,num_sample);
            Score = zeros(num_sample,1);
            for num = 1:num_sample
                sample_aug = [rand_samples(num,:), CV_norm(i)];
                CosineAngle = 1-pdist2(sample_aug,aW,'cosine');
                SineAngle   = sqrt(1 - CosineAngle .^ 2);
                Score(num,1) = minPDM-(theta*sqrt(sum(sample_aug.^2,2))*SineAngle+mean(sample_aug,2));
            end
            inScore(i,1) = mean(Score);
        end

        Pof = prod(normcdf((0-PopCon)./CMSE),2);
        CV_pred = sum(max(0,PopCon),2);
        CEPDM = Pof.*inScore./(1+CV_pred);
        [~,best] = max(CEPDM);
    else
        Pof = prod(normcdf((0-PopCon)./CMSE),2);
        CV_pred = sum(max(0,PopCon),2);
        Score = Pof./(1+CV_pred);
        [~,best] = max(Score);
    end

    xbest = PopDec(best,:);
end
