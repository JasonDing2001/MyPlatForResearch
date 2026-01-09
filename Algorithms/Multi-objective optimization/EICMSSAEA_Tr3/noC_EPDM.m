function xbest = noC_EPDM(TSDec,TSObj,PopDec,PopObj,Popreal,W,Problem,Model)
        beta = 0;
        [PopObj,~,MSE,~] = GP_estimate(PopDec,Model,Problem.M,beta);
        %PopObj = PopObj+2.*MSE;
        Lia = ismember(PopDec,Popreal,'rows');
        TSDec = [TSDec;PopDec(Lia==1,:)];
        TSObj = [TSObj;PopObj(Lia==1,:)];
        
        PopDec = PopDec(Lia==0,:);
        PopObj = PopObj(Lia==0,:);
        
        MSE    = MSE(Lia==0,:);
        
        l = randperm(size(W,1),1);
        aW  = W(l,:);
       
        N = size(PopDec,1);
        theta = 5;
        % EPDM
        CosineAngle = 1-pdist2(TSObj,aW,'cosine');
        SineAngle   = sqrt(1 - CosineAngle .^ 2);
        PDM = theta*sqrt(sum(TSObj.^2,2)).*SineAngle+mean(TSObj,2);
        minPDM = min(PDM,[],1);
        num_sample = 1000;
        Score = zeros(num_sample,1);
        inScore = zeros(N,1);
        for i = 1:N
            rand_samples = mvnrnd(PopObj(i,:),diag(MSE(i,:).^2),num_sample);
            for num = 1:num_sample
                CosineAngle = 1-pdist2(rand_samples(num,:),aW,'cosine');
                SineAngle   = sqrt(1 - CosineAngle .^ 2);
                Score(num,1) = minPDM-(theta*sqrt(sum(rand_samples(num,:).^2,2))*SineAngle+mean(rand_samples(num,:),2));
            end
            inScore(i,1) = sum(Score)/num_sample;
        end
        [~,best] = max(inScore);
        xbest = PopDec(best,:);

end

