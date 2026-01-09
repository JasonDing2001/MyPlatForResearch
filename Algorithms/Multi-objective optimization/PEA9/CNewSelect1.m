function PopNew = CNewSelect1(PopDec,PopObj,PopCon,ObjMSE,ConMSE,A1,A2,mu,Problem,PopMSE)
    %% Preparing Data
    index   = ismember(PopDec,A1.decs,'rows');
    PopNew  = [];
    if sum(index) == size(PopDec,1)
        return;
    elseif sum(~index) <= mu
        PopNew_ = PopDec(~index,:);
        for i = 1:size(PopNew_,1)
            dist2 = pdist2(real(PopNew_(i,:)),real(A1.decs));
            if min(dist2) > 1e-5
                PopNew = [PopNew;PopNew_(i,:)];
            end
        end
        return;
    end
  
    
    % Normalize the objectives
    PopDec = PopDec(~index,:);
    PopObj = PopObj(~index,:); ObjMSE  = ObjMSE(~index,:);
    PopCon = PopCon(~index,:); ConMSE  = ConMSE(~index,:);
    A1Obj  = A1.objs;
    A1Con  = A1.cons;
    A2Obj = A2.objs;
    zmin   = min([PopObj],[],1); zmax = max([PopObj],[],1);
    A1Obj  = (A1Obj - zmin)./max(zmax - zmin,10e-10);
    PopObj = (PopObj - zmin)./max(zmax - zmin,10e-10);
    
    %% Perform K-Means Clustering on PopObj
    % Assuming you want to cluster PopObj into K clusters
    K = min(mu, size(PopObj, 1));  % Choose K based on the available individuals or set a fixed value
    [Idx, Centers] = kmeans(PopObj, K);
    

    LB = -0.5;
    UB = 1.2;
    nSample = 10000;
    randp = 0.3;
    %  Sample and calculate hypervolume improvement
    Lowb = LB .*ones(1,Problem.M);
    Upb  = UB .*ones(1,Problem.M);
    
    S    = UniformPoint(nSample,Problem.M,'Latin');
    S    =S.*repmat(Upb-Lowb,nSample,1)+repmat(Lowb,nSample,1);
    
    S_S = zeros(nSample,nSample);
    for i = 1 : nSample
        y        = sum(repmat(S(i,:),nSample,1)-S<=0,2) == Problem.M;  
        S_S(i,y) = 1 ;
    end 

    RealFirstObj = A2.best.objs;
    if isempty(RealFirstObj)
        [RealFrontNo, ~] = NDSort(A2.objs, 1); 
        Objs = A2.objs;
        RealFirstObj = Objs(RealFrontNo == 1,:); 
    end
    % Efficent EHVI calcualtion using importance sampling       
    EHVI=CalEHVI(RealFirstObj,PopObj,PopMSE,S,S_S);
    
    
    for k = 1 : K
        tempDec = PopDec(Idx==k,:);  % 当前簇的决策变量
        tempEHVI = EHVI(Idx==k,:);   % 当前簇的 EHVI 值
        tempCon = PopCon(Idx==k,:);  % 当前簇的约束值
        LFP = Feasible_Probability(tempCon,ConMSE(Idx==k,:));
        tempCEHVI = tempEHVI .* LFP;
        [~,index] = max(tempCEHVI);
        PopNew = [PopNew;tempDec(index,:)];
            
        
    end
    
end

function LFP = Feasible_Probability(PopCon,ConMSE)
    [N,M] = size(PopCon);
    LFP   = ones(N,1);
    for i = 1 : N
        for j = 1 : M
             LFP(i) = LFP(i) * normcdf((0-(PopCon(i,j)+0*sqrt(ConMSE(i,j))))/sqrt(ConMSE(i,j)));
        end
    end
end