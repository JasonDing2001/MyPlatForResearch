function Popreal_Dec2 = SelectByDeltaPBI(ND_TSObj,ND_TSDec,ND_PopObj,ND_PopDec,W)
% PBI 改善最大（参考向量驱动收敛/分布）

    if isempty(ND_PopDec)
        Popreal_Dec2 = [];
        return;
    end

    AllPopObj = [ND_TSObj;ND_PopObj];
    AllPopDec = [ND_TSDec;ND_PopDec];
    [z,znad]  = deal(min(AllPopObj),max(AllPopObj));
    [AllPopObj,~,~] = NormalizationCC(AllPopObj,z,znad);

    N1 = size(ND_TSObj,1);
    ND_TSObj = AllPopObj(1:N1,:);
    ND_TSDec = AllPopDec(1:N1,:);
    normT  = sqrt(sum(ND_TSObj.^2,2));
    Cosine = 1 - pdist2(ND_TSObj,W,'cosine');
    d1     = repmat(normT,1,size(W,1)).*Cosine;
    d2     = repmat(normT,1,size(W,1)).*sqrt(1-Cosine.^2);
    [~,class_T] = min(d2,[],2);

    theta = 5;
    [~,ia,~] = unique(class_T);
    unique_class_T = class_T(ia);
    Archive_PBI = -1*10000*ones(size(W,1),1);
    for i = unique_class_T'
        current = find(class_T==i);
        PBI = d1(current,i)+theta*d2(current,i);
        [bestPBI,~] = min(PBI);
        Archive_PBI(i,1) = bestPBI;
    end

    activeW = W(unique_class_T,:);
    ND_PopObj = AllPopObj(N1+1:end,:);
    ND_PopDec = AllPopDec(N1+1:end,:);
    normP  = sqrt(sum(ND_PopObj.^2,2));
    Cosine = 1 - pdist2(ND_PopObj,activeW,'cosine');
    Pd1     = repmat(normP,1,size(activeW,1)).*Cosine;
    Pd2     = repmat(normP,1,size(activeW,1)).*sqrt(1-Cosine.^2);
    [~,class_P] = min(Pd2,[],2);

    [~,ia,~] = unique(class_P);
    unique_class_C = class_P(ia);
    Candidate_PBI = 10000*ones(size(activeW,1),1);
    Archive_PBI_new = -1*10000*ones(size(W,1),1);
    Next = zeros(1,size(activeW,1));
    class_CC = [];
    for i = unique_class_C'
        current = find(class_P==i);
        PBI = Pd1(current,i)+theta*Pd2(current,i);
        [bestPBI,bestPBI_index] = min(PBI);
        Candidate_PBI(i,1) = bestPBI;
        class_C = unique_class_T(i);
        Archive_PBI_new(class_C,1) = Candidate_PBI(i,1);
        Next(i)  = current(bestPBI_index);
        class_CC = cat(1,class_CC,class_C);
    end
    index = Next(Next~=0);
    PBI_old = Archive_PBI(class_CC,:);
    PBI_new = Archive_PBI_new(class_CC,:);
    delta = PBI_new-PBI_old;
    [~,best] = sort(delta);
    Popreal_Dec2 = ND_PopDec(index(best(1)),:);
end
