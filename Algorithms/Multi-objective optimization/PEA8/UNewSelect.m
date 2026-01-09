function PopNew = NewSelect1(PopDec,PopObj,ObjMSE,A1,A2,mu,Problem)
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
    
    
    % Generate the reference vectors
    [W, numVectors] = UniformPoint(mu,Problem.M);
    
    
    % Normalize the objectives
    PopDec = PopDec(~index,:);
    PopObj = PopObj(~index,:); ObjMSE  = ObjMSE(~index,:);
    A1Obj  = A1.objs;
    A1Con  = A1.cons;
    zmin   = min([A1Obj;PopObj],[],1); zmax = max([A1Obj;PopObj],[],1);
    A1Obj  = (A1Obj - zmin )./max(zmax - zmin,10e-10);
    PopObj = (PopObj - zmin)./max(zmax - zmin,10e-10);

    % Assign each point to its corresponding reference vector
    numPoints = size(PopObj,1);
    assignIndex = zeros(numPoints,1);
    
    for i = 1:numPoints
        angles = zeros(1, numVectors);
        for j = 1:numVectors
            cos_theta = dot(PopObj(i, :), W(j, :)) / (norm(PopObj(i, :)) * norm(W(j, :)));
            angles(j) = acos(cos_theta); % 计算夹角
        end
        [~, minIdx] = min(angles);
        assignIndex(i) = minIdx;
    end
   
    for regId = 1:numVectors
        regDec = PopDec(assignIndex == regId, :);
        regObj = PopObj(assignIndex == regId, :);
        regObjMSE = ObjMSE(assignIndex == regId, :);

        if  isempty(regDec)
            tch = max(abs(PopObj*W(regId, :)'),[],1);
            bestidx = find(tch == min(tch));
            PopNew = [PopNew; PopDec(bestidx, :)];
        else
    
            [FrontNo,MaxFNo] = NDSort(regObj,1);
            Next = FrontNo <= MaxFNo;
            PopNew_ = regDec(Next,:);
            newObj = regObj(Next,:);
            
            % 计算 PopNew_ 中每个解到 A1.decs 中所有解的距离
            dist2_all = pdist2(newObj, real(A2.objs));  
            min_distances = min(dist2_all, [], 2); 
            % 选择距离最远的候选解
            [~, farthest_idx] = max(min_distances);
            Popnew_ = PopNew_(farthest_idx, :);  % 最远解对应的决策变量

            dist2 = pdist2(real(PopNew_(1,:)),real(A2.decs));
            if min(dist2) > 1e-5
                PopNew = [PopNew;PopNew_(1,:)];
            end
           
        end
    end
    % disp(PopNew);
end

