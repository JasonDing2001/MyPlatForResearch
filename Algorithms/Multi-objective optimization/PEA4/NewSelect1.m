function PopNew = NewSelect1(PopDec,PopObj,PopCon,ObjMSE,ConMSE,A1,mu,Problem,stage)
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
            [~, nearest_idx] = min(pdist2(PopObj, W(regId, :)));
            PopNew = [PopNew; PopDec(nearest_idx, :)];
        else
            if stage == 1
                regCon = PopCon(assignIndex == regId,:);
                regConMSE  = ConMSE(assignIndex == regId,:);
                [FrontNo,MaxFNo] = NDSort_CPPD(regObj,regObjMSE,regCon,regConMSE,1);
            elseif stage == 0 
                [FrontNo,MaxFNo] = UNDSort(regObj,regObjMSE,1);
            end
            Next = FrontNo <= MaxFNo;
            PopNew_ = regDec(Next,:);
            newObj = regObj(Next,:);
            
            % 计算 PopNew_ 中每个解到 A1.decs 中所有解的距离
            dist2_all = pdist2(newObj, real(A1.objs));  
            min_distances = min(dist2_all, [], 2); 
            % 选择距离最远的候选解
            [~, farthest_idx] = max(min_distances);
            Popnew_ = PopNew_(farthest_idx, :);  % 最远解对应的决策变量

            dist2 = pdist2(real(PopNew_(1,:)),real(A1.decs));
            if min(dist2) > 1e-5
                PopNew = [PopNew;PopNew_(1,:)];
            end
           
        end
    end
    disp(PopNew);
end


function index = EucDistance_Select(PopObj,ALL_Obj)
    N1 = size(PopObj,1);
    N2 = size(ALL_Obj,1);
    Distance = zeros(N1,N2);
    %% Calculate the shifted distance between each two solutions
    for i = 1 : N1
        for j = 1 : N2
            Distance(i,j) = norm(PopObj(i,:)-ALL_Obj(j,:),2);    
        end
    end
    Distance = sort(Distance,2);
    [~,index] = max(Distance(:,1));
end
