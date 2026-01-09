function PopNew = CNewSelect(PopDec,PopObj,PopCon,ObjMSE,ConMSE,A1,A2,mu,Problem,stage)
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
    A2Obj = A2.objs;
    zmin   = min([PopObj],[],1); zmax = max([PopObj],[],1);
    A1Obj  = (A1Obj - zmin)./max(zmax - zmin,10e-10);
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
   

    if  isempty(regDec)
        
        %%%%  有点问题，当没有解的时候，选择在该方向上可能最好的解？
        %%% 用Tchebycheff 标量调整多样性？
        %%%
        
        tch = max(abs(PopObj*W(regId, :)'),[],1);
        bestidx = find(tch == min(tch));
        PopNew = [PopNew; PopDec(bestidx, :)];


    else

        regCon = PopCon(assignIndex == regId,:);
        regConMSE  = ConMSE(assignIndex == regId,:);
        

        if stage == 2
            [FrontNo,MaxFNo] = NDSort([regObj,sum(max(0,regCon),2)],1);
        elseif stage == 3
            [FrontNo,MaxFNo] = NDSort(regObj,regCon,1);
        end
        Next = FrontNo <= MaxFNo;
        PopNew_ = regDec(Next,:);
        newObj = regObj(Next,:);

        
        Distance = pdist2(newObj,PopObj);
        Distance = sort(Distance,2);
        D = Distance(:,floor(sqrt(size(PopObj,1))));
        [~,index] = max(D);
        PopNew_ = PopNew_(index,:);

        dist2 = pdist2(real(PopNew_(1,:)),real(A2.decs));
        if min(dist2) > 1e-5
            PopNew = [PopNew;PopNew_(1,:)];
        end
       
    end

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
