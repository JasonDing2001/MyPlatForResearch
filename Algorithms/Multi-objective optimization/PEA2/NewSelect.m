function PopNew = NewSelect(PopDec,PopObj,PopCon,ObjMSE,ConMSE,A1,mu,Problem)
    %% Preparing Data
    index   = ismember(PopDec,A1.decs,'rows');
    if sum(index) == size(PopDec,1)
        PopNew  = [];
        return;
    elseif sum(~index) <= mu
        PopNew_ = PopDec(~index,:);
        PopNew  = [];
        for i = 1:size(PopNew_,1)
            dist2 = pdist2(real(PopNew_(i,:)),real(A1.decs));
            if min(dist2) > 1e-5
                PopNew = [PopNew;PopNew_(i,:)];
            end
        end
        return;
    end
    
    PopDec = PopDec(~index,:);
    PopObj = PopObj(~index,:); ObjMSE  = ObjMSE(~index,:);
    PopCon = PopCon(~index,:); ConMSE  = ConMSE(~index,:);

    A2Obj  = A1.objs;
    A2Con  = A1.cons;
    zmin   = min([A2Obj;PopObj],[],1); zmax = max([A2Obj;PopObj],[],1);
    A2Obj  = (A2Obj - zmin )./max(zmax - zmin,10e-10);
    PopObj = (PopObj - zmin)./max(zmax - zmin,10e-10);
    ObjMSE = ObjMSE./(max(zmax - zmin,10e-10).^2);
    
    %% Reference Set 
 
    num = length(find(all(A2Con<=0,2)));
    [FrontNo,~] = NDSort(A2Obj,A2Con,inf);
    if num >= Problem.N
        A2Obj = A2Obj(FrontNo==1,:);
    else
        i = 1;
        Next = FrontNo == i;
        while sum(Next) < Problem.N
            Next(FrontNo == i) = true;
            i = i + 1;
        end
        A2Obj = A2Obj(Next,:);
    end

    %% Select mu points
    [FrontNo,MaxFNo] = NDSort_CPPD(PopObj,ObjMSE,PopCon,ConMSE,mu);
    Next = FrontNo < MaxFNo;
    Last = find(FrontNo == MaxFNo);

    if length(Last) == mu - sum(Next)
        Next(Last) = true;
    elseif length(Last) > mu - sum(Next)
        A2Obj = [A2Obj;PopObj(Next,:)];
        for i = 1 : mu - sum(Next)
            index = EucDistance_Select(PopObj(Last,:),A2Obj);
            Next(Last(index)) = true;
            A2Obj = [A2Obj;PopObj(Last(index),:)];
            Last(index)=[];
        end
    end
    
    %% Output
    PopNew_ = PopDec(Next,:);
    PopNew  = [];
    for i = 1:size(PopNew_,1)
        dist2 = pdist2(real(PopNew_(i,:)),real(A1.decs));
        if min(dist2) > 1e-5
            PopNew = [PopNew;PopNew_(i,:)];
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
