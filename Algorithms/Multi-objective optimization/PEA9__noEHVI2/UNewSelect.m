function PopNew = NewSelect1(PopDec,PopObj,ObjMSE,A1,A2,mu,Problem)
    %% Preparing Data
    index = ismember(PopDec, A1.decs, 'rows');
    PopNew = [];
    if sum(index) == size(PopDec, 1)
        return;
    elseif sum(~index) <= mu
        PopNew_ = PopDec(~index, :);
        for i = 1:size(PopNew_, 1)
            dist2 = pdist2(real(PopNew_(i, :)), real(A1.decs));
            if min(dist2) > 1e-5
                PopNew = [PopNew; PopNew_(i, :)];
            end
        end
        return;
    end

    % Normalize the objectives
    PopDec = PopDec(~index, :);
    PopObj = PopObj(~index, :); ObjMSE = ObjMSE(~index, :);
    A1Obj = A1.objs;
    A2Obj = A2.objs;
    zmin = min([PopObj], [], 1); zmax = max([PopObj], [], 1);
    A1Obj = (A1Obj - zmin) ./ max(zmax - zmin, 10e-10);
    PopObj = (PopObj - zmin) ./ max(zmax - zmin, 10e-10);

    [FrontNo, MaxFNo] = NDSort(PopObj, mu);
    Next = FrontNo < MaxFNo;
    % Calculate crowding distance

    CrowdDis = CrowdingDistance(PopObj, FrontNo);
    % Select the solutions in the last front

    % Select the solutions in the last front
    Last = find(FrontNo == MaxFNo);
    [~, Rank] = sort(CrowdDis(Last), 'descend');
    Next(Last(Rank(1:mu - sum(Next)))) = true;

    % Combine the selected solutions
    PopNew = PopDec(Next, :);
end

