% function PopNew = CNewSelect2(PopDec, PopObj, PopCon, ObjMSE, ConMSE, A1, A2, mu, Problem, PopMSE)
%     %% Preparing Data
%     index = ismember(PopDec, A1.decs, 'rows');
%     PopNew = [];
%     if sum(index) == size(PopDec, 1)
%         return;
%     elseif sum(~index) <= mu
%         PopNew_ = PopDec(~index, :);
%         for i = 1:size(PopNew_, 1)
%             dist2 = pdist2(real(PopNew_(i, :)), real(A1.decs));
%             if min(dist2) > 1e-5
%                 PopNew = [PopNew; PopNew_(i, :)];
%             end
%         end
%         return;
%     end

%     % 排除A1中的个体
%     PopDec = PopDec(~index, :);
%     PopObj = PopObj(~index, :);
%     PopCon = PopCon(~index, :);

%     % 归一化目标值
%     zmin = min(PopObj, [], 1);
%     zmax = max(PopObj, [], 1);
%     PopObj = (PopObj - zmin) ./ max(zmax - zmin, 1e-10);

%     % 计算约束违反总和
%     CV = sum(max(0, PopCon), 2);
%     feasible = CV == 0;

%     % 分离可行和不可行解
%     PopDecFeasible = PopDec(feasible, :);
%     PopObjFeasible = PopObj(feasible, :);
%     PopDecInfeasible = PopDec(~feasible, :);
%     CVInfeasible = CV(~feasible);

%     % 处理可行解的非支配排序和拥挤度选择
%     SelectedFeasible = [];
%     if ~isempty(PopDecFeasible)
%         [FrontNo, MaxFNo] = NDSort(PopObjFeasible, inf);
%         currentCount = 0;
%         Selected = [];
%         for i = 1:MaxFNo
%             currentFront = find(FrontNo == i);
%             numCurrent = length(currentFront);
%             try
%                 % if currentCount + numCurrent <= mu
%                 %     Selected = [Selected; currentFront];
%                 %     currentCount = currentCount + numCurrent;
%                 % else
%                 %     CrowdDis = CrowdingDistance(PopObjFeasible(currentFront, :));
%                 %     [~, idx] = sort(CrowdDis, 'descend');
%                 %     needed = mu - currentCount;
                   
%                 %    Selected = [Selected; currentFront(idx(1:needed))];
    
%                 %     break;
%                 % end
%                 if currentCount + numCurrent <= mu
%                     Selected = [Selected; currentFront];
%                     currentCount = currentCount + numCurrent;
%                 else
%                     CrowdDis = CrowdingDistance(PopObjFeasible(currentFront, :));
%                     [~, idx] = sort(CrowdDis, 'descend');
%                     needed = mu - currentCount;
%                     selectedFromCurrent = currentFront(idx(1:needed));
%                     selectedFromCurrent = selectedFromCurrent(:); % 确保列向量
%                     Selected = [Selected; selectedFromCurrent];
%                     break;
%             catch
%                 keyboard;
%             end
%         end
%         SelectedFeasible = PopDecFeasible(Selected, :);
%     end

%     % 处理不可行解的选择
%     numSelected = size(SelectedFeasible, 1);
%     if numSelected < mu
%         numNeeded = mu - numSelected;
%         [~, idx] = sort(CVInfeasible, 'ascend');
%         if numNeeded > length(idx)
%             numNeeded = length(idx);
%         end
%         SelectedInfeasible = PopDecInfeasible(idx(1:numNeeded), :);
%         PopNew = [SelectedFeasible; SelectedInfeasible];
%     else
%         PopNew = SelectedFeasible;
%     end

%     % 检查是否与A1中的解重复
%     tempPopNew = [];
%     for i = 1:size(PopNew, 1)
%         dist2 = pdist2(real(PopNew(i, :)), real(A1.decs));
%         if min(dist2) > 1e-5
%             tempPopNew = [tempPopNew; PopNew(i, :)];
%         end
%     end
%     PopNew = tempPopNew;

%     % 确保返回不超过mu个解
%     if size(PopNew, 1) > mu
%         PopNew = PopNew(1:mu, :);
%     end
% end

% function CrowdDis = CrowdingDistance(PopObj)
%     [N, M] = size(PopObj);
%     if N == 0
%         CrowdDis = [];
%         return;
%     end
%     CrowdDis = zeros(N, 1);
%     for m = 1:M
%         [sorted, rank] = sort(PopObj(:, m));
%         CrowdDis(rank(1)) = Inf;
%         CrowdDis(rank(end)) = Inf;
%         denominator = sorted(end) - sorted(1);
%         if denominator == 0
%             continue;
%         end
%         for i = 2:(N-1)
%             CrowdDis(rank(i)) = CrowdDis(rank(i)) + (sorted(i+1) - sorted(i-1)) / denominator;
%         end
%     end
% end


function PopNew = CNewSelect2(PopDec, PopObj, PopCon, ObjMSE, ConMSE, A1, A2, mu, Problem, PopMSE)
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
    PopCon = PopCon(~index,:); ConMSE = ConMSE(~index,:);
    A1Obj = A1.objs;
    A1Con = A1.cons;
    A2Obj = A2.objs;
    zmin = min([PopObj], [], 1); zmax = max([PopObj], [], 1);
    A1Obj = (A1Obj - zmin) ./ max(zmax - zmin, 10e-10);
    PopObj = (PopObj - zmin) ./ max(zmax - zmin, 10e-10);

    [FrontNo, MaxFNo] = NDSort(PopObj, PopCon, mu);
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
