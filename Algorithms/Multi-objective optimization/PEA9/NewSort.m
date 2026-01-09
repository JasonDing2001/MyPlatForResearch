function [FrontNo, MaxFNo] = NewSort(PopObj, PopCon,nSort)
    % 1. 将解分为可行解和不可行解
    CV = sum(max(PopCon, 0), 2);
    feas = CV == 0;   % 可行解：约束违反为0
    infeas = CV > 0;  % 不可行解：约束违反>0

    % 2. 对可行解进行非支配排序
    ObjFeas = PopObj(feas, :); % 可行解的目标值
    [FrontFeas, MaxFeas] = NDSort(ObjFeas, 20); % 非支配排序（可行解）
    NonDomSolutions = ObjFeas(FrontFeas <= MaxFeas, :);

    % 3. 找出不被非支配可行解支配的不可行解
    ObjInfeas = PopObj(infeas, :); % 不可行解的目标值
    ConInfeas = PopCon(infeas, :); % 不可行解的约束违反
    nonDomInfeas = false(size(ObjInfeas, 1), 1);
    for i = 1:size(ObjInfeas, 1)
        if ~any(FrontFeas < inf & domination(ObjInfeas(i,:), NonDomSolutions))
            nonDomInfeas(i) = true;
        end
    end

    % 4. 对不被支配的不可行解进行带约束的非支配排序
    ObjNonDomInfeas = ObjInfeas(nonDomInfeas, :);
    ConNonDomInfeas = ConInfeas(nonDomInfeas, :);
    [FrontNonDomInfeas, MaxNonDomInfeas] = NDSort(ObjNonDomInfeas, ConNonDomInfeas, inf);

    % 5. 初始化FrontNo
    FrontNo = inf(size(PopObj, 1), 1);

    % 6. 将非支配可行解放在最前面
    FrontNo(feas) = FrontFeas;
    MaxFNo = MaxFeas;

    % 7. 将不被非支配可行解支配的不可行解放在后面
    infeasIdx = find(infeas); % 获取不可行解的索引
    FrontNo(infeasIdx(nonDomInfeas)) = FrontNonDomInfeas + MaxFNo;
    MaxFNo = MaxFNo + MaxNonDomInfeas;

    % 8. 将其他可行解放在后面
    remFeas = find(feas & FrontNo == inf);
    if ~isempty(remFeas)
        [FrontRemFeas, MaxRemFeas] = NDSort(PopObj(remFeas, :), inf);
        FrontNo(remFeas) = FrontRemFeas + MaxFNo;
        MaxFNo = MaxFNo + MaxRemFeas;
    end

    % 9. 将其他不可行解放在最后
    remInfeas = find(infeas); % 获取所有不可行解的索引
    remInfeas = remInfeas(~nonDomInfeas); % 获取其他不可行解的索引
    if ~isempty(remInfeas)
        [~, infeasOrderRem] = sort(CV(remInfeas), 'ascend');
        FrontNo(remInfeas(infeasOrderRem)) = MaxFNo + (1:length(remInfeas))';
        MaxFNo = MaxFNo + length(remInfeas);
    end
    % 只保留前 nSort 个解
    if nSort < length(FrontNo)
        [~, sortedIdx] = sort(FrontNo);
        thresholdFront = FrontNo(sortedIdx(nSort)); % 找到第 nSort 个解所处的支配层级
        FrontNo(FrontNo > thresholdFront) = inf; % 将大于该支配层级的解置为 inf
        MaxFNo = max(FrontNo(FrontNo < inf));
    end
end

function isDom = domination(sol1, solSet)
    % 检查sol1是否被solSet中的任何一个解支配
    isDom = any(all(solSet <= sol1, 2) & any(solSet < sol1, 2));
end