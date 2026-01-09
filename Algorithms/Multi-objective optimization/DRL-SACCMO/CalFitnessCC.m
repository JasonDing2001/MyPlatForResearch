function Fitness = CalFitnessCC(PopObj,PopCon)
% 计算 SPEA2 适应度（可选约束）

    N = size(PopObj,1);
    if nargin == 1
        CV = zeros(N,1);
    else
        CV = sum(max(0,PopCon),2);
    end

    % 支配关系
    Dominate = false(N);
    for i = 1 : N-1
        for j = i+1 : N
            if CV(i) < CV(j)
                Dominate(i,j) = true;
            elseif CV(i) > CV(j)
                Dominate(j,i) = true;
            else
                k = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
                if k == 1
                    Dominate(i,j) = true;
                elseif k == -1
                    Dominate(j,i) = true;
                end
            end
        end
    end

    % 强度
    S = sum(Dominate,2);

    % 原始适应度
    R = zeros(1,N);
    for i = 1 : N
        R(i) = sum(S(Dominate(:,i)));
    end

    % 密度
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Distance = sort(Distance,2);
    D = 1./(Distance(:,floor(sqrt(N)))+2);

    Fitness = R + D';
end
