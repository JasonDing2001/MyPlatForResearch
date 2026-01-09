function Pop = EnvironmentalSelectionCC(Pop,N,isConstrained)
% 基于 SPEA2 的环境选择（复用 CCMO 思路）

    % 删除重复决策向量
    [~,idx] = unique(Pop.Dec,'rows');
    Pop.Dec = Pop.Dec(idx,:);
    Pop.Obj = Pop.Obj(idx,:);
    Pop.Con = Pop.Con(idx,:);

    if isConstrained && ~isempty(Pop.Con)
        Fitness = CalFitnessCC(Pop.Obj,Pop.Con);
    else
        Fitness = CalFitnessCC(Pop.Obj);
    end

    Next = Fitness < 1;
    if sum(Next) < N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    elseif sum(Next) > N
        Del  = TruncationCC(Pop.Obj(Next,:),sum(Next)-N);
        Temp = find(Next);
        Next(Temp(Del)) = false;
    end

    Pop.Dec = Pop.Dec(Next,:);
    Pop.Obj = Pop.Obj(Next,:);
    Pop.Con = Pop.Con(Next,:);
end
