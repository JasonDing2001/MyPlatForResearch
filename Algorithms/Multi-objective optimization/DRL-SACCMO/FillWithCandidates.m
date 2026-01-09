function NewDec = FillWithCandidates(NewDec,PopUDec,PopCDec,mu,Problem)
% 兜底补点，优先使用候选池中的剩余个体

    need = mu - size(NewDec,1);
    if need <= 0
        return;
    end

    Candidate = unique([PopUDec;PopCDec],'rows');
    Candidate = Candidate(~ismember(Candidate,NewDec,'rows'),:);

    if ~isempty(Candidate)
        if size(Candidate,1) >= need
            idx = randperm(size(Candidate,1),need);
            NewDec = [NewDec;Candidate(idx,:)];
            return;
        else
            NewDec = [NewDec;Candidate];
            need = mu - size(NewDec,1);
        end
    end

    if need > 0
        P = lhsamp(need,Problem.D);
        RandDec = repmat(Problem.upper-Problem.lower,need,1).*P + repmat(Problem.lower,need,1);
        NewDec = [NewDec;RandDec];
        % TODO: 可加入修复/去重步骤
    end
end
