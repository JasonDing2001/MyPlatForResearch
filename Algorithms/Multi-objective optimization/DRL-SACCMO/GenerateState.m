function [state,action,nextstate,ter,reward] = GenerateState(Problem,action,LastPop,Pop)
% DRL-SAEA 同款状态与奖励

    ter = 0;
    if isempty(LastPop) || isempty(Pop)
        state = zeros(1,4);
        nextstate = state;
        reward = 0;
        if action == -1
            action = 3;
        end
        return;
    end

    Ref = max([LastPop.objs;Pop.objs],[],1);
    LastHV = HV(LastPop,Ref);
    CurrHV = HV(Pop,Ref);
    reward1 = (CurrHV-LastHV)/max(LastHV,eps);
    if isnan(reward1)
        reward1 = 0;
    end

    LastObjsVar = sum(var(LastPop.objs));
    LastObjsCon = sum(LastPop.objs,'all');
    LastCV      = sum(max(0,LastPop.cons),'all');
    LastRatio   = Problem.FE/Problem.maxFE;

    CurrObjsVar = sum(var(Pop.objs));
    CurrObjsCon = sum(Pop.objs,'all');
    CurrCV      = sum(max(0,Pop.cons),'all');
    CurrRatio   = Problem.FE/Problem.maxFE;

    state     = [LastObjsVar,LastObjsCon,LastCV,LastRatio];
    nextstate = [CurrObjsVar,CurrObjsCon,CurrCV,CurrRatio];

    if LastCV == 0
        reward2 = 0;
    elseif CurrCV < LastCV
        reward2 = abs((CurrCV-LastCV)/LastCV);
    else
        reward2 = -abs((CurrCV-LastCV)/LastCV);
    end
    reward = reward1 + reward2;

    if action == -1
        action = 3;
    end
end
