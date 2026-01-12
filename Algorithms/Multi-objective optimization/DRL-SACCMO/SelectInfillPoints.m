function NewDec = SelectInfillPoints(action,Pop,ArchiveAll,W,Model,CModel,Problem,mu,isConstrained)
% 根据动作从候选池中选择真实评估补点

    NewDec = [];
    PopDec = Pop.Dec;
    PopObj = Pop.Obj;
    PopCon = Pop.Con;

    TSDec = ArchiveAll.decs;
    TSObj = ArchiveAll.objs;
    TSCon = ArchiveAll.cons;

    Popreal = [];
    for k = 1 : mu
        if isempty(PopDec)
            break;
        end

        if isConstrained
            [ND_TSDec,ND_TSObj,ND_TSCon,ND_PopDec,ND_PopObj,ND_PopCon] = GetNondominatedFeasible(TSDec,TSObj,TSCon,PopDec,PopObj,PopCon);
        else
            [ND_TSDec,ND_TSObj,ND_PopDec,ND_PopObj] = GetNondominated(TSDec,TSObj,PopDec,PopObj);
        end

        if isempty(ND_PopDec)
            if isConstrained
                selDec = SelectByMinCV(PopDec,PopCon);
            else
                selDec = PopDec(randi(size(PopDec,1)),:);
            end
        else
            switch action
                case 1
                    if isConstrained
                        [F_TSDec,F_TSObj,F_PopDec,F_PopObj] = FilterByFeasibility(ND_TSDec,ND_TSObj,ND_TSCon,ND_PopDec,ND_PopObj,ND_PopCon);
                        selDec = SelectByAngle(F_TSObj,F_PopDec,F_PopObj);
                    else
                        selDec = SelectByAngle(ND_TSObj,ND_PopDec,ND_PopObj);
                    end
                case 2
                    if isConstrained
                        [F_TSDec,F_TSObj,F_PopDec,F_PopObj] = FilterByFeasibility(ND_TSDec,ND_TSObj,ND_TSCon,ND_PopDec,ND_PopObj,ND_PopCon);
                        selDec = SelectByDeltaPBI(F_TSObj,F_TSDec,F_PopObj,F_PopDec,W);
                    else
                        selDec = SelectByDeltaPBI(ND_TSObj,ND_TSDec,ND_PopObj,ND_PopDec,W);
                    end
                case 3
                    if isConstrained
                        selDec = SelectByEPDM_C(TSDec,TSObj,TSCon,PopDec,PopObj,PopCon,Popreal,W,Model,CModel,Problem);
                    else
                        selDec = SelectByEPDM_U(TSDec,TSObj,PopDec,PopObj,Popreal,W,Model,Problem);
                    end
                otherwise
                    selDec = PopDec(randi(size(PopDec,1)),:);
            end
            if isempty(selDec)
                if isConstrained
                    selDec = SelectByMinCV(PopDec,PopCon);
                else
                    selDec = PopDec(randi(size(PopDec,1)),:);
                end
            end
        end

        NewDec = [NewDec;selDec];
        Popreal = [Popreal;selDec];

        idx = find(ismember(PopDec,selDec,'rows'),1);
        if ~isempty(idx)
            PopDec(idx,:) = [];
            PopObj(idx,:) = [];
            PopCon(idx,:) = [];
        end
    end
end

function [ND_TSDec,ND_TSObj,ND_TSCon,ND_PopDec,ND_PopObj,ND_PopCon] = GetNondominatedFeasible(TSDec,TSObj,TSCon,PopDec,PopObj,PopCon)
% 基于可行性优先原则的非支配筛选

    AllDec = [TSDec;PopDec];
    AllObj = [TSObj;PopObj];
    AllCon = [TSCon;PopCon];
    CV     = sum(max(0,AllCon),2);

    N = size(AllObj,1);
    dominated = false(N,1);
    for i = 1:N
        if dominated(i)
            continue;
        end
        for j = 1:N
            if i == j || dominated(j)
                continue;
            end
            if DominatesDeb(AllObj(i,:),CV(i),AllObj(j,:),CV(j))
                dominated(j) = true;
            elseif DominatesDeb(AllObj(j,:),CV(j),AllObj(i,:),CV(i))
                dominated(i) = true;
                break;
            end
        end
    end

    ND = ~dominated;
    nTS = size(TSDec,1);
    idxTS = find(ND(1:nTS));
    idxPop = find(ND(nTS+1:end)) + nTS;

    ND_TSDec = AllDec(idxTS,:);
    ND_TSObj = AllObj(idxTS,:);
    ND_TSCon = AllCon(idxTS,:);
    ND_PopDec = AllDec(idxPop,:);
    ND_PopObj = AllObj(idxPop,:);
    ND_PopCon = AllCon(idxPop,:);
end

function flag = DominatesDeb(objA,cvA,objB,cvB)
% Deb 约束可行性优先支配关系

    if cvA == 0 && cvB > 0
        flag = true;
    elseif cvA > 0 && cvB == 0
        flag = false;
    elseif cvA == 0 && cvB == 0
        flag = all(objA <= objB) && any(objA < objB);
    else
        flag = cvA < cvB;
    end
end

function selDec = SelectByMinCV(PopDec,PopCon)
% 不可行时优先选择 CV 最小的候选

    if isempty(PopDec)
        selDec = [];
        return;
    end
    CV = sum(max(0,PopCon),2);
    [~,idx] = min(CV);
    selDec = PopDec(idx,:);
end

function [TSDec,TSObj,PopDec,PopObj] = FilterByFeasibility(TSDec,TSObj,TSCon,PopDec,PopObj,PopCon)
% 先可行，若全不可行则保留最小 CV

    if ~isempty(TSCon)
        CV_TS = sum(max(0,TSCon),2);
        if any(CV_TS == 0)
            maskTS = CV_TS == 0;
        else
            minCV = min(CV_TS);
            maskTS = CV_TS == minCV;
        end
        TSDec = TSDec(maskTS,:);
        TSObj = TSObj(maskTS,:);
    end

    if ~isempty(PopCon)
        CV_Pop = sum(max(0,PopCon),2);
        if any(CV_Pop == 0)
            maskPop = CV_Pop == 0;
        else
            minCV = min(CV_Pop);
            maskPop = CV_Pop == minCV;
        end
        PopDec = PopDec(maskPop,:);
        PopObj = PopObj(maskPop,:);
    end
end
