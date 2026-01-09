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
            [TSObjC,PopObjC] = AppendNormalizedCV(TSObj,TSCon,PopObj,PopCon);
            [ND_TSDec,ND_TSObj,ND_PopDec,ND_PopObj] = GetNondominated(TSDec,TSObjC,PopDec,PopObjC);
        else
            [ND_TSDec,ND_TSObj,ND_PopDec,ND_PopObj] = GetNondominated(TSDec,TSObj,PopDec,PopObj);
        end

        if isempty(ND_PopDec)
            selDec = PopDec(randi(size(PopDec,1)),:);
        else
            switch action
                case 1
                    selDec = SelectByAngle(ND_TSObj,ND_PopDec,ND_PopObj);
                case 2
                    selDec = SelectByDeltaPBI(ND_TSObj,ND_TSDec,ND_PopObj,ND_PopDec,W);
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
                selDec = PopDec(randi(size(PopDec,1)),:);
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
