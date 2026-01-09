function [ND_TSDec,ND_TSObj,ND_PopDec,ND_PopObj] = GetNondominated(TSDec,TSObj,PopDec,PopObj)
% 提取非支配解集合

    if isempty(TSObj)
        ND_TSDec = TSDec;
        ND_TSObj = TSObj;
    else
        [FrontNo,~] = NDSort(TSObj,inf);
        ND_TSDec = TSDec(FrontNo==1,:);
        ND_TSObj = TSObj(FrontNo==1,:);
    end

    if isempty(PopObj)
        ND_PopDec = [];
        ND_PopObj = [];
    else
        [FrontNo,~] = NDSort(PopObj,inf);
        ND_PopDec = PopDec(FrontNo==1,:);
        ND_PopObj = PopObj(FrontNo==1,:);
    end
end
