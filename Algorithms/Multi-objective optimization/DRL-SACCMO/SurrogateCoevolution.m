function [PopU,PopC] = SurrogateCoevolution(Archive,Problem,Model,CModel,wmax)
% 内层代理演化：U 子任务忽略约束，C 子任务强调可行性

    PopU.Dec = Archive.decs;
    PopU.Obj = Archive.objs;
    PopU.Con = Archive.cons;

    PopC = PopU;

    for gen = 1 : wmax
        % U 子任务选择
        FitU = CalFitnessCC(PopU.Obj);
        MatingPoolU = TournamentSelection(2,Problem.N,FitU);
        OffDecU     = OperatorGA(Problem,PopU.Dec(MatingPoolU,:));

        % C 子任务选择
        FitC = CalFitnessCC(PopC.Obj,PopC.Con);
        MatingPoolC = TournamentSelection(2,Problem.N,FitC);
        OffDecC     = OperatorGA(Problem,PopC.Dec(MatingPoolC,:));

        % 代理预测（目标与约束）
        [OffObjU,~,OffConU,~] = PredictSurrogate(OffDecU,Model,CModel);
        [OffObjC,~,OffConC,~] = PredictSurrogate(OffDecC,Model,CModel);

        % CCMO 协同更新：两个子任务共享后代
        CandU.Dec = [PopU.Dec;OffDecU;OffDecC];
        CandU.Obj = [PopU.Obj;OffObjU;OffObjC];
        CandU.Con = [PopU.Con;OffConU;OffConC];

        CandC.Dec = [PopC.Dec;OffDecU;OffDecC];
        CandC.Obj = [PopC.Obj;OffObjU;OffObjC];
        CandC.Con = [PopC.Con;OffConU;OffConC];

        PopU = EnvironmentalSelectionCC(CandU,Problem.N,false);
        PopC = EnvironmentalSelectionCC(CandC,Problem.N,true);
    end

    % 输出时统一为代理预测值
    [PopU.Obj,~,PopU.Con,~] = PredictSurrogate(PopU.Dec,Model,CModel);
    [PopC.Obj,~,PopC.Con,~] = PredictSurrogate(PopC.Dec,Model,CModel);
end
