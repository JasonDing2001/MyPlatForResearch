function [Model,CModel,THETA_OBJ,THETA_CON] = TrainKrigingModels(Population,Problem,THETA_OBJ,THETA_CON)
% 训练目标与约束的 Kriging 模型（dacefit）

    PopDec = Population.decs;
    PopObj = Population.objs;
    PopCon = Population.cons;
    numC   = size(PopCon,2);

    % TODO: 可根据样本规模限制训练集大小，提高效率
    TrainDec = PopDec;
    TrainObj = PopObj;
    TrainCon = PopCon;

    Model = cell(1,Problem.M);
    for i = 1 : Problem.M
        dmodel         = dacefit(TrainDec,TrainObj(:,i),'regpoly0','corrgauss',THETA_OBJ(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
        Model{i}       = dmodel;
        THETA_OBJ(i,:) = dmodel.theta;
    end

    if numC > 0
        if isempty(THETA_CON) || size(THETA_CON,1) ~= numC
            THETA_CON = 5.*ones(numC,Problem.D);
        end
        CModel = cell(1,numC);
        for i = 1 : numC
            cdmodel         = dacefit(TrainDec,TrainCon(:,i),'regpoly0','corrgauss',THETA_CON(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
            CModel{i}       = cdmodel;
            THETA_CON(i,:)  = cdmodel.theta;
        end
    else
        CModel    = {};
        THETA_CON = [];
    end
end
