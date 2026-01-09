classdef PEA2< ALGORITHM
% <multi/many> <real> <constrained/none>
% Surrogate-assisted RVEA
% wmax  --- 20 --- Number of generations before updating Kriging models
% mu    ---  5 --- Number of re-evaluated solutions at each generation

%------------------------------- Reference --------------------------------
% T. Chugh, Y. Jin, K. Miettinen, J. Hakanen, and K. Sindhya, A surrogate-
% assisted reference vector guided evolutionary algorithm for
% computationally expensive many-objective optimization, IEEE Transactions
% on Evolutionary Computation, 2018, 22(1): 129-142.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [wmax,mu] = Algorithm.ParameterSet(20,5);
            %% Generate the reference points and population
            NI             = Problem.N;
            P              = UniformPoint(NI,Problem.D,'Latin');
            A2             = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));
            A1             = A2; 
            THETA_obj      = 5.*ones(Problem.M,Problem.D);
            THETA_con      = 5.*ones(size(A2.cons,2),Problem.D);
            Model_obj      = cell(1,Problem.M);
            Model_con      = cell(1,size(A2.cons,2));
            sample_success = 1;
            %% Optimization
            while Algorithm.NotTerminated(A2.best)
                % Refresh surrogate models
                if sample_success
                    [Model_obj,Model_con,THETA_obj,THETA_con] = model_train(A2,THETA_obj,THETA_con);
                end
                
                %Optimization
                [PopDec,PopObj,PopCon,ObjMSE,ConMSE] = optimizaiton(A1,wmax,Model_obj,Model_con,Problem);
                 
                % Select mu solutions for re-evaluation
                PopNew = NewSelect(PopDec,PopObj,PopCon,ObjMSE,ConMSE,A2,mu,Problem);
                
                sample_success = 0;
                if isempty(PopNew) == 0
                    PopNew         = Problem.Evaluation(PopNew);
                    A2             = [A2,PopNew];
                    index          = EnvironmentalSelection(A2.objs,A2.cons,NI);
                    A1             = A2(index);
                    sample_success = 1;
                end
            end
        end
    end
end

function [Model_obj,Model_con,THETA_obj,THETA_con] = model_train(A2,THETA_obj,THETA_con)
    Dec = A2.decs;
    Obj = A2.objs;
    Con = A2.cons;
    Len_dec = size(Dec,2);
    Len_obj = size(Obj,2);
    Len_con = size(Con,2);
    for i = 1 : Len_obj
        [~,distinct1] = unique(round(Dec*1e100)/1e100,'rows');
        [~,distinct2] = unique(round(Obj(:,i)*1e100)/1e100,'rows');
        distinct = intersect(distinct1,distinct2);
        
        dmodel     = dacefit(Dec(distinct,:),Obj(distinct,i),'regpoly1','corrgauss',THETA_obj(i,:),1e-5.*ones(1,Len_dec),100.*ones(1,Len_dec));
        Model_obj{i}   = dmodel;
        THETA_obj(i,:) = dmodel.theta;
    end
    for i = 1 : Len_con
        [~,distinct1] = unique(round(Dec*1e100)/1e100,'rows');
        [~,distinct2] = unique(round(Con(:,i)*1e100)/1e100,'rows');
        distinct = intersect(distinct1,distinct2);
        
        dmodel     = dacefit(Dec(distinct,:),Con(distinct,i),'regpoly1','corrgauss',THETA_con(i,:),1e-5.*ones(1,Len_dec),100.*ones(1,Len_dec));
        Model_con{i}   = dmodel;
        THETA_con(i,:) = dmodel.theta;
    end
end

function [OffObj,Off_ObjMSE,OffCon,Off_ConMSE] = model_predict(Model_obj,Model_con,OffDec)
    N          = size(OffDec,1);
    Len_obj    = length(Model_obj);
    Len_con    = length(Model_con);
    OffObj     = zeros(N,Len_obj);
    OffCon     = zeros(N,Len_con);
    Off_ObjMSE = zeros(N,Len_obj);
    Off_ConMSE = zeros(N,Len_con);
      
    for i = 1 : N
        for j = 1 : Len_obj
            [OffObj(i,j),~,Off_ObjMSE(i,j)] = predictor(OffDec(i,:),Model_obj{j});
        end
        for j = 1 : Len_con
            [OffCon(i,j),~,Off_ConMSE(i,j)] = predictor(OffDec(i,:),Model_con{j});
        end
    end
end

function [PopDec,PopObj,PopCon,ObjMSE,ConMSE] = optimizaiton(A1,wmax,Model_obj,Model_con,Problem)

    PopDec = A1.decs;
    PopObj = A1.objs;
    PopCon = A1.cons;
    ObjMSE = zeros(Problem.N,Problem.M);
    ConMSE = zeros(Problem.N,size(PopCon,2));
    w      = 1;
    while w <= wmax
        drawnow();
        OffDec = OperatorGA(Problem,PopDec);
        [OffObj,Off_ObjMSE,OffCon,Off_ConMSE] = model_predict(Model_obj,Model_con,OffDec);

        PopDec = [PopDec;OffDec];
        PopObj = [PopObj;OffObj];
        PopCon = [PopCon;OffCon];
        ObjMSE = [ObjMSE;Off_ObjMSE];
        ConMSE = [ConMSE;Off_ConMSE];
        
        index  = SEnvironmentalSelection(PopObj,ObjMSE,PopCon,ConMSE,length(A1));
      
        PopDec = PopDec(index,:);
        PopObj = PopObj(index,:);
        PopCon = PopCon(index,:);
        ObjMSE = ObjMSE(index,:);
        ConMSE = ConMSE(index,:);
        w = w + 1;
    end
end

