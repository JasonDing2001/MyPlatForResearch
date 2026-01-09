classdef PEA3< ALGORITHM
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
            stage = 1;
            %% Optimization
            while Algorithm.NotTerminated(A1)
                           
                % Refresh surrogate models
                if sample_success
                    [Model_obj,Model_con,THETA_obj,THETA_con] = model_train(A2,THETA_obj,THETA_con,stage);
                end


                %Two-stage Optimization
                if stage == 0
                    [PopDec,PopObj,ObjMSE] = Uoptimization(A1,wmax,Model_obj,Problem);
                    PopNew = NewSelect(PopDec,PopObj,[],ObjMSE,[],A2,mu,Problem,stage);
                else
                    [PopDec,PopObj,PopCon,ObjMSE,ConMSE] = Coptimization(A1,wmax,Model_obj,Model_con,Problem);
                    % Select mu solutions for re-evaluation
                    PopNew = NewSelect(PopDec,PopObj,PopCon,ObjMSE,ConMSE,A2,mu,Problem,stage);
                end

                sample_success = 0;
                if isempty(PopNew) == 0
                    PopNew         = Problem.Evaluation(PopNew);
                    Q1 = Cal_Q(PopNew.objs,A2.objs);
                    Q2 = Cal_Q(A2.objs);
                    A2             = [A2,PopNew];
                    index          = EnvironmentalSelection(A2.objs,A2.cons,NI);
                    A1             = A2(index);
                    sample_success = 1;
                end
                
                minQ1 = min(Q1(:));
                minQ2 = min(Q2(:));
                if minQ1 >= minQ2 || Problem.FE<=0.5*Problem.maxFE
                    stage = 1;

                else 
                    stage = 0;
                end
                disp('Minimum value of Q1:');
                disp(minQ1);
                disp('Minimum value of Q2:');
                disp(minQ2);
                disp('Stage:');
                disp(stage);
                disp('FE:');
                disp(Problem.FE);
            end
        end
    end
end

function [Model_obj,Model_con,THETA_obj,THETA_con] = model_train(A2,THETA_obj,THETA_con,flag)
    Dec = A2.decs;
    Obj = A2.objs;
    Len_dec = size(Dec,2);
    Len_obj = size(Obj,2);
    

    for i = 1 : Len_obj
        [~,distinct1] = unique(round(Dec*1e100)/1e100,'rows');
        [~,distinct2] = unique(round(Obj(:,i)*1e100)/1e100,'rows');
        distinct = intersect(distinct1,distinct2);
        
        dmodel     = dacefit(Dec(distinct,:),Obj(distinct,i),'regpoly1','corrgauss',THETA_obj(i,:),1e-5.*ones(1,Len_dec),100.*ones(1,Len_dec));
        Model_obj{i}   = dmodel;
        THETA_obj(i,:) = dmodel.theta;
    end

    if flag == 1                
        Con = A2.cons;
        Len_con = size(Con,2); 
        for i = 1 : Len_con
            [~,distinct1] = unique(round(Dec*1e100)/1e100,'rows');
            [~,distinct2] = unique(round(Con(:,i)*1e100)/1e100,'rows');
            distinct = intersect(distinct1,distinct2);
            
            dmodel     = dacefit(Dec(distinct,:),Con(distinct,i),'regpoly1','corrgauss',THETA_con(i,:),1e-5.*ones(1,Len_dec),100.*ones(1,Len_dec));
            Model_con{i}   = dmodel;
            THETA_con(i,:) = dmodel.theta;
        end
    elseif flag == 0
        Model_con = [];
        % THETA_con = THETA_con;
    end
end

