classdef PEA9nonewsort< ALGORITHM
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
            
            %% A2 for Archive, A1 for Population.
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

                Q          = Cal_Q(A2.objs); % calculate the convergence contribution of all evaluated solutions
                [Q, index] = sort(Q,'descend');
                CV         = sum(max(A2.cons,0),2); % calculate the CV of all evaluated solutions
                CV         = CV(index);
                coef = corrcoef(Q(end-mu:end), CV(end-mu:end)); % calculate the correlation coefficient
                r_coef = coef(2);
                phi = -0.2;
                tau = 0.6;


                % adaptive switching
                if (r_coef < phi || Problem.FE >= 0.7 * Problem.maxFE) && (Problem.FE - NI) > Problem.maxFE * 0.1
                    stage = 3;   %  Constrained surrogate-assisted evolutionary search
                else
                    if r_coef <  tau && Problem.FE < 0.7 * Problem.maxFE && (Problem.FE - NI) > Problem.maxFE * 0.1
                        stage = 2; 
                    else
                        stage = 1;
                    end
                end
                %stage = 1;

                % Refresh surrogate models
                if sample_success
                    [Model_obj,Model_con,THETA_obj,THETA_con] = model_train(A2,THETA_obj,THETA_con,stage);
                end
                
                % Three-stage Optimization
                if stage == 1 
                    [PopDec,PopObj,ObjMSE] = Uoptimization(A1,wmax,Model_obj,Problem);
                    PopNew = UNewSelect(PopDec,PopObj,ObjMSE,A1,A2,mu,Problem);
                elseif stage == 2
                    [PopDec,PopObj,PopCon,ObjMSE,ConMSE] = Coptimization(A1,wmax,Model_obj,Model_con,Problem,stage);
                    % PopNew = CNewSelect(PopDec,PopObj,PopCon,ObjMSE,ConMSE,A1,A2,mu,Problem,stage);
                    PopNew = CNewSelect(PopDec,PopObj,PopCon,ObjMSE,ConMSE,A1,A2,mu,Problem,stage);
                elseif stage == 3
                    [PopDec,PopObj,PopCon,ObjMSE,ConMSE] = Coptimization(A1,wmax,Model_obj,Model_con,Problem,stage);
                    % PopNew = CNewSelect(PopDec,PopObj,PopCon,ObjMSE,ConMSE,A1,A2,mu,Problem,stage);
                    PopNew = CNewSelect1(PopDec,PopObj,PopCon,ObjMSE,ConMSE,A1,A2,mu,Problem,ObjMSE);
                end
                % Select mu solutions for re-evaluation
                %PopNew = NewSelect(PopDec,PopObj,PopCon,ObjMSE,ConMSE,A2,mu,Problem);
                
                sample_success = 0;
                if isempty(PopNew) == 0
                    PopNew         = Problem.Evaluation(PopNew);
                    A2             = [A2,PopNew];
                    % A1             = [A1,PopNew];
                    index          = EnvironmentalSelection(A2.objs,A2.cons,NI,stage);
                    A1             = A2(index);
                    sample_success = 1;
                end

                disp(stage);
                disp(Problem.FE);
            end
        end
    end
end

function value = Cal_Q(Obj)
    N = size(Obj,1);
    Obj = (Obj-repmat(min(Obj),N,1))./(repmat(max(Obj)-min(Obj),N,1));
    I = zeros(N);
    for i = 1 : N
        for j = 1 : N
            I(i,j) = max(Obj(i,:)-Obj(j,:));
        end
    end
    C = max(abs(I));
    F = sum(-exp(-I./repmat(C,N,1)/0.05)) + 1;
    value = 1./F;
end