classdef EICMSSAEA_Tr3 < ALGORITHM
    % <multi/many> <real/integer/label/binary/permutation> <constrained><expensive>
    % Surrogate-assisted RVEA
    % alpha ---  2 --- The parameter controlling the rate of change of penalty
    % wmax  --- 20 --- Number of generations before updating Kriging models
    % mu    ---  5 --- Number of re-evaluated INDIVIDUALs at each generation
    % This function is written by Cheng He
    methods
        function main(Algorithm, Problem)
            %% Generate the reference points and population
            [alpha,~,mu] = Algorithm.ParameterSet(5,20,5);
            [V0,Problem.N] = UniformPoint(Problem.N,Problem.M);
            [V01,Problem.N] = UniformPoint(Problem.N,Problem.M+1);
            NI    = 11*Problem.D-1;
            P     = lhsamp(NI,Problem.D);
            % P  = lhsdesign(NI, Problem.D,'criterion','maximin','iterations',1000) ;
            % Initial expensive evaluations
            A2    = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));
            A1    = A2;
            numC = size(A2.cons,2);
            % GP hyper-parameters for objectives and constraints
            THETA = 5.*ones(Problem.M,Problem.D);
            CTHETA = 5.*ones(numC,Problem.D);
            Model = cell(1,Problem.M);
            CModel = cell(1,numC);
            %[~,FR] = Feasible_rate(A1);
            Arc_CV = [];Arc_meanCV= [];
            tr = 10;iter=1;% LIRCMOP 10 MW 20
            Stageflag=1;t = 1;
            while Algorithm.NotTerminated(A2)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Train a GP model for each objective %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Clean duplicates/invalid data and build training set
                [A1Dec,A1Obj,A1Con,A1]=DataUpdate(A1,A2);
                TSDec = A1Dec;
                TSObj = A1Obj; 
                TSCon = A1Con;
                maxnumData = 1000;%11*Problem.D-1+5;
                numTS = size(TSDec,1);
                if size(TSDec,1)>=maxnumData
                trainX = TSDec(numTS-maxnumData+1:end,:);
                trainY = TSObj(numTS-maxnumData+1:end,:);
                trainC = TSCon(numTS-maxnumData+1:end,:);
                else
                trainX = TSDec;
                trainY = TSObj;
                trainC = TSCon;
                end
                for i = 1 : Problem.M
                    % Objective GP
                    dmodel     = dacefit(trainX,trainY(:,i),'regpoly1','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                    Model{i}   = dmodel;
                    THETA(i,:) = dmodel.theta;
                end
                for i = 1 : numC
                    % Constraint GP
                    cdmodel     = dacefit(trainX,trainC(:,i),'regpoly1','corrgauss',CTHETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                    CModel{i}   = cdmodel;
                    CTHETA(i,:) = cdmodel.theta;
                end

                Stageflag=1;
                %% Stage1
                if Stageflag==1 &&  Problem.FE<=0.5*Problem.maxFE
                    % Surrogate-based NSGA-III on objectives, then infill by UPF
                    [PopDec,PopObj,MSE]=NSGAIII_opt(A1,Model,V0,Problem);
                    [A1,A2,New] = realFE_UPF(PopDec,PopObj,A1Dec,A1Obj,V0,Problem,Model,A2);
                end
                if iter == tr
                last_t = size(Arc_meanCV,1);
                mean_lastCV = mean(Arc_meanCV(floor(last_t/2):last_t,:));
                mean_preCV  = mean(Arc_meanCV(1:floor(last_t/2),:));
                deltaCV = (mean_lastCV-mean_preCV)/max(Arc_meanCV,[],1);
                %dCV = [dCV;deltaCV];
                if deltaCV<0
                    Stageflag=1;
                else
                    preCV = Arc_CV(1:floor(last_t/2),:);
                    %lastCV = Arc_CV(floor(last_t/2):end,:);
                    prefr = sum(preCV==0)/size(preCV,1);
                    %lastfr = sum(lastCV==0)/size(lastCV,1);
                    if prefr>0
                        Stageflag=1;
                    else
                        Stageflag=3;
                    end
                end
                t = t+1;
                end
                if Problem.FE>0.5*Problem.maxFE
                    if Stageflag==1
                        % Switch to constraint-focused stage after mid budget
                        Stageflag=2;
                    end
                end
                %% Stage2
                if  Stageflag==2 
                    % Focus on the most difficult constraint (lowest feasible rate)
                    [priority,~] = Constraint_priority(A1);
                    con_index = priority(1);
                    [PopDec,PopObj,MSE,PopCon,CMSE]=subCNSGAIII_opt(A1,Model,CModel,V01,Problem,numC,con_index);
                    [A1,A2] = realFE_CPF1(PopDec,PopObj,PopCon,A2,numC,Problem,con_index,A1Dec,A1Obj,A1Con,V0,Model,CModel);
                    if Problem.FE >= 0.7*Problem.maxFE
                        Stageflag = 3;
                    end
                elseif  Stageflag==3
                %% Stage3
                    % Full constrained NSGA-III with objective and constraint GPs
                    [PopDec,PopObj,~,PopCon,~]=CNSGAIII_opt(A1,A2,Model,CModel,V0,Problem,numC);
                    [A1,A2] = realFE_CPF(PopDec,PopObj,PopCon,A2,numC,Problem,A1Dec,A1Obj,A1Con,V0,Model,CModel);
                end
                if iter<=tr
                    if ~isempty(New)
                        % Track constraint violation trend for stage decision
                        CV = sum(max(New.cons,0),2);
                        currentCV = min(CV,[],1);
                        Arc_meanCV = [Arc_meanCV;currentCV];
                        Arc_CV = [Arc_CV;CV];
                    end
                    iter = iter+1;
                end
                X = sprintf('FE:%d;StageFlag:%d;',Problem.FE,Stageflag);
                disp(X);
            end
        end
    end
end
