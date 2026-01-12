classdef CCMO_DRL_Kriging < ALGORITHM
% <multi> <real> <expensive> <constrained>
% DRL-SACCMO：基于协同结构与深度强化学习的昂贵约束多目标优化
% wmax   --- 20 --- 每次真实评估前的代理演化代数
% mu     ---  5 --- 每轮真实评估补点数量

    methods
        function main(Algorithm, Problem)
            %% 参数设置
            [wmax,mu] = Algorithm.ParameterSet(20,5);

            %% 初始真实评估样本
            NI      = 11*Problem.D - 1;
            P       = lhsamp(NI,Problem.D);
            PopInit = Problem.Evaluation(repmat(Problem.upper-Problem.lower,NI,1).*P + repmat(Problem.lower,NI,1));

            % ArchiveAll 保存全部真实评估样本，Archive 用于终止与输出
            ArchiveAll = PopInit;
            Archive    = UpdateArchiveCC(PopInit,Problem.N);

            % Kriging 超参数初始化
            numC      = size(PopInit.cons,2);
            THETA_OBJ = 5.*ones(Problem.M,Problem.D);
            if numC > 0
                THETA_CON = 5.*ones(numC,Problem.D);
            else
                THETA_CON = [];
            end

            % 参考向量（PBI 动作使用）
            [W,~]   = UniformPoint(Problem.N,Problem.M);

            %% U/C 子任务归档
            ArchiveU = Archive;
            ArchiveC = SelectConstrainedArchive(Archive,Problem.N);

            %% 强化学习代理初始化
            LastArchiveU = ArchiveU;
            LastArchiveC = ArchiveC;
            [stateU,num_actions,~,~,~] = GenerateState(Problem,-1,LastArchiveU,ArchiveU);
            [stateC,~,~,~,~]           = GenerateState(Problem,-1,LastArchiveC,ArchiveC);
            num_states = length(stateU);

            exp_replay_freq   = 8;
            copy_weights_freq = exp_replay_freq + 3;
            batch_size        = 16;
            stepU             = 0;
            stepC             = 0;

            agentU = DDQN(num_states,num_actions,[5,5],100);
            agentC = DDQN(num_states,num_actions,[5,5],100);

            %% 外层真实评估迭代
            while Algorithm.NotTerminated(ArchiveAll)
                % 两个代理分别选择动作
                actionU = agentU.action(stateU);
                actionC = agentC.action(stateC);
                stepU   = stepU + 1;
                stepC   = stepC + 1;

                % 训练 Kriging 模型（目标与约束分开建模）
                [Model,CModel,THETA_OBJ,THETA_CON] = TrainKrigingModels(ArchiveAll,Problem,THETA_OBJ,THETA_CON);

                % 内层代理演化（CCMO 协同结构）
                [PopU,PopC] = SurrogateCoevolution(Archive,Problem,Model,CModel,wmax);

                % 依据动作选择真实评估补点
                muU     = ceil(mu/2);
                muC     = mu - muU;
                NewDecU = SelectInfillPoints(actionU,PopU,ArchiveAll,W,Model,CModel,Problem,muU,false);
                NewDecC = SelectInfillPoints(actionC,PopC,ArchiveAll,W,Model,CModel,Problem,muC,true);
                NewDec  = unique([NewDecU;NewDecC],'rows');

                % 兜底补点，保证每轮数量
                if size(NewDec,1) < mu
                    NewDec = FillWithCandidates(NewDec,PopU.Dec,PopC.Dec,mu,Problem);
                end
                if isempty(NewDec)
                    % TODO: 候选池为空时可加入更稳健的备用采样策略
                    break;
                end
                if size(NewDec,1) > mu
                    NewDec = NewDec(1:mu,:);
                end

                % 真实评估与归档更新
                New        = Problem.Evaluation(NewDec);
                ArchiveAll = [ArchiveAll,New];
                Archive    = UpdateArchiveCC([Archive,New],Problem.N);

                % 更新 U/C 子任务归档
                ArchiveU = Archive;
                ArchiveC = SelectConstrainedArchive(Archive,Problem.N);

                % 更新 U 代理经验
                [stateU,actionU,next_stateU,~,rewardU] = GenerateState(Problem,actionU,LastArchiveU,ArchiveU);
                agentU.store(stateU,actionU,rewardU,next_stateU);
                stateU       = next_stateU;
                LastArchiveU = ArchiveU;

                if mod(stepU,exp_replay_freq) == 0
                    agentU.experience_replay(batch_size);
                end
                if mod(stepU,copy_weights_freq) == 0
                    agentU.copy_weights_agent_to_target();
                end

                % 更新 C 代理经验
                [stateC,actionC,next_stateC,~,rewardC] = GenerateState(Problem,actionC,LastArchiveC,ArchiveC);
                agentC.store(stateC,actionC,rewardC,next_stateC);
                stateC       = next_stateC;
                LastArchiveC = ArchiveC;

                if mod(stepC,exp_replay_freq) == 0
                    agentC.experience_replay(batch_size);
                end
                if mod(stepC,copy_weights_freq) == 0
                    agentC.copy_weights_agent_to_target();
                end
            end
        end
    end
end
