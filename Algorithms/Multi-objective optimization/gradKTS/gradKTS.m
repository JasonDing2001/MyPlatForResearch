classdef gradKTS < ALGORITHM
% <multi/many> <real> <expensive> <constrained>
% alpha  ---  0.35 --- yuzhi
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% 1. Parameter setting
            phi1  = 0.1;
            wmax1 = 10;
            mu1   = 5;
            % [tau,phi,mu] 这些旧参数已经不再需要，因为我们换了决策核心
            % 但为了兼容可能的后续调用，保留变量定义
            [alpha] = Algorithm.ParameterSet(0.35);
            
            %% 2. Initialization
            p      = 1/Problem.M;
            CAsize = Problem.N;
            N      = Problem.N;
            P_latin = UniformPoint(N, Problem.D, 'Latin');
            Population = Problem.Evaluation(repmat(Problem.upper-Problem.lower,N,1).*P_latin+repmat(Problem.lower,N,1));
            A1     = Population;
            
            % Archives Initialization
            CA  = UpdateCA([],Population,CAsize);
            DA1 = UpdateDA(Population,[],Problem.N,p);
            P1  = Population;
            P2  = Population;
            DA  = DA1;
            
            % Initialize Model Parameters (Persistent Theta)
            THETA = 5.*ones((Problem.M + size(DA.cons,2)),Problem.D);
            Model = cell(1,(Problem.M + size(DA.cons,2)));
            
            
            %% 2. Main Loop
            while Algorithm.NotTerminated(A1)
                
                % =========================================================
                % Step 1: Mandatory Model Training (每一代都强制训练)
                % =========================================================
                
                % 1.1 Train Objectives (Using A1)
                Dec = A1.decs; Obj = A1.objs;
                for i = 1 : Problem.M
                    dmodel = dacefit(Dec,Obj(:,i),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                    Model{i} = dmodel; THETA(i,:) = dmodel.theta;
                end
                
                % 1.2 Train Constraints (Using DA+CA+P2)
                Dec = [DA.decs; CA.decs; P2.decs];
                Obj = [DA.cons; CA.cons; P2.cons];
                % 去重处理
                try
                    [~,index] = unique(roundn(Dec,-4),'rows');
                catch
                    [~,index] = unique(round(Dec*1e4)/1e4,'rows');
                end
                Dec = Dec(index,:); Obj = Obj(index,:);
                
                for i = Problem.M+1 : (Problem.M + size(DA.cons,2))
                    dmodel = dacefit(Dec,Obj(:,i-Problem.M),'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));
                    Model{i} = dmodel; THETA(i,:) = dmodel.theta;
                end
                
                % =========================================================
                % Step 2: Geometric Conflict Detection (你的创新点)
                % =========================================================
                
                % 2.1 筛选探针 (Probes): 当前种群 A1 中的 Rank 1 解
                [FrontNo, ~] = NDSort(A1.objs, 1); 
                
                % 提取 Rank 1 的掩码
                Rank1_Mask = (FrontNo == 1);
                
                % 获取决策变量
                Probes_Dec = A1(Rank1_Mask).decs;
                
                % [建议增加] 对探针进行去重 (Unique)
                % 因为 A1 中可能存在极其相近的解，或者纯目标Rank1包含重复解
                % 避免对同一个位置重复计算梯度，节省算力
                [~, unique_idx] = unique(round(Probes_Dec*1e6)/1e6, 'rows');
                Probes_Dec      = Probes_Dec(unique_idx, :);
                
                n_probes = size(Probes_Dec, 1);
                
               % 初始化拼接向量
                SuperVec_Obj = [];
                SuperVec_Con = [];
                epsilon      = 1e-8; 
                
                % [新增] 初始化统计计数器
                Count_Synergy  = 0; % 协同计数 (Cos > alpha)
                Count_Conflict = 0; % 冲突计数 (Cos < -alpha)
                Count_Neutral  = 0; % 中立/弱相关计数 (中间地带)
                
                % [新增] 记录每个探针的具体余弦值(可选，用于调试分析)
                Probe_Cosines  = zeros(n_probes, 1);
                
                % 2.2 遍历探针计算梯度并统计
                for k = 1 : n_probes
                    x_curr = Probes_Dec(k, :);
                    
                    % --- A. 计算综合目标下降方向 ---
                    sum_grad_obj = zeros(1, Problem.D);
                    for m = 1 : Problem.M
                        [~, g] = predictor(x_curr, Model{m});
                        g = g(:)';
                        sum_grad_obj = sum_grad_obj + (g / (norm(g) + epsilon));
                    end
                    vec_obj_unit = sum_grad_obj / (norm(sum_grad_obj) + epsilon);
                    
                    % --- B. 计算综合约束违反方向 (修正版：只看违反的约束) ---
                    sum_grad_con = zeros(1, Problem.D);
                    n_con_funcs = size(DA.cons, 2);
                    
                    Has_Violation = false;
                    
                    for c = 1 : n_con_funcs
                        % 同时获取预测值 y (即 g_c(x)) 和 梯度 g (即 ∇g_c(x))
                        [pred_val, g] = predictor(x_curr, Model{Problem.M + c});
                        
                        % [关键逻辑] 只有当约束被违反 (pred_val > 0) 时，才计入梯度
                        % 或者为了保险，设一个微小的阈值比如 -1e-5，把边界附近的也算上
                        if pred_val > 1e-6 
                            g = g(:)';
                            sum_grad_con = sum_grad_con + (g / (norm(g) + epsilon));
                            Has_Violation = true;
                        end
                    end
                    
                    % 如果所有约束都满足 (Has_Violation = false)
                    % 那么约束梯度理论上是 0，或者你可以认为任何方向都是“可行方向”
                    if ~Has_Violation
                        vec_con_unit = zeros(1, Problem.D); % 或者设为特定值
                    else
                        vec_con_unit = sum_grad_con / (norm(sum_grad_con) + epsilon);
                    end
                    
                    % --- [新增] C. 计算当前个体的余弦相似度 ---
                    % 因为 vec_obj_unit 和 vec_con_unit 已经是单位向量
                    % 所以点积结果直接就是 cos theta
                    cos_k = dot(vec_obj_unit, vec_con_unit);
                    Probe_Cosines(k) = cos_k; % 记录下来
                    
                    % --- [新增] D. 统计归类 ---
                    if cos_k > 0
                        Count_Synergy = Count_Synergy + 1;
                    elseif cos_k < 0
                        Count_Conflict = Count_Conflict + 1;
                    end
                    
                    % --- E. 拼接 (保留用于计算原来的全局Index，或者你可以选择弃用) ---
                    SuperVec_Obj = [SuperVec_Obj, vec_obj_unit];
                    SuperVec_Con = [SuperVec_Con, vec_con_unit];
                end
                
                % 计算全局 Index (保留原逻辑)
                Conflict_Index = dot(SuperVec_Obj, SuperVec_Con) / ...
                                 (norm(SuperVec_Obj) * norm(SuperVec_Con) + epsilon);
                
                % [新增] 打印统计信息 (调试用)
                fprintf('FE:%d | Total Probes:%d | Synergy:%d | Conflict:%d| Global Index:%.4f\n', ...
                        Problem.FE, n_probes, Count_Synergy, Count_Conflict, Conflict_Index);
                
                % 2.3 决策切换 (这里你可以选择用 Global Index，也可以改用上面的 Count 比例)
                % 推荐逻辑：既然统计了，可以用比例辅助决策，或者依旧用 Index + 随机
                
                if Conflict_Index > alpha
                    search_mode = 0; % 协同 -> 目标优先
                elseif Conflict_Index < -alpha
                    search_mode = 1; % 冲突 -> 约束优先
                else
                    % 中间地带随机
                    if rand < 0.5, search_mode = 0; else, search_mode = 1; end
                end
                
                % 2.3 计算几何冲突指标 (Cosine Similarity)
                % Index > 0: 协同 (Cooperative)
                % Index < 0: 冲突 (Conflicting)
                Conflict_Index = dot(SuperVec_Obj, SuperVec_Con) / ...
                                 (norm(SuperVec_Obj) * norm(SuperVec_Con) + epsilon);
                
                % 2.4 决策切换
                if Conflict_Index > alpha
                    search_mode = 0; % 目标优先 (Unconstrained/Cooperative)
                elseif Conflict_Index < -alpha
                    search_mode = 1; % 约束优先 (Constrained/Conflicting)
                else
                    if rand < 0.5
                        search_mode = 0;
                    else
                        search_mode = 1;
                    end
                end
                
                % 选择对应存档
                if search_mode == 0
                    DA = DA1;
                else
                    DA = P1;
                end
                
                % =========================================================
                % Step 3: Evolutionary Optimization (不变)
                % =========================================================
                
                CCA.obj = CA.objs; CCA.dec = CA.decs; CCA.con = CA.cons; CCA.MSE = zeros(size(CCA.con,1),Problem.M+size(CCA.con,2));
                CP2.obj = P2.objs; CP2.dec = P2.decs; CP2.con = P2.cons; CP2.MSE = zeros(size(CP2.con,1),Problem.M+size(CP2.con,2));
                CDA.obj = DA.objs; CDA.dec = DA.decs; CDA.con = DA.cons; CDA.MSE = zeros(size(CDA.con,1),Problem.M+size(CDA.con,2));
                
                w = 1;
                while w <= wmax1
                    if search_mode == 0
                        % Mode 0: KTA2 Strategy
                        [~,ParentCdec,~,ParentMdec] = MatingSelection_KTA2(CCA.obj,CCA.dec,CDA.obj,CDA.dec,Problem.N);
                        OffspringDec = [OperatorGA(Problem,ParentCdec,{1,20,0,0});OperatorGA(Problem,ParentMdec,{0,0,1,20})];
                    else
                        % Mode 1: KCCMO Strategy
                        Fitness1 = CalFitness(CP2.obj,CP2.con);
                        Fitness2 = CalFitness(CDA.obj);
                        MatingPool1 = TournamentSelection(2,Problem.N,Fitness1);
                        MatingPool2 = TournamentSelection(2,Problem.N,Fitness2);
                        OffspringDec = [OperatorGA(Problem,CP2.dec(MatingPool1,:)); OperatorGA(Problem,CDA.dec(MatingPool2,:))];
                    end
                    
                    % Predictor (using updated models)
                    Pop.dec = OffspringDec;
                    N_off   = size(Pop.dec,1);
                    PopObj  = zeros(N_off, Problem.M + size(DA.cons,2));
                    Pop.MSE = zeros(N_off, Problem.M + size(DA.cons,2));
                    
                    for i = 1 : N_off
                        for j = 1 : (Problem.M + size(DA.cons,2))
                            [PopObj(i,j),~,Pop.MSE(i,j)] = predictor(Pop.dec(i,:), Model{j});
                        end
                    end
                    
                    Pop.obj = PopObj(:,1:Problem.M);
                    Pop.con = PopObj(:,Problem.M+1:end);
                    
                    PopC = cat_struct(CCA,Pop); CCA = K_UpdateCA(PopC,CAsize);
                    PopD = cat_struct(CDA,Pop);
                    
                    if search_mode == 0
                        CDA = K_UpdateDA(PopD,Problem.N,p);
                    else
                        CDA = K_UpdateP(PopD,Problem.N,false);
                    end
                    
                    PopC1 = cat_struct(CP2,Pop); [CP2,~] = K_UpdateP(PopC1,Problem.N,true);
                    w = w + 1;
                end
                
                % =========================================================
                % Step 4: Infill Sampling (不变)
                % =========================================================
                if search_mode == 0
                    [~,ia,~] = setxor(CCA.dec,A1.decs,'rows'); CCA = givevalue(CCA,ia);
                    [~,ia,~] = setxor(CDA.dec,A1.decs,'rows'); CDA = givevalue(CDA,ia);
                    Offspring01 = Adaptive_sampling(CCA.obj,CDA.obj,CCA.dec,CDA.dec,CDA.MSE,DA,P2,mu1,p,phi1);
                else
                    [CCA2,~]    = KCCMO_sampling(CP2,P2,mu1);
                    Offspring01 = CCA2.dec;
                end
                
                [~,index] = unique(Offspring01 ,'rows');
                PopNew = Offspring01(index,:);
                
                if ~isempty(PopNew)
                    Offspring = Problem.Evaluation(PopNew);
                    
                    temp = A1.decs;
                    for i = 1 : size(Offspring,2)
                        dist2 = pdist2(Offspring(i).decs,temp);
                        if min(dist2) > 1e-5
                            A1 = [A1,Offspring(i)];
                        end
                        temp = A1.decs;
                    end
                    
                    CA     = UpdateCA(CA,Offspring,CAsize);
                    DA1    = UpdateDA(DA1,Offspring,Problem.N,p);
                    [P1,~] = Update_P([P1,Offspring],Problem.N,false);
                    [P2,~] = Update_P([P2,Offspring],Problem.N,true);
                    
                    % 打印日志，观察你的Idea是否生效
                    fprintf('FE: %d | Mode: %d | Conflict_Index: %.4f\n', Problem.FE, search_mode, Conflict_Index);
                end
            end
        end
    end
end