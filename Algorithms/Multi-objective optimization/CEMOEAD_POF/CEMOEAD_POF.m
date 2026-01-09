classdef CEMOEAD_POF < ALGORITHM
% <multi/many> <real/integer> <expensive> <constrained>
% Expected direction-based hypervolume improvement
% batch_size --- 5 --- number of true function evaluations per iteration

    methods 
        function main(Algorithm,Problem)
            %% Parameter setting
            batch_size = Algorithm.ParameterSet(5);
            n_init = 11*Problem.D-1;
            
            %% Generate initial design using LHS or other DOE methods
            x_lhs   = lhsdesign(n_init, Problem.D,'criterion','maximin','iterations',1000);
            x_init  = Problem.lower +  (Problem.upper - Problem.lower).*x_lhs;  
            Archive = Problem.Evaluation(x_init);     

            con_theta = repmat({(n_init ^ (-1 ./ n_init)) .* ones(1, Problem.D)}, 1, size(Archive.cons,2));
            obj_theta = repmat({(n_init ^ (-1 ./ n_init)) .* ones(1, Problem.D)}, 1, Problem.M);
            
            % find non-dominated solutions
            FrontNo = NDSort(Archive.objs,1); 


            %% Optimization
            while Algorithm.NotTerminated(Archive)
                Batch_size = min(batch_size, Problem.maxFE - Problem.FE);
                if Problem.FE <= 0.5 * Problem.maxFE
                stage = 1;
               % elseif Problem.FE <= 0.7 * Problem.maxFE
               % stage = 2;
                else
                stage = 3;
                end

                % Train surrogate models for objective function;
                [objModels ,obj_theta] = trainObj(obj_theta, Archive.decs, Archive.objs, Problem);
                if stage == 1
                    train_y = NormalizeObj(Archive.objs);
                    train_y_nds = train_y(FrontNo==1,:);
                    new_x = Opt_DirHV_EI(Problem.M,Problem.D,Problem.lower,Problem.upper,objModels,train_y_nds,Batch_size);  
                    
                elseif stage == 2

                elseif stage == 3
                    [conModels, con_theta] = trainCon(con_theta, Archive.decs , Archive.cons, Problem,size(Archive.cons,2));
                   
                    if isempty(Archive.best)
                        train_y_nds = NormalizeObj(Archive(FrontNo==1).objs);
                    else
                        train_y_nds = NormalizeObj(Archive.best.objs);
                    end
                    
                    new_x = Opt_DirHV_CEI(Problem.M,Problem.D,Problem.lower,Problem.upper,objModels,conModels,train_y_nds,Batch_size);

                end
                Archive = [Archive,Problem.Evaluation(new_x)];
                FrontNo = NDSort(Archive.objs,1);

            end
        end

end

end