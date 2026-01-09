classdef CDirHVEICV < ALGORITHM
% <multi/many> <real/integer> <expensive>
% Expected direction-based hypervolume improvement
% batch_size --- 5 --- number of true function evaluations per iteration

%------------------------------- Reference --------------------------------
% L. Zhao and Q. Zhang, Hypervolume-guided decomposition for parallel
% expensive multiobjective optimization, IEEE Transactions on Evolutionary
% Computation, 28(2): 432-444, 2024.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function was written by Liang Zhao.
% https://github.com/mobo-d/DirHV-EGO

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            batch_size = Algorithm.ParameterSet(5); 
            % number of initial samples
            n_init = 11*Problem.D-1;
            % Initial hyperparameters for GP
            theta = repmat({(n_init ^ (-1 ./ n_init)) .* ones(1, Problem.D)}, 1, Problem.M);
           
            
            %% Generate initial design using LHS or other DOE methods
            x_lhs   = lhsdesign(n_init, Problem.D,'criterion','maximin','iterations',1000);
            x_init  = Problem.lower +  (Problem.upper - Problem.lower).*x_lhs;  
            Archive = Problem.Evaluation(x_init);     
            theta_c = repmat({(n_init ^ (-1 ./ n_init)) .* ones(1, Problem.D)}, 1, size(Archive.cons,2));
           
            % find non-dominated solutions
            FrontNo = NDSort(Archive.objs,1); 

            %% Optimization
            while Algorithm.NotTerminated(Archive)
              %% Scale the objective values 
                train_x = Archive.decs; 
                % original objective values
                ori_objs = Archive.objs;

                ymin    = min(ori_objs,[],1); 
                ymax = max(ori_objs,[],1);
		            
                train_y = (ori_objs-ymin)./(ymax - ymin);
                tempobj = Archive.best.objs;
                
                if isempty(tempobj)
                  tempobj = Archive(FrontNo==1).objs;
                end
         
                train_y_nds = tempobj./(ymax - ymin);

                train_c = Archive.cons;
                train_CV = sum(train_c,2);
                %% Build GP model for each objective function 
                GPModels = cell(1,Problem.M);
                for i = 1 : Problem.M
                    GPModels{i} = Dacefit(train_x,train_y(:,i),'regpoly0','corrgauss',theta{i},1e-6*ones(1,Problem.D),20*ones(1,Problem.D));
                    theta{i}    = GPModels{i}.theta;
                end 

              %% Build GP model for constraints
               GPmodels_con = cell(1,size(train_CV,2)); 
               for i = 1 : size(train_CV,2)
                   GPmodels_con{i} = Dacefit(train_x,train_CV(:,i),'regpoly0','corrgauss',theta_c{i},1e-6*ones(1,Problem.D),20*ones(1,Problem.D));
                   theta_c{i}       = GPmodels_con{i}.theta;
               end

                
              %% Maximize DirHV-CEI using the MOEA/D framework and select multiple candidate points
                Batch_size = min(Problem.maxFE - Problem.FE,batch_size); % the total budget is  Problem.maxFE 
                new_x = Opt_DirHV_CEI(Problem.M,Problem.D,Problem.lower,Problem.upper,GPModels,GPmodels_con,train_y_nds,Batch_size);  
              
              %% Expensive Evaluation
                Archive = [Archive,Problem.Evaluation(new_x)];
                FrontNo = NDSort(Archive.objs,1);
            end
        end
    end
end