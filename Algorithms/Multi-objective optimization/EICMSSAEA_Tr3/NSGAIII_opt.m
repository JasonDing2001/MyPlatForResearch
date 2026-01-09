function [PopDec,PopObj,MSE] = NSGAIII_opt(A1,Model,V0,Problem)
          PopDec = A1.decs;
          PopObj = A1.objs;
          Zmin   = min(PopObj,[],1);
          if size(PopDec,1) >= Problem.N
             Next = nsga3EnvironmentalSelection(PopDec,PopObj,Problem.N,V0,Zmin); 
             PopDec = PopDec(Next,:);
             PopObj = PopObj(Next,:);
          end
          w = 1;wmax = 20;
          beta = 0;
          while w <= wmax
                OffDec = generateOffing(PopDec,Problem,A1);
                PopDec = cat(1,PopDec,OffDec);  
                [PopObj,~,MSE,~] = GP_estimate(PopDec,Model,Problem.M,beta);
                Zmin  = min([Zmin; PopObj],[],1);
                Choose = nsga3EnvironmentalSelection(PopDec,PopObj,Problem.N,V0,Zmin); 
                PopDec = PopDec(Choose,:);
                PopObj = PopObj(Choose,:);
                MSE = MSE(Choose,:);
                w = w+1;
          end
end

