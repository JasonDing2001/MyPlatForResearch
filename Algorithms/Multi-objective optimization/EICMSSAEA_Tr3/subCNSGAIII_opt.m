function [PopDec,PopObj,MSE,PopCon,CMSE] = subCNSGAIII_opt(A1,Model,CModel,V01,Problem,numC,con_index)
          PopDec = A1.decs;
          PopObj = A1.objs;
          PopCon = A1.cons;
          PopCV  = sum(max(0,PopCon(:,con_index)),2);
          ALLObj = [PopObj,PopCV];
          
          Zmin   = min(ALLObj,[],1);
          if size(PopDec,1) >= Problem.N
             Next = nsga3EnvironmentalSelection(PopDec,ALLObj,Problem.N,V01,Zmin); 
             PopDec = PopDec(Next,:);
             ALLObj = ALLObj(Next,:);
          end
          w = 1;wmax = 20;
          beta = 0;cbeta = 0;
          while w <= wmax
                OffDec = generateOffing(PopDec,Problem,A1);
                PopDec = cat(1,PopDec,OffDec);  
                [PopObj,~,MSE,~] = GP_estimate(PopDec,Model,Problem.M,beta);
                [PopCon1,~,CMSE1,~] = CGP_estimate(PopDec,CModel,numC,cbeta);
                PopCon = PopCon1(:,con_index);
                CMSE = CMSE1(:,con_index);
                PopCV = sum(max(PopCon,0),2);
                ALLObj = [PopObj,PopCV];
                ALLMSE = [MSE,CMSE];
                Zmin  = min([Zmin; ALLObj],[],1);
                Choose = nsga3EnvironmentalSelection(PopDec,ALLObj,Problem.N,V01,Zmin); 
                PopDec = PopDec(Choose,:);
                ALLObj = ALLObj(Choose,:);
                ALLMSE = ALLMSE(Choose,:);
                w = w+1;
          end
          PopObj = ALLObj(:,1:Problem.M);
          PopCon = ALLObj(:,Problem.M+1:end);
          MSE    = ALLMSE(:,1:Problem.M);
          CMSE   = ALLMSE(:,Problem.M+1:end);

end

