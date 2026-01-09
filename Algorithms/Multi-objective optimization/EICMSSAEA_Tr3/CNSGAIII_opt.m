function [PopDec,PopObj,MSE,PopCon,CMSE] = CNSGAIII_opt(A1,A2,Model,CModel,V0,Problem,numC)
          PopDec = A1.decs;
          PopObj = A1.objs;
          PopCon = A1.cons;
          %PopCV  = sum(max(0,PopCon(:,con_index)),2);
          Zmin   = min(PopObj,[],1);
          if size(PopDec,1) >= Problem.N
             Next = cnsga3EnvironmentalSelection(PopDec,PopObj,PopCon,Problem.N,V0,Zmin); 
             PopDec = PopDec(Next,:);
%              PopObj = PopObj(Next,:);
%              PopCon = PopCon(Next,:);
          end
          w = 1;wmax = 20;
          beta = 0;cbeta = 0;
          while w <= wmax
                OffDec = generateOffing(PopDec,Problem,A1);
                PopDec = cat(1,PopDec,OffDec);  
                [PopObj,~,MSE,~] = GP_estimate(PopDec,Model,Problem.M,beta);
                [PopCon,~,CMSE,~] = CGP_estimate(PopDec,CModel,numC,cbeta);
                Zmin  = min([Zmin; PopObj],[],1);
                 % Nodominated sorting considering the uncertainty           
                %Choose = SelectionMSE(ALLObj,ALLMSE,RealFirstObj,Problem.N,V01,Zmin);  
                Choose = cnsga3EnvironmentalSelection(PopDec,PopObj,PopCon,Problem.N,V0,Zmin); 
                PopDec = PopDec(Choose,:);
                PopObj = PopObj(Choose,:);
                PopCon = PopCon(Choose,:);
                MSE    = MSE(Choose,:);
                CMSE   = CMSE(Choose,:);
                w = w+1;
          end


end

