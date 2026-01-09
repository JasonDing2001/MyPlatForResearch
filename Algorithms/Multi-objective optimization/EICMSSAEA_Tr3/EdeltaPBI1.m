function Popreal = EdeltaPBI1(PopDec,PopObj,TSDec,TSObj,W,Problem,Model)
         [FrontNo,~] = NDSort(TSObj,inf);
         ND_TSObj = TSObj(FrontNo==1,:);
         ND_TSDec = TSDec(FrontNo==1,:);
         N1 = size(ND_TSObj,1);
         
         [FrontNo1,~] = NDSort(PopObj,inf);
         ND_PopObj = PopObj(FrontNo1==1,:);
         ND_PopDec = PopDec(FrontNo1==1,:);
        %%    Angle  
         num=1;
         for i = 1:size(ND_PopObj,1)
            mu = ND_PopObj(i,:);
            R{i} =  repmat(mu,num,1);
         end
         allR =  cell2mat(R(1:end)');
         Angle  = acos(1-pdist2((max(allR,0)),max(ND_TSObj,0),'cosine'));

         [angle,~] = min(Angle,[],2);
         temp = reshape(angle,num,length(R));
         dd = mean(temp,1);
         dd = dd';
         dd = dd./max(dd);
         [~,dbest] = max(dd,[],1);
         Popreal_Dec1 = ND_PopDec(dbest,:);
         Popreal_Obj1 = ND_PopObj(dbest,:);
         
         ND_TSObj  = cat(1,ND_TSObj,Popreal_Obj1);
         ND_TSDec  = cat(1,ND_TSDec,Popreal_Dec1);
         ND_PopObj(dbest,:) = [];
         ND_PopDec(dbest,:) = [];
         %N2 = size(ND_PopObj,1);
         N1 = N1+1;
         Popreal = Popreal_Dec1;
         %% delta PBI
         if ~isempty(ND_PopDec)
             AllPopObj = [ND_TSObj;ND_PopObj];
             AllPopDec = [ND_TSDec;ND_PopDec];
             [z,znad]      = deal(min(AllPopObj),max(AllPopObj));
             %        
             [AllPopObj,~,~] = Normalization(AllPopObj,z,znad);
             ND_TSObj = AllPopObj(1:N1,:);
             %ND_TSDec = AllPopDec(1:N1,:);
             normT  = sqrt(sum(ND_TSObj.^2,2));
             Cosine = 1 - pdist2(ND_TSObj,W,'cosine');
             d1     = repmat(normT,1,size(W,1)).*Cosine;
             d2     = repmat(normT,1,size(W,1)).*sqrt(1-Cosine.^2);
             [~,class_T] = min(d2,[],2);% class_T:  reference index of each Archive pop corresponding 

             theta = 5;
             [~,ia,~] =  unique(class_T);
             unique_class_T = class_T(ia);
             Archive_PBI = -1*10000*ones(size(W,1),1);
             for i = unique_class_T'
                 % i is active reference index
                 % current is pop index of i reference 
                 current = find(class_T==i);
                 PBI = d1(current,i)+theta*d2(current,i); 
                 [bestPBI,~] = min(PBI);
                 Archive_PBI(i,1) = bestPBI;
             end

             activeW = W(unique_class_T,:);
             ND_PopObj = AllPopObj(N1+1:end,:);
             ND_PopDec = AllPopDec(N1+1:end,:);
             normP  = sqrt(sum(ND_PopObj.^2,2));
             Cosine = 1 - pdist2(ND_PopObj,activeW,'cosine');
             Pd1     = repmat(normP,1,size(activeW,1)).*Cosine;
             Pd2     = repmat(normP,1,size(activeW,1)).*sqrt(1-Cosine.^2);
             [~,class_P] = min(Pd2,[],2); %class_P:  reference index of each Candidate pop corresponding  

             [~,ia,~] =  unique(class_P);
             unique_class_C = class_P(ia);
             Candidate_PBI = 10000*ones(size(activeW,1),1);
             Archive_PBI_new = -1*10000*ones(size(W,1),1);
             Next = zeros(1,size(activeW,1));
             class_CC = [];
             for i = unique_class_C'
                 current = find(class_P==i);% 
                 PBI = Pd1(current,i)+theta*Pd2(current,i); % act index of current
                 [bestPBI,bestPBI_index] = min(PBI);
                 Candidate_PBI(i,1) = bestPBI;
                 class_C = unique_class_T(i);% unique_class_T(i)
                 Archive_PBI_new(class_C,1) = Candidate_PBI(i,1);
                 Next(i)  = current(bestPBI_index);
                 class_CC = cat(1,class_CC,class_C);
             end
             index = Next(Next~=0);
             PBI_old = Archive_PBI(class_CC,:);
             PBI_new = Archive_PBI_new(class_CC,:);
             delta = PBI_new-PBI_old;
             [~,best] = sort(delta);
             %num = min(length(best),batchsizes);
             Popreal_Dec2 = ND_PopDec(index(best(1)),:);
             Popreal_Obj2 = ND_PopObj(index(best(1)),:);

             ND_TSDec  = cat(1,ND_TSDec,Popreal_Dec2);
             ND_TSObj  = cat(1,ND_TSObj,Popreal_Obj2);

             ND_PopDec(index(best(1)),:) = [];
             ND_PopObj(index(best(1)),:) = [];
             Popreal = cat(1,Popreal_Dec1,Popreal_Dec2);
         end
         if ~isempty(PopDec)
             Popreal_Dec3 = noC_EPDM(TSDec,TSObj,PopDec,PopObj,Popreal,W,Problem,Model);
             Popreal = [Popreal;Popreal_Dec3];
         end
         

         

        
   
         
         
         

         
end


