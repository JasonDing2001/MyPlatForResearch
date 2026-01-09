function numRV = NoactivateRV(A2,V)
           
           PopObj = A2.objs;
           PopCon = A2.cons;
           [FrontNo,~] = NDSort(PopObj,PopCon,inf);
           PopObj = PopObj(FrontNo==1,:);
           [N,~] = size(PopObj);
          %% Translate the population
           PopObj = PopObj - repmat(min(PopObj,[],1),N,1);
          %% Calculate the smallest angle value between each vector and others
           cosine = 1 - pdist2(V,V,'cosine');
           cosine(logical(eye(length(cosine)))) = 0;
           Angle = acos(1-pdist2(PopObj,V,'cosine'));
           [~,associate] = min(Angle,[],2);
           numRV = length(unique(associate));
end

