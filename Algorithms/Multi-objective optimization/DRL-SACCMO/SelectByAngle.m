function Popreal_Dec1 = SelectByAngle(ND_TSObj,ND_PopDec,ND_PopObj)
% Angle 最大化（多样性扩展）

    if isempty(ND_PopDec)
        Popreal_Dec1 = [];
        return;
    end
    if isempty(ND_TSObj)
        Popreal_Dec1 = ND_PopDec(randi(size(ND_PopDec,1)),:);
        return;
    end

    num = 1;
    for i = 1 : size(ND_PopObj,1)
        mu = ND_PopObj(i,:);
        R{i} = repmat(mu,num,1);
    end
    allR = cell2mat(R(1:end)');
    Angle = acos(1 - pdist2(max(allR,0),max(ND_TSObj,0),'cosine'));

    [angle,~] = min(Angle,[],2);
    temp = reshape(angle,num,length(R));
    dd = mean(temp,1);
    dd = dd';
    dd = dd./max(dd);
    [~,dbest] = max(dd,[],1);
    Popreal_Dec1 = ND_PopDec(dbest,:);
end
