function [PopObj,z,znad] = NormalizationCC(PopObj,z,znad)
% RVEA 风格归一化

    [N,M] = size(PopObj);

    % 更新理想点
    z = min(z,min(PopObj,[],1));

    % 更新极端点/截距
    W = zeros(M) + 1e-6;
    W(logical(eye(M))) = 1;
    ASF = zeros(N,M);
    for i = 1 : M
        ASF(:,i) = max(abs((PopObj-repmat(z,N,1))./(repmat(znad-z,N,1)))./repmat(W(i,:),N,1),[],2);
    end
    [~,extreme] = min(ASF,[],1);
    Hyperplane = (PopObj(extreme,:)-repmat(z,M,1))\ones(M,1);
    a = (1./Hyperplane)' + z;
    if any(isnan(a)) || any(a<=z)
        a = max(PopObj,[],1);
    end
    znad = a;

    % 归一化
    PopObj = (PopObj-repmat(z,N,1))./(repmat(znad-z,N,1));
end
