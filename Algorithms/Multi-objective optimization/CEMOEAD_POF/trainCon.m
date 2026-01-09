function [conModels, con_theta] = trainCon(con_theta, x ,y, Problem,size)
    %scale
    ymin = min(y,[],1);
    ymax = max(y,[],1);
    y = (y - ymin) ./ (ymax - ymin);
    conModels = cell(1, size);
    for i = 1 : size
        conModels{i} = Dacefit(x, y(:,i), 'regpoly0', 'corrgauss', con_theta{i}, 1e-6*ones(1, Problem.D), 20*ones(1, Problem.D));
        con_theta{i} = conModels{i}.theta;
    end
end
