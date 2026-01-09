function [ObjModels,obj_theta]= trainObj(obj_theta, x, y, Problem)
    % scale the objective values
    ymin = min(y,[],1);
    ymax = max(y,[],1);
    train_y = (y-ymin)./(ymax - ymin);
    
    % Build GP model for each objective function
    ObjModels = cell(1, Problem.M);
    for i = 1 : Problem.M
        ObjModels{i} = Dacefit(x, train_y(:,i), 'regpoly0', 'corrgauss', obj_theta{i}, 1e-6*ones(1, Problem.D), 20*ones(1, Problem.D));
        obj_theta{i} = ObjModels{i}.theta;
    end
    
end