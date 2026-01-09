function [Model_obj,Model_con,THETA_obj,THETA_con] = model_train(A2,THETA_obj,THETA_con,stage)
   try
% Train the surrogate models for objectives
    Dec = A2.decs;
    Obj = A2.objs;
    Len_dec = size(Dec,2);
    Len_obj = size(Obj,2);
        for i = 1 : Len_obj
        [~,distinct1] = unique(round(Dec*1e100)/1e100,'rows');
        [~,distinct2] = unique(round(Obj(:,i)*1e100)/1e100,'rows');
        distinct = intersect(distinct1,distinct2);
        
        dmodel     = dacefit(Dec(distinct,:),Obj(distinct,i),'regpoly1','corrgauss',THETA_obj(i,:),1e-5.*ones(1,Len_dec),100.*ones(1,Len_dec));
        Model_obj{i}   = dmodel;
        THETA_obj(i,:) = dmodel.theta;
    end

    % Train the surrogate models for constraints
    if stage ~= 1
        Con = A2.cons;
        Len_con = size(Con,2);
        for i = 1 : Len_con
            [~,distinct1] = unique(round(Dec*1e100)/1e100,'rows');
            [~,distinct2] = unique(round(Con(:,i)*1e100)/1e100,'rows');
            distinct = intersect(distinct1,distinct2);
            
            dmodel     = dacefit(Dec(distinct,:),Con(distinct,i),'regpoly1','corrgauss',THETA_con(i,:),1e-5.*ones(1,Len_dec),100.*ones(1,Len_dec));
            Model_con{i}   = dmodel;
            THETA_con(i,:) = dmodel.theta;
        end
    else
        % If stage == 1, then Model_con is empty
        Model_con = [];
    end
   catch 
       keyboard;
   end

end