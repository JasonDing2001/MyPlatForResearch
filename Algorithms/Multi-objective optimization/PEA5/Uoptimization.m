function [PopDec,PopObj,ObjMSE] = Uoptimization(A1,wmax,Model_obj,Problem)
    %% Optimization1 for unconstrained optimization
    PopDec = A1.decs;
    PopObj = A1.objs;
    ObjMSE = zeros(Problem.N,Problem.M);
    w      = 1;
    while w <= wmax
        drawnow();
        OffDec = OperatorGA(Problem,PopDec);
        [OffObj,Off_ObjMSE] = model_predict1(Model_obj,OffDec);

        PopDec = [PopDec;OffDec];
        PopObj = [PopObj;OffObj];
        ObjMSE = [ObjMSE;Off_ObjMSE];
        
        index  = SEnvironmentalSelection(PopObj,[],length(A1),1);
      
        PopDec = PopDec(index,:);
        PopObj = PopObj(index,:);
        ObjMSE = ObjMSE(index,:);
        w = w + 1;
    end
end

function [OffObj,Off_ObjMSE] = model_predict1(Model_obj,OffDec)
    N          = size(OffDec,1);
    Len_obj    = length(Model_obj);
    OffObj     = zeros(N,Len_obj);
    Off_ObjMSE = zeros(N,Len_obj);
      
    for i = 1 : N
        for j = 1 : Len_obj
            [OffObj(i,j),~,Off_ObjMSE(i,j)] = predictor(OffDec(i,:),Model_obj{j});
        end
    end
end