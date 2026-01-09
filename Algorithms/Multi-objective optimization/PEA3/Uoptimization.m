function [PopDec,PopObj,ObjMSE] = Uoptimization(A1,wmax,Model_obj,Problem)

    PopDec = A1.decs;
    PopObj = A1.objs;
    ObjMSE = zeros(Problem.N,Problem.M);

    w      = 1;
    while w <= wmax
        drawnow();
        OffDec = OperatorGA(Problem,PopDec);
        [OffObj,Off_ObjMSE, ~, ~] = model_predict(Model_obj, [], OffDec, 0);

        PopDec = [PopDec;OffDec];
        PopObj = [PopObj;OffObj];
        ObjMSE = [ObjMSE;Off_ObjMSE];
        
        index  = USEnvironmentalSelection(PopObj,ObjMSE,length(A1));
      
        PopDec = PopDec(index,:);
        PopObj = PopObj(index,:);
        ObjMSE = ObjMSE(index,:);
        w = w + 1;
    end
end

function [OffObj,Off_ObjMSE,OffCon,Off_ConMSE] = model_predict(Model_obj,Model_con,OffDec,flag)
    N          = size(OffDec,1);
    Len_obj    = length(Model_obj);
    OffObj     = zeros(N,Len_obj);
    Off_ObjMSE = zeros(N,Len_obj);
    
    for i = 1 : N
        for j = 1 : Len_obj
            [OffObj(i,j),~,Off_ObjMSE(i,j)] = predictor(OffDec(i,:),Model_obj{j});
        end
    end

    OffCon = [];
    Off_ConMSE=[] ;

end
