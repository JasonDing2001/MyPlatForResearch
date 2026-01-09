function [PopDec,PopObj,PopCon,ObjMSE,ConMSE] = Coptimization(A1,wmax,Model_obj,Model_con,Problem,stage)

    PopDec = A1.decs;
    PopObj = A1.objs;
    PopCon = A1.cons;
    ObjMSE = zeros(Problem.N,Problem.M);
    ConMSE = zeros(Problem.N,size(PopCon,2));
    w      = 1;
    while w <= wmax
        drawnow();
        OffDec = OperatorGA(Problem,PopDec);
        [OffObj,Off_ObjMSE,OffCon,Off_ConMSE] = model_predict2(Model_obj,Model_con,OffDec);

        PopDec = [PopDec;OffDec];
        PopObj = [PopObj;OffObj];
        PopCon = [PopCon;OffCon];
        ObjMSE = [ObjMSE;Off_ObjMSE];
        ConMSE = [ConMSE;Off_ConMSE];
        
        index  = SEnvironmentalSelection(PopObj,PopCon,length(A1),stage);
      
        PopDec = PopDec(index,:);
        PopObj = PopObj(index,:);
        PopCon = PopCon(index,:);
        ObjMSE = ObjMSE(index,:);
        ConMSE = ConMSE(index,:);
        w = w + 1;
    end
end

function [OffObj,Off_ObjMSE,OffCon,Off_ConMSE] = model_predict2(Model_obj,Model_con,OffDec)
    N          = size(OffDec,1);
    Len_obj    = length(Model_obj);
    Len_con    = length(Model_con);
    OffObj     = zeros(N,Len_obj);
    OffCon     = zeros(N,Len_con);
    Off_ObjMSE = zeros(N,Len_obj);
    Off_ConMSE = zeros(N,Len_con);
        
    for i = 1 : N
        for j = 1 : Len_obj
            [OffObj(i,j),~,Off_ObjMSE(i,j)] = predictor(OffDec(i,:),Model_obj{j});
        end
        for j = 1 : Len_con
            [OffCon(i,j),~,Off_ConMSE(i,j)] = predictor(OffDec(i,:),Model_con{j});
        end
    end
end