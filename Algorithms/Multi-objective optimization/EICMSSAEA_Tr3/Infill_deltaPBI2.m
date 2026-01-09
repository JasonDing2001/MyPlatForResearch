function Popreal = Infill_deltaPBI2(PopDec,PopObj,PopCon,TSDec,TSObj,TSCon,W,Model,CModel,Problem,numC)
         [FrontNo,~] = NDSort(TSObj,TSCon,inf);
         ND_TSDec = TSDec(FrontNo==1,:);
         ND_TSObj = TSObj(FrontNo==1,:);
         ND_TSCon = TSCon(FrontNo==1,:);
         
         % 
         [FrontNo1,~] = NDSort(PopObj,PopCon,inf);
         ND_PopDec = PopDec(FrontNo1==1,:);
         ND_PopObj = PopObj(FrontNo1==1,:);
         ND_PopCon = PopCon(FrontNo1==1,:);
         %
         % Angle-based selection
         [Popreal_Dec1, Popreal_Obj1, ND_PopDec, ND_PopObj] = SelectByAngle(ND_TSObj, ND_PopDec, ND_PopObj);
         
         ND_TSObj  = cat(1,ND_TSObj,Popreal_Obj1);
         ND_TSDec  = cat(1,ND_TSDec,Popreal_Dec1);
         Popreal = Popreal_Dec1;
         if ~isempty(ND_PopDec)
             Popreal_Dec2 = SelectByDeltaPBI(ND_TSObj, ND_TSDec, ND_PopObj, ND_PopDec, W);
             Popreal = cat(1,Popreal_Dec1,Popreal_Dec2);
         end
         Popreal_Dec3 = SelectByEPDM(TSDec,TSObj,TSCon,PopDec,PopObj,PopCon,Popreal,W,Model,CModel,Problem,numC);
         Popreal = [Popreal;Popreal_Dec3];

end


