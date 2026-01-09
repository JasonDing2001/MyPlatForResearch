function Popreal_Dec3 = SelectByEPDM(TSDec,TSObj,TSCon,PopDec,PopObj,PopCon,Popreal,W,Model,CModel,Problem,numC)
% Select one point by EPDM criterion

    Popreal_Dec3 = [];
    if ~isempty(PopDec)
        Popreal_Dec3 = EPDM(TSDec,TSObj,TSCon,PopDec,PopObj,PopCon,Popreal,W,Model,CModel,Problem,numC);
    end
end

