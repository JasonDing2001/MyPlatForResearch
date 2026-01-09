function [A1,A2] = realFE_CPF1(PopDec,PopObj,PopCon,A2,numC,Problem,con_index,A1Dec,A1Obj,A1Con,V0,Model,CModel)
                
                PopNew = Infill_deltaPBI2_sub(PopDec,PopObj,PopCon,A1Dec,A1Obj,A1Con,V0,Model,CModel,Problem,numC,con_index);
                if ~isempty(PopNew)
                    [~,ib]= intersect(PopNew,A2.decs,'rows');
                    PopNew(ib,:)=[];

                    if ~isempty(PopNew)
                        New       = Problem.Evaluation(PopNew);
                    else
                        New = [];
                    end
                else
                    New = [];
                end
                A2        = cat(2,A2,New);
                A1 = A2;
end

