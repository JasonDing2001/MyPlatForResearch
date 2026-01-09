function [A1Dec,A1Obj,A1Con,A1] = DataUpdate(A1,A2)
                A1Dec = A1.decs;
                A1Obj = A1.objs;
                A1Con = A1.cons;
                [c,ia,ic] = unique(A1Obj,'rows','stable');
                A1Obj = A1Obj(ia,:);
                A1Dec = A1Dec(ia,:);
                A1Con = A1Con(ia,:);
                A1 = A1(ia);

                A1(find(sum(isnan(A1Obj),2)>0))=[];
                A1Dec(find(sum(isnan(A1Obj),2)>0),:)=[];
                A1Obj(find(sum(isnan(A1Obj),2)>0),:)=[];
                A1Con(find(sum(isnan(A1Obj),2)>0),:)=[];

                A1(find(sum(isinf(A1Obj),2)>0))=[];
                A1Dec(find(sum(isinf(A1Obj),2)>0),:)=[];
                A1Obj(find(sum(isinf(A1Obj),2)>0),:)=[];
                A1Con(find(sum(isinf(A1Obj),2)>0),:)=[];

                [c,ia,ic] = unique(A1Obj,'rows');
                A1Obj = A1Obj(ia,:);
                A1Dec = A1Dec(ia,:);
                A1Con = A1Con(ia,:);
                A1 = A1(ia);
                [c,ib,ic] = unique(A1Dec,'rows');
                A1Obj = A1Obj(ib,:);
                A1Dec = A1Dec(ib,:);
                A1Con = A1Con(ib,:);
                A1 = A1(ib);
                
                
%                 tr_x = A2.decs;
%                 tr_y = A2.objs;
%                 [~,distinct1] = unique(round(A2.decs*1e6)/1e6,'rows');
%                 [~,distinct2] = unique(round(A2.objs*1e6)/1e6);
%                 distinct = intersect(distinct1,distinct2);
%                 tr_xx     = tr_x(distinct,:);
%                 tr_yy     = tr_y(distinct,:);
end

