function OffDec = generateOffing(PopDec,Problem,A1)
                    if size(PopDec,1)<Problem.N
                        MatingPool = randi(size(PopDec,1),1,Problem.N);
                        OffDec = OperatorGA(Problem,PopDec(MatingPool,:));
                    else
                        OffDec = OperatorGA(Problem,PopDec);
                    end
                    pop_candi = [];
                    NP = size(OffDec,1);
                    for ii = 1 : NP
                        if min(sqrt(sum((OffDec(ii,:) - [A1.decs;pop_candi]).^2,2)))>1E-6
                            pop_candi = cat(1,pop_candi,OffDec(ii,:));
                        end
                    end
                    OffDec = pop_candi;
end

