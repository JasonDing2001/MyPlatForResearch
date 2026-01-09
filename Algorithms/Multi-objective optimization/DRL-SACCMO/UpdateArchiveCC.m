function Archive = UpdateArchiveCC(Archive,N)
% 更新归档并控制规模

    [~,Unduplicated] = unique(Archive.objs,'rows');
    Archive = Archive(Unduplicated);

    if length(Archive) > N
        Fitness = CalFitnessCC(Archive.objs,Archive.cons);
        Next = Fitness < 1;
        if sum(Next) < N
            [~,Rank] = sort(Fitness);
            Next(Rank(1:N)) = true;
        elseif sum(Next) > N
            Del  = TruncationCC(Archive(Next).objs,sum(Next)-N);
            Temp = find(Next);
            Next(Temp(Del)) = false;
        end
        Archive = Archive(Next);
    end
end
