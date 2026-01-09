function ArchiveC = SelectConstrainedArchive(Archive,N)
% 从归档中构造 C 子任务的可行/低违约子集

    if isempty(Archive)
        ArchiveC = Archive;
        return;
    end

    CV = sum(max(0,Archive.cons),2);
    feasible = CV == 0;
    if any(feasible)
        ArchiveC = Archive(feasible);
        if length(ArchiveC) > N
            ArchiveC = UpdateArchiveCC(ArchiveC,N);
        end
    else
        [~,idx] = sort(CV);
        idx = idx(1:min(N,length(idx)));
        ArchiveC = Archive(idx);
    end
end
