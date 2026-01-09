function [TSObjC,PopObjC] = AppendNormalizedCV(TSObj,TSCon,PopObj,PopCon)
% 将 CV 归一化后拼接到目标空间（用于 C 子任务动作评分）

    CV_TS  = sum(max(0,TSCon),2);
    CV_Pop = sum(max(0,PopCon),2);

    CV_all = [CV_TS;CV_Pop];
    if isempty(CV_all)
        CV_TS_N  = zeros(size(CV_TS));
        CV_Pop_N = zeros(size(CV_Pop));
    else
        minCV = min(CV_all);
        maxCV = max(CV_all);
        if maxCV > minCV
            CV_TS_N  = (CV_TS - minCV)/(maxCV - minCV);
            CV_Pop_N = (CV_Pop - minCV)/(maxCV - minCV);
        else
            CV_TS_N  = zeros(size(CV_TS));
            CV_Pop_N = zeros(size(CV_Pop));
        end
    end

    TSObjC  = [TSObj,CV_TS_N];
    PopObjC = [PopObj,CV_Pop_N];
end
