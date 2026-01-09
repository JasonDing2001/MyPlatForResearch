function index = CKEnvironmentalSelection1(PopObj,PopCon,V,theta)
% The environmental selection of K-RVEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Cheng He
    [N,~] = size(PopObj);
    NV    = size(V,1);
    %% Translate the population
    CV = sum(max(0,PopCon),2);
    PopObj = [PopObj,CV];
    PopObj = PopObj - repmat(min(PopObj,[],1),N,1);
    %% Associate each solution to a reference vector
    Angle = acos(1-pdist2(PopObj,V,'cosine'));
    CosineAngle = 1-pdist2(PopObj,V,'cosine');
    SineAngle   = sqrt(1 - CosineAngle .^ 2);
    [~,associate] = min(Angle,[],2);
    
    %% Select one solution for each reference vector
    Next = zeros(1,NV);
    for i = unique(associate)'
        current1 = find(associate==i);
        if ~isempty(current1)
            % Calculate the APD value of each solution
            %APD = (1+M*theta*Angle(current1,i)/gamma(i)).*sqrt(sum(PopObj(current1,:).^2,2));
            PDM = theta*sqrt(sum(PopObj(current1,:).^2,2)).*SineAngle(current1,i)+mean(PopObj(current1,:),2);
            % Select the one with the minimum APD value
            [~,best] = min(PDM);
            Next(i)  = current1(best);
        end
    end
    % Population for next generation
    index = Next(Next~=0);
end