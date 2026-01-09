function Q_output = Cal_Q(Obj1, Obj2)
    % Cal_Q - 计算Q值矩阵，可以处理单个或两个解集的情况
    % 输入:
    %   Obj1 - 单个或第一个解集的目标值矩阵 (N1 x M)
    %   Obj2 - 第二个解集的目标值矩阵 (N2 x M)，可选
    % 输出:
    %   Q_output - 若输入一个解集，则返回单个解集的Q值向量
    %              若输入两个解集，则返回两个解集间的 Q 值矩阵 (N1 x N2)

    if nargin == 1
        N = size(Obj1,1);
        Obj1 = (Obj1-repmat(min(Obj1),N,1))./(repmat(max(Obj1)-min(Obj1),N,1));
        I = zeros(N);
        for i = 1 : N
            for j = 1 : N
                I(i,j) = max(Obj1(i,:)-Obj1(j,:));
            end
        end
        C = max(abs(I));
        F = sum(-exp(-I./repmat(C,N,1)/0.05)) + 1;
        Q_output = 1./F;

    elseif nargin == 2
        N1 = size(Obj1,1);
        N2 = size(Obj2,1);
        
        Obj = [Obj1;Obj2];
        N = size(Obj,1);
        Obj1 = (Obj1-repmat(min(Obj),N1,1))./(repmat(max(Obj)-min(Obj),N1,1));
        Obj2 = (Obj2-repmat(min(Obj),N2,1))./(repmat(max(Obj)-min(Obj),N2,1));


        I = zeros(N1,N2);
        for i = 1 : N1
            for j = 1 : N2
                I(i,j) = max(Obj1(i,:)-Obj2(j,:));
            end
        end
        C = max(abs(I));
        F = sum(-exp(-I./repmat(C,N1,1)/0.05)) + 1;
        Q_output = 1./F;
    end
end
