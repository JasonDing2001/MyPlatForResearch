classdef Experience_Replay < handle
    % 经验回放：循环缓冲区
    properties
        buffer;     % 存储经验的矩阵
        index = 1;  % 当前插入位置
        overflow = 0; % 是否已填满
        N;          % 最大容量
        nS;         % 状态维度
    end

    methods
        function obj = Experience_Replay(n_S, maxER)
            % 构造函数：初始化缓冲区
            obj.N = maxER;
            obj.buffer = zeros(maxER, (2*n_S + 2));
            obj.nS = n_S;
        end

        function insert_experience(obj,S,A,R,Snew)
            % 插入一条经验
            if obj.index > obj.N
                obj.index = 1;
                obj.overflow = 1;
            end
            nS = obj.nS;
            insert = zeros(1,(2*nS + 2));
            for i = 1:(2*nS + 2)
                if i < nS+1
                    insert(i) = S(i);
                elseif i ==  nS + 1
                    insert(i) = A;
                elseif i == nS + 2
                    insert(i) = R;
                elseif i > (nS + 2)
                    insert(i) = Snew(i - nS - 2);
                end
            end
            obj.buffer(obj.index, 1:(2*nS + 2)) = insert;
            obj.index = obj.index + 1;
        end

        function [Sold, A, R, Snew] = get_batch(obj, batch_size)
            % 随机采样一个批次
            if ~obj.overflow
                rand_index = obj.index - 1;
            else
                rand_index = obj.N;
            end
            all_i = randi(rand_index, [batch_size, 1]);
            batch = obj.buffer(all_i, 1:(2*obj.nS + 2));
            Sold = batch(1:batch_size, 1:obj.nS);
            A = batch(1:batch_size, obj.nS+1);
            R = batch(1:batch_size, obj.nS+2);
            Snew = batch(1:batch_size, obj.nS+3:2*obj.nS+2);
        end
    end
end
