function X = repair_stayAtBoundary(X,xrange)
% 该函数用于修复矩阵 X 中的元素，使其保持在指定的范围内
    % 输入:
    %   X: 需要修复的矩阵，每行表示一个向量
    %   xrange: 每个维度的取值范围，由两行组成的矩阵
    %           xrange(1,:) 表示每个维度的最小值
    %           xrange(2,:) 表示每个维度的最大值
    % 输出:
    %   X: 修复后的矩阵，确保所有元素在指定范围内

    for i = 1 : size(X,1)
        % 检查矩阵 X 第 i 行中哪些元素小于对应维度的最小值
        flag = X(i,:) < xrange(1,:);
        % 对于这些元素，将其修复为最小值

        X(i,flag) = xrange(1,flag);

         % 检查矩阵 X 第 i 行中哪些元素大于对应维度的最大值
        flag = X(i,:) > xrange(2,:);


         % 对于这些元素，将其修复为最大值
        X(i,flag) = xrange(2,flag);
    end
end