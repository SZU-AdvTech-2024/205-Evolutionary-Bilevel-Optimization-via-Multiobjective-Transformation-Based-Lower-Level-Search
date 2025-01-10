function W  = generateweight(UX,BI)   
%      UX = reshape(UX,num,BI.u_dim);
%      num = size(UX,1);
%     cofficentMatrix = tril(ones(M,M),-1);
%     X = repmat(x(1:M),M,1);
%     X = X.*cofficentMatrix;
%     for i = 1:M-1
%         X(i,i) = 1-x(i);
%     end
%     X(M,M) = 1;
%     X = X+(tril(ones(M,M),-1))';
%     y = cumprod(X,2);
% 每个个体生成一个权重向量。其基本思想是利用个体的属性（UX）
% 来计算一个基于正弦和余弦函数的加权矩阵，
% 然后通过累积乘积（cumprod）来得到最终的权重。这些权重反映了每个个体的相对重要性或者影响力





num = size(UX,2)/BI.u_dim;   % 计算个体的数量 num


for j = 1 : num      % 对每个个体进行操作

     % 规范化 UX 中的每个个体的属性值（x） 将每个个体的属性值（UX）进行 归一化
    x = [(UX((j-1)*BI.u_dim+1:j*BI.u_dim)-BI.ux_lb)./max([ones(1,BI.u_dim)*1e-6;BI.ux_ub-BI.ux_lb]),0];


    if any(x>1) || any(x<0) 
        disp('error'); % 如果 x 超出 [0, 1] 范围，则输出错误信息
    end
%     y = x/sqrt(sum(x.^2));
    M = BI.u_dim+1;  % 设置 M 的值，表示矩阵的维度


    cofficentMatrix = tril(ones(M,M),-1); % 生成下三角矩阵（不包括对角线）
    
     % 计算 X 矩阵，其中每个元素是基于 x 中的值来生成的
    X =cos( pi/2*repmat(x(1:M),M,1)); % 每个元素由余弦函数计算得出
        X = X.*cofficentMatrix; % 将 X 与下三角矩阵进行逐元素相乘

        % 修改 X 中的对角线元素
        for i = 1:M-1
            X(i,i) = sin(pi/2*x(i));  % 对角线元素由 x 的值决定
        end 
        X(M,M) = 1; % 最后一个对角线元素设为 1
        

           % 将 X 与下三角矩阵的转置矩阵相加
        X = X+(tril(ones(M,M),-1))';

          % 将 y 的最后一列作为权重
        y = cumprod(X,2);  % 每行元素的累积乘积


%     W(:,j) = y';   
     W(:,j) =y(:,end);  % 存储最后一列的值作为个体 j 的权重
end 
end