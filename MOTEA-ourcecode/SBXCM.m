function cld_x = SBXCM(parent1,parent2,Xmin,Xmax,BI)
  d = BI.dim; % 获取决策变量的维度
  p_cross = 1;  % 交叉概率，设置为1，表示每次都会进行交叉
  epsilon = 1e-12;  % 用于判断两个父代是否相等的极小值
  eta_c = 30; % 交叉操作中的指数，控制交叉强度
  r = rand; % 生成一个随机数，用于判断是否进行交叉
if r <= p_cross    % 遍历每个决策变量（维度）    
        for i=1:d % 遍历每个决策变量（维度）
            if rand <= 0.5 % 50% 的概率决定是否进行交叉
                if abs( parent1(i)-parent2(i) ) > epsilon % 如果父代的差异足够大
                    if parent1(i) < parent2(i) % 如果第一个父代的值小于第二个父代的值
                        y1 = parent1(i);
                        y2 = parent2(i);
                    else % 如果第二个父代的值小于第一个父代的值
                        y1 = parent2(i);
                        y2 = parent1(i);
                    end
                    yl = Xmin(i); % 获取该维度的下界
                    yu = Xmax(i); % 获取该维度的上界
                    beta = 1.0 + (2.0*(y1-yl)/(y2-y1));  % 计算 beta 和 alpha，用于交叉操作的强度和范围
                    alpha = 2.0 - beta^(-(eta_c+1.0));
                    rand_var = rand;
                    if rand_var <= (1.0/alpha) % 生成一个随机数，用于决定交叉的方式 % 如果随机数小于等于 1/alpha，执行一个交叉方式
                        betaq = (rand_var*alpha)^(1.0/(eta_c+1.0));
                    else % 否则执行另一种交叉方式
                        betaq = (1.0/(2.0 - rand_var*alpha))^(1.0/(eta_c+1.0));
                    end
                    c1 = 0.5*((y1+y2) - betaq*(y2-y1)); % 计算第一个子代
                    beta = 1.0 + (2.0*(yu-y2)/(y2-y1)); % 计算第二个子代的 beta 和 alpha
                    alpha = 2.0 - beta^(-(eta_c+1.0));
                    rand_var = rand; % 生成新的随机数
                    if rand_var <= (1.0/alpha) % 决定第二个子代的交叉方式
                        betaq = (rand_var*alpha)^(1.0/(eta_c+1.0));
                    else
                        betaq = (1.0/(2.0 - rand_var*alpha))^(1.0/(eta_c+1.0));
                    end
                    c2 = 0.5*((y1+y2)+betaq*(y2-y1)); % 计算第二个子代


                     % 确保子代的值在合法范围内（下界和上界之间）
                    if (c1 < yl), c1 = yl; end
                    if (c2 < yl), c2 = yl; end
                    if (c1 > yu), c1 = yu; end
                    if (c2 > yu), c2 = yu; end

                     % 随机选择将 c1 和 c2 分配给子代
                    if rand <= 0.5
                        child1(i) = c1;
                        child2(i) = c2;
                    else
                        
                        child1(i) = c2;
                        child2(i) = c1;
                    end
                else   % 如果父代在该维度上的差异小于 epsilon，则不进行交叉，直接复制父代
                    child1(i) = parent1(i);
                    child2(i) = parent2(i);
                end
            else  % 如果不进行交叉，直接复制父代的值到子代
                child1(i) = parent1(i);
                child2(i) = parent2(i);
            end
        end
else   % 如果不进行交叉，直接复制父代
        child1 = parent1;
        child2 = parent2;
end

 % 对子代进行多项式变异
child1 = polymut(child1,Xmin,Xmax,BI);
child2 = polymut(child2,Xmin,Xmax,BI);

 % 最终选择一个子代作为结果
   cld_x = child1;
  if max(abs(parent1-cld_x))<min(1e-7,max(abs(parent1-child2)))  % 如果 child2 更接近 parent1，则选择 child2
    cld_x = child2;
  end  
%  if sum(abs(parent1-cld_x))==0
%     cld_x = polymut(cld_x,Xmin,Xmax,d);
%   end    
end
    