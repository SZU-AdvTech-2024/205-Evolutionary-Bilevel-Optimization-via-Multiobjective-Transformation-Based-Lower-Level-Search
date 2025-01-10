function cld_x = polymut(cld_x,Xmin,Xmax,BI)
   d = BI.l_dim;  % 获取决策变量的维度
   pmut = 1/d; % 变异的概率设置为 1/d，每个决策变量变异的概率
   eta_m = 20; % 多项式变异的控制因子，控制变异的强度，值越大变异幅度越小
    for i = BI.u_dim+1:BI.u_dim+d % 遍历决策变量的每一个维度
        if rand <= pmut % 如果生成的随机数小于等于变异概率，则进行变异
            y = cld_x(i);  % 当前决策变量的值
            yl = Xmin(i);  % 当前决策变量的最小值（上下界的下界）
            yu = Xmax(i);  % 当前决策变量的最大值（上下界的上界）
            delta1 = (y-yl) / (yu-yl); % 当前值与最小值的归一化差值
            delta2 = (yu-y) / (yu-yl);  % 当前值与最大值的归一化差值
            rand_var = rand;   % 生成一个新的随机数，用于变异操作的概率控制
            mut_pow = 1.0/(eta_m+1.0); % 变异的幂指数，用来调节变异的幅度


            if rand_var <= 0.5 % 如果随机数小于等于 0.5，应用第一个变异公式
                xy = 1.0 - delta1;  % 计算变异系数的一个部分
                val = 2.0*rand_var + (1.0 - 2.0*rand_var) * xy^(eta_m+1.0);   % 计算变异值 适用于当 rand_var 较小（即小于等于 0.5）时，变异幅度较小，趋向于更接近最小
                deltaq =  val^mut_pow - 1.0; % 最终变异差值
            else % 如果随机数大于 0.5，应用第二个变异公式
                xy = 1.0 - delta2;  % 计算变异系数的另一个部分
                val = 2.0*(1.0 - rand_var) + 2.0*(rand_var-0.5) * xy^(eta_m+1.0);  % 计算变异值 适用于当 rand_var 较大（即大于 0.5）时，变异幅度较大，趋向于更接近最大值。
                deltaq = 1.0 - val^mut_pow;  % 最终变异差值
            end
            y = y + deltaq*(yu - yl); % 根据变异差值调整决策变量的值


            if (y<yl), y = yl; end  % 如果变异后的值小于最小值，则修正为最小值
            if (y>yu), y = yu; end   % 如果变异后的值大于最大值，则修正为最大值
            cld_x(i) = y;             % 将变异后的值赋回原数组
        end
    end