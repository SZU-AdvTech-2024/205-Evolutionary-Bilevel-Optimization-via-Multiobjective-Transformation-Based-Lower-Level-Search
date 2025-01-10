% @Author: hxy
% @Date:   2017-10-16 18:09:13
% @Last Modified by:   hxy
% @Last Modified time: 2018-03-08 11:56:29
function BI = getBLOPinfo(benchmark,fno,dim)
% 函数定义，返回一个结构体 BI，包含基准测试问题的相关信息。
% benchmark: 基准函数名，确定使用哪一类基准问题。 smd or tp
% fno: 目标函数编号，用于选择特定的目标函数。如tp1-1 tp2-》
% dim: 问题的维度。
% 设置并返回最终的结构体 BI，包含所有信息

%返回的信息如下
% BI.u_ieqcon_num = ul_ieqcon_num; % 设置上层不等式约束的数量
% BI.l_ieqcon_num = ll_ieqcon_num;
% BI.u_eqcon_num = ul_eqcon_num;
% BI.l_eqcon_num = ll_eqcon_num; % 设置下层等式约束的数量
% 
% 
% BI.u_dim = ulDim; % 设置上层的决策变量维度
% BI.l_dim = llDim;
% BI.dim = ulDim + llDim; % 设置总的决策变量维度（上层 + 下层）
% 
% BI.u_lb = ulDimMin; % 设置上层决策变量的下界（最小值）
% BI.u_ub = ulDimMax;
% BI.l_lb = llDimMin;
% BI.l_ub = llDimMax; % 设置下层决策变量的上界（最大值）
% 
% 
% 
% % 设置上下层决策变量的范围，合并了上层和下层的最小值与最大值
% BI.xrange = [ulDimMin llDimMin;ulDimMax llDimMax];
% 
% 
% % 设置容忍度，用于停止准则判定
% BI.u_ftol = 1e-6; % 上层目标函数的容忍度
% BI.l_ftol = 1e-6;
% 
% % from the original BLEAQ2 source code  % 设置停止准则
% BI.ulStoppingCriteria = ulStoppingCriteria;
% BI.llStoppingCriteria = llStoppingCriteria;
% 
% 
% % 设置最优函数值，通常是在基准问题中预设的最优值
% BI.u_fopt = ulOpt;
% BI.l_fopt = llOpt;
% 
% BI.fn = fn;
% 
% % from BLEAQ2% 设置上下层的种群大小
% BI.u_N = ulPopSize;
% BI.l_N = llPopSize;
% BI.u_maxGen = 2000;% 设置最大迭代次数
% BI.l_maxGen = 2000;
% % 设置上下层是否包含下层约束（用于优化过程的约束处理）
% BI.isLowerLevelConstraintsIncludedInUpperLevel = isLowerLevelConstraintsIncludedInUpperLevel;

% BI.UmaxFEs = 2000;
% BI.LmaxFEs = 500;
% BI.UmaxImprFEs = 400;
% BI.LmaxImprFEs = 100;


% 如果输入参数大于1，处理多目标函数的情况
if nargin > 1
	fno = fno(:);% 将目标函数编号转换为列向量
	if length(fno)>1 
		BI = []; % 初始化输出结构体
		for i = fno'
	        if nargin == 3
                % 如果传入了维度参数，递归调用，合并多个目标函数信息
	            BI = [BI;getBLOPinfo(benchmark,i,dim)];
            else
                % 如果没有传入维度参数，仍然递归调用
	            BI = [BI;getBLOPinfo(benchmark,i)];
	        end
		end
		return;% 返回合并的结构体
	end
end

% 根据基准函数的名称进行不同的处理
if strcmp(benchmark,'TP') % 如果基准问题是 'TP'
	fn = sprintf('tp%d',fno); % 构造目标函数的名称
	assert(ismember(fno,1:10));  % 断言目标函数编号在1到10之间
	ulStoppingCriteria = 1e-4; % 上层停止准则
	llStoppingCriteria = 1e-5; % 下层停止准则

	%% check where the constraints are from

     % 根据目标函数编号，决定约束条件的来源
	switch fno
		case {2,3,4,5,6,7,8}

            % 当两个层次具有相同约束时，去除重复的约束
			% this happens when the two levels have the same constraints
			% thus the duplications should be removed
			isLowerLevelConstraintsIncludedInUpperLevel = true;
        otherwise

            % 否则，认为下层约束不包含在上层
			isLowerLevelConstraintsIncludedInUpperLevel = false;
	end
	%% population size


        % 设置种群大小
	if fno == 10
	    ulPopSize=200; % 上层种群大小
	    llPopSize=200;  % 下层种群大小
	else
	    ulPopSize=50;
	    llPopSize=50;
	end
	% dimensionality

    % 设置维度
	switch fno
		case {1,2,3,5,6,7,8}
			ulDim = 2; llDim = 2;   % 上层和下层的维度都是2
		case 4
			ulDim = 2; llDim = 3;   % 上层维度是2，下层维度是3
		case 9
			ulDim = 5; llDim = 5;   % 上层和下层维度都是5
		case 10
			ulDim = 10; llDim = 10;  % 上层和下层维度都是10
    end

	%% boundaries and number of constraints
	switch fno

        % 根据目标函数编号设置边界和约束数量， tp1-》case1

		case 1
			ulDimMin = [-30 -30]; % 上层的维度最小值和最大值
			ulDimMax = [30 15];
			llDimMin = 0*ones(1,llDim); % 下层的维度最小值和最大值
			llDimMax = 10*ones(1,llDim); 
			ul_ieqcon_num = 2;% 上下层不等式约束数量
			ll_ieqcon_num = 0;
		case 2
			ulDimMin = zeros(1,ulDim);
			ulDimMax = 50*ones(1,ulDim);
			llDimMin = -10*ones(1,llDim);
			llDimMax = 20*ones(1,llDim);
			ul_ieqcon_num = 3;
			ll_ieqcon_num = 2;
		case 3
			ulDimMin = zeros(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = zeros(1,llDim);
			llDimMax = 10*ones(1,llDim);
			ul_ieqcon_num = 3;
			ll_ieqcon_num = 2;
		case 4
			ulDimMin = zeros(1,ulDim);
			ulDimMax = 1*ones(1,ulDim);
			llDimMin = zeros(1,llDim);
			llDimMax = 1*ones(1,llDim);
			ul_ieqcon_num = 3;
			ll_ieqcon_num = 3;
		case 5
			ulDimMin = zeros(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = zeros(1,llDim);
			llDimMax = 10*ones(1,llDim);
			ul_ieqcon_num = 2;
			ll_ieqcon_num = 2;
		case 6
			ulDimMin = zeros(1,ulDim);
			ulDimMax = 2*ones(1,ulDim);
			llDimMin = zeros(1,llDim);
			llDimMax = 2*ones(1,llDim);
			ul_ieqcon_num = 4;
			ll_ieqcon_num = 4;
		case 7
			ulDimMin = zeros(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = zeros(1,llDim);
			llDimMax = [1 10];
			ul_ieqcon_num = 4;
			ll_ieqcon_num = 2;
		case 8
			ulDimMin = zeros(1,ulDim);
			ulDimMax = 50*ones(1,ulDim);
			llDimMin = -10*ones(1,llDim);
			llDimMax = 20*ones(1,llDim);
			ul_ieqcon_num = 3;
			ll_ieqcon_num = 2;
		case 9
			ulDimMin = -1*ones(1,ulDim);
			ulDimMax = 1*ones(1,ulDim);
			llDimMin = -pi*ones(1,llDim);
			llDimMax = pi*ones(1,llDim);
			ul_ieqcon_num = 0;
			ll_ieqcon_num = 0;
		case 10
			ulDimMin = -1*ones(1,ulDim);
			ulDimMax = 1*ones(1,ulDim);
			llDimMin = -pi*ones(1,llDim);
			llDimMax = pi*ones(1,llDim);
			ul_ieqcon_num = 0;
			ll_ieqcon_num = 0;
	end
	ul_eqcon_num = 0;
	ll_eqcon_num = 0;

	%% optimal values




    % 设置目标函数的最优值


	% from BLEAQ2
	% ulBestKnownFunctionValue = [225, 0, -18.6787, -29.2, -3.6, -1.2091, -1.96, 0, 0, 0];
    % llBestKnownFunctionValue = [100, 100, -1.0156, 3.2, -2.0, 7.6145, 1.96, 100, 1.0, 1.0];
    % my results
%     does TP3 has better results?



    % 设置目标函数的最优值
	ulBestKnownFunctionValue = [225, 0, -18.6787, -29.2, -3.6, -1.20987, -1.96146, 0, 0, 0];
    llBestKnownFunctionValue = [100, 100, -1.0156, 3.2, -2.0, 7.61728, 1.96146, 100, 1.0, 1.0];
    

    % 获取目标函数的最优值
    ulOpt = -ulBestKnownFunctionValue(fno);  % 上层最优值
    llOpt = -llBestKnownFunctionValue(fno);    % 下层最优值
end




 
if strcmp(benchmark,'SMD')  % 如果基准函数是 'SMD'
	fn = sprintf('smd%d',fno); % 如果基准函数是 'SMD'
	assert(ismember(fno,1:12));  % 断言目标函数编号在1到12之间
	ulStoppingCriteria = 1e-4;
	llStoppingCriteria = 1e-5;
	eps = 0.00001;



      % 处理维度的设置
    if nargin == 3
    	if dim == 5
    		ulDim = 2;
    		llDim = 3;
		elseif dim == 10 || dim == 20
            ulDim = dim / 2;
            llDim = dim / 2;
        else
            error('no supported dimensionality');
        end
        % if dim >= 5 && mod(dim,5) == 0
        %     ulDim = dim * 2 / 5;
        %     llDim = dim - ulDim;
        % else
        %     error('no supported dimensionality');
        % end
    else
        ulDim=2; llDim=3;
    end


    %% check the constraints
	isLowerLevelConstraintsIncludedInUpperLevel = false;
    %% population size

      % 设定种群大小
    ulPopSize=20;
    llPopSize=30;
    %% decision variable decomposition


    % 设定决策变量的分解
    switch fno
    	case 6
			r = floor(ulDim/2);
			p = ulDim - r;
			q = floor((llDim - r)/2 - eps);
			s = ceil((llDim - r)/2 + eps);
		case {1,2,3,4,5,7,8,9,10,11,12}
			r = floor(ulDim/2); 
			p = ulDim - r; 
			q = llDim - r;
    end


    % 获取最优解  --->步入文件getOptimalSolutionSMD
    %% optimum function values
	[~,~,ulOpt,llOpt]=getOptimalSolutionSMD(ulDim,llDim,sprintf('smd%d',fno));




% 根据目标函数编号设置边界
    %% boundaries
	switch fno
		case 1 %ok
			ulDimMin = -5*ones(1,ulDim);                    
			ulDimMax = 10*ones(1,ulDim);                    
			llDimMin = [-5*ones(1,q) -pi/2*ones(1,r)+eps];  
			llDimMax = [10*ones(1,q)  pi/2*ones(1,r)-eps];  
		case 2 %ok
			ulDimMin = -5*ones(1,ulDim); 
			ulDimMax = [10*ones(1,p) 1*ones(1,r)];
			llDimMin = [-5*ones(1,q) eps*ones(1,r)];
			llDimMax = [10*ones(1,q) exp(1)*ones(1,r)];
		case 3 %ok
			ulDimMin = -5*ones(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = [-5*ones(1,q) -pi/2*ones(1,r)+eps];
			llDimMax = [10*ones(1,q)  pi/2*ones(1,r)-eps];
		case 4 %ok
			ulDimMin = [-5*ones(1,p) -1*ones(1,r)];
			ulDimMax = [10*ones(1,p) 1*ones(1,r)];
			llDimMin = [-5*ones(1,q) zeros(1,r)];
			llDimMax = [10*ones(1,q) exp(1)*ones(1,r)];
		case 5 %ok
			ulDimMin = -5*ones(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = [-5*ones(1,q) -5*ones(1,r)];
			llDimMax = [10*ones(1,q) 10*ones(1,r)];
		case 6 %ok
			ulDimMin = -5*ones(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = -5*ones(1,llDim);
			llDimMax = 10*ones(1,llDim);
		case 7 %ok
			ulDimMin = -5*ones(1,ulDim);
			ulDimMax = [10*ones(1,p) 1*ones(1,r)];
			llDimMin = [-5*ones(1,q) eps*ones(1,r)];
			llDimMax = [10*ones(1,q) exp(1)*ones(1,r)];
		case 8 %ok
			ulDimMin = -5*ones(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = [-5*ones(1,q) -5*ones(1,r)];
			llDimMax = [10*ones(1,q) 10*ones(1,r)];
		case 9 %ok
			ulDimMin = -5*ones(1,ulDim);
			ulDimMax = [10*ones(1,p) 1*ones(1,r)];
			llDimMin = [-5*ones(1,q) -1+eps*ones(1,r)];
			llDimMax = [10*ones(1,q) -1+exp(1)*ones(1,r)];
		case 10	%ok		
			ulDimMin = -5*ones(1,ulDim);
			ulDimMax = 10*ones(1,ulDim);
			llDimMin = [-5*ones(1,q) -pi/2*ones(1,r)+eps];
			llDimMax = [10*ones(1,q)  pi/2*ones(1,r)-eps];
		case 11 %ok
			ulDimMin = [-5*ones(1,p) -1*ones(1,r)];
			ulDimMax = [10*ones(1,p) 1*ones(1,r)];
			llDimMin = [-5*ones(1,q) 1/exp(1)*ones(1,r)];
			llDimMax = [10*ones(1,q) exp(1)*ones(1,r)];
		case 12 %???? inconsistent with the paper
			ulDimMin = [-5*ones(1,p) -1*ones(1,r)];
			ulDimMax = [10*ones(1,p) 1*ones(1,r)];
			llDimMin = [-5*ones(1,q) -pi/4*ones(1,r)+eps];
			llDimMax = [10*ones(1,q)  pi/4*ones(1,r)-eps];
    end



	%% number of constraints
	switch fno  %约束的数量 
		case {1,2,3,4,5,6,7,8}
			ul_ieqcon_num = 0;
			ll_ieqcon_num = 0;
		case 9
			ul_ieqcon_num = 1;
			ll_ieqcon_num = 1;
		case 10
			ul_ieqcon_num = p+r;
			ll_ieqcon_num = q;
		case 11
			ul_ieqcon_num = r;
			ll_ieqcon_num = 1;
		case 12
			ul_ieqcon_num = p+2*r;
			ll_ieqcon_num = q+1;
	end
	ul_eqcon_num = 0;
	ll_eqcon_num = 0;
end



if strcmp(benchmark,'DecisionMaking')
	fn = benchmark;
	ulStoppingCriteria = 1e-4;
	llStoppingCriteria = 1e-5;

	isLowerLevelConstraintsIncludedInUpperLevel = false;

	%% population size
    ulPopSize=50; llPopSize=50;
	% dimensionality
	ulDim = 3; llDim = 3;
	%% boundaries and number of constraints
	ulDimMin = [0 0 0];
	ulDimMax = [250 250 1];
	llDimMin = [0 0 0];
	llDimMax = [70 70 70];
	ul_ieqcon_num = 4;
	ll_ieqcon_num = 5;
	ul_eqcon_num = 0;
	ll_eqcon_num = 0;

	%% optimal values
	ulBestKnownFunctionValue = 0;
    llBestKnownFunctionValue = 0;

    ulOpt = -ulBestKnownFunctionValue(1);
    llOpt = -llBestKnownFunctionValue(1);
end
if strcmp(benchmark,'GoldMining')
	fn = benchmark;
	ulStoppingCriteria = 1e-4;
	llStoppingCriteria = 1e-5;

	isLowerLevelConstraintsIncludedInUpperLevel = false;

	%% population size
    ulPopSize=50; llPopSize=50;
	% dimensionality
	ulDim = 2; llDim = 1;
	%% boundaries and number of constraints
	ulDimMin = [0 0];
	ulDimMax = [100 1];
	llDimMin = 0;
	llDimMax = 100;
	ul_ieqcon_num = 0;
	ll_ieqcon_num = 1;
	ul_eqcon_num = 0;
	ll_eqcon_num = 0;

	%% optimal values
	ulBestKnownFunctionValue = 0;
    llBestKnownFunctionValue = 0;

    ulOpt = -ulBestKnownFunctionValue(1);
    llOpt = -llBestKnownFunctionValue(1);
end


% 设置并返回最终的结构体 BI，包含所有信息
BI.u_ieqcon_num = ul_ieqcon_num; % 设置上层不等式约束的数量
BI.l_ieqcon_num = ll_ieqcon_num;
BI.u_eqcon_num = ul_eqcon_num;
BI.l_eqcon_num = ll_eqcon_num; % 设置下层等式约束的数量


BI.u_dim = ulDim; % 设置上层的决策变量维度
BI.l_dim = llDim;
BI.dim = ulDim + llDim; % 设置总的决策变量维度（上层 + 下层）

BI.u_lb = ulDimMin; % 设置上层决策变量的下界（最小值）
BI.u_ub = ulDimMax;
BI.l_lb = llDimMin;
BI.l_ub = llDimMax; % 设置下层决策变量的上界（最大值）



% 设置上下层决策变量的范围，合并了上层和下层的最小值与最大值
BI.xrange = [ulDimMin llDimMin;ulDimMax llDimMax];


% 设置容忍度，用于停止准则判定
BI.u_ftol = 1e-6; % 上层目标函数的容忍度
BI.l_ftol = 1e-6;

% from the original BLEAQ2 source code  % 设置停止准则
BI.ulStoppingCriteria = ulStoppingCriteria;
BI.llStoppingCriteria = llStoppingCriteria;


% 设置最优函数值，通常是在基准问题中预设的最优值
BI.u_fopt = ulOpt;
BI.l_fopt = llOpt;

BI.fn = fn;

% from BLEAQ2% 设置上下层的种群大小
BI.u_N = ulPopSize;
BI.l_N = llPopSize;
BI.u_maxGen = 2000;% 设置最大迭代次数
BI.l_maxGen = 2000;
% 设置上下层是否包含下层约束（用于优化过程的约束处理）
BI.isLowerLevelConstraintsIncludedInUpperLevel = isLowerLevelConstraintsIncludedInUpperLevel;