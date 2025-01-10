
% @Author: hxy
% @Date:   2017-11-15 11:48:12
% @Last Modified by:   hxy
% @Last Modified time: 2018-03-06 20:00:16

% i --- injection correction
% e --- elite preservation
% c --- constrained handling
% s --- advanced stopping criteria

% 2017/11/15 modify some termination criteria

function ins = MOTEA(BI)

% MOTEA 多目标优化算法实现，基于CMA-ES和上下层函数评估
% 输入:
%   BI - 包含优化问题相关参数的结构体
% 输出:
%   ins - 包含优化结果和记录的结构体

 % 声明全局变量，用于记录上层和下层函数的评估次数
global ulFunctionEvaluations;
global llFunctionEvaluations;
ulFunctionEvaluations = 0; % 初始化上层函数评估计数器
llFunctionEvaluations = 0;  % 初始化下层函数评估计数器

elite = []; % 初始化精英个体
record = []; % 初始化记录结构，用于存储每次迭代的精英个体
CMA = initCMAES(BI); % 初始化CMA-ES参数，调用自定义初始化函数


   % 设置上下界，将上层和下层的上下界组合在一起
BI.ub = [BI.u_ub,BI.l_ub];
BI.lb = [BI.u_lb,BI.l_lb];



    % 计算最大迭代次数和判断平稳状态所需的迭代次数
maxIter = ceil(BI.UmaxFEs/CMA.lambda);
ImprIter = ceil(BI.UmaxImprFEs/CMA.lambda);


% %--------------------------------------------------------------
% % 初始化计数器
global cmaes_wins temp_wins;
cmaes_wins = 0;   % CMA-ES 生成的个体胜利次数
temp_wins = 0;    % Temp 个体胜利次数
% 

global cmaes_wins_record temp_wins_record;

% % 创建一个空的数组，用来存储每轮的胜利次数
cmaes_wins_record = [];
temp_wins_record = [];
% %---------------------------------------------------------------

  % 主循环，最多迭代maxIter次  350代
for iter = 1 : maxIter  
     % 采样阶段，生成新一代种群（种群数量为Lambda）
    % sampling
    for i = 1 : CMA.lambda


        % 根据CMA-ES的当前均值和协方差矩阵生成新的个体U
        U = CMA.xmean + CMA.sigma * randn(1,BI.dim) .* CMA.D * CMA.B';   %第一步根据CMA-ES生产新的个体并修正


         % 检查并修正U是否超过定义的范围上界
        flag = U > BI.xrange(2,:);
        U(flag) = (CMA.xmean(flag) + BI.xrange(2,flag))/2;% 修复方法为，通过cmaxmean对应flag位置的中心均值和上界对应的值求和除2得到


         % 检查并修正U是否低于定义的范围下界
        flag = U < BI.xrange(1,:);
        U(flag) = (CMA.xmean(flag) + BI.xrange(1,flag))/2;


         % 创建一个新的种群个体，并赋值相关属性

% 
%             case 'popclass' 
%         str = struct(...
%             ...%%%%%%%%%%%%%%%%%%%%%%%%%%%公有参数
%             'X',[],...               上下层合并的变量，称为种群个体
%             'UX',[],...                上层决策变量
%             'UF',[],...                %inter是子种群
%             'UC',[],...              % num_ind子种群的规模  
%             'LX',[],...                下层决策变量
%             'LF',[],...                %inter是子种群
%             'LC',[],...              % num_ind子种群的规模  
%             'UFEs',[],...                %center是子区域的中心
%             'LFEs',[],...                %inter是子种群
%             'RF',[],...              % num_ind子种群的规模  
%             'W',[],...
%             'fit',[]);               % 每个子种群的权重
        


        %pop结构体，属性如上
        POP(i)  = get_structure('popclass');
        POP(i).X  = U;
        POP(i).UX = U(1:BI.u_dim);
        POP(i).LX = U(BI.u_dim+1:end);    
    end  



            % 并行执行下层搜索，更新种群个体的下层变量
    % parallel lower level search
    POP = parallelLLS(POP,CMA,BI,iter/maxIter);   %根据CMA-ES生成的新个体进行下层搜索







    % 评估每个种群个体的上层函数
    for i = 1 : CMA.lambda
        [POP(i).UF,POP(i).UC] = UL_evaluate(POP(i).UX,POP(i).LX,BI.fn);
        POP(i).UFEs = ulFunctionEvaluations;
        POP(i).LFEs = llFunctionEvaluations;
    end



      % 分配适应度并进行精炼
    % fitness assignment && refinement
    POP = assignUpperFitness(POP,BI);



          % 查找当前最佳和最优解
    % find current best and optimal solution
    RF_idx = find([POP.RF]);  % 找到标记为RF的个体索引



    if isempty(RF_idx) 

        RF_idx = 1:length(POP);  % 如果没有标记，则考虑所有个体
    end



    [~,bestIdx] = max([POP(RF_idx).fit]); % 找到适应度最高的个体索引
    bestIdx = RF_idx(bestIdx);  % 获取实际索引
    bestIndv = POP(bestIdx);  % 获取最佳个体




        % 更新精英个体
    % replace the elite   %rand没弄，只用看前面那个即可
    if upperLevelComparator(bestIndv,elite,BI) || rand > 0.5
    	bestIndv = refine(bestIndv,CMA,BI);  % 精炼最佳个体
    	POP(bestIdx) = bestIndv;     % 更新种群中的最佳个体
    	if upperLevelComparator(bestIndv,elite,BI)
    		elite = bestIndv;    % 如果最佳个体优于精英，则更新精英
    	end
    else
    	elite = refine(elite,CMA,BI);  % 否则精炼当前精英
    end



    % 重新分配适应度
    % re-assign fitness
    POP = assignUpperFitness(POP,BI);




      % 记录当前精英个体的状态
    % recording
    elite.UFEs = ulFunctionEvaluations;    
    elite.LFEs = llFunctionEvaluations;
    record = [record;elite];  % 将当前精英添加到记录中





      %% 终止条件检查
    %% termination check
    reachTheTarget_b = abs(bestIndv.UF-BI.u_fopt) < BI.u_ftol;  % 检查最佳个体是否达到目标
    reachTheTarget_e = abs(elite.UF-BI.u_fopt) < BI.u_ftol;    % 检查精英个体是否达到目标
    reachTheMaxFEs = ulFunctionEvaluations >= BI.UmaxFEs;   % 检查是否达到最大函数评估次数
    reachTheFlat = iter>ImprIter && abs(record(iter).UF-record(iter-ImprIter+1).UF)/(abs(record(1).UF)+abs(record(iter).UF)) < BI.u_ftol;   % 检查是否达到平稳状态



    reachTheFlat = reachTheFlat && abs(record(iter).UF-record(iter-ImprIter+1).UF) < BI.u_ftol; % 进一步检查变化量




    % 如果满足任一终止条件，则退出循环
    if reachTheTarget_b || reachTheTarget_e || reachTheMaxFEs || reachTheFlat
        if reachTheTarget_b
            elite = bestIndv;  % 如果最佳个体达到目标，则更新精英为最佳个体
        end
        break;
    end

    % 更新CMA-ES参数以进行下一代迭代
    CMA = updateCMAES(CMA,POP,BI);
end




% 保存胜利次数记录到文件
save('victory_counts.mat', 'cmaes_wins_record', 'temp_wins_record');

  % 将优化结果存储在输出结构体中
ins.UF = elite.UF;  % 精英个体的上层函数值
ins.LF = elite.LF;  % 精英个体的下层函数值
ins.UX = elite.UX;  % 精英个体的上层决策变量
ins.LX = elite.LX;  % 精英个体的下层决策变量
ins.UFEs = ulFunctionEvaluations; % 总上层函数评估次数
ins.LFEs = llFunctionEvaluations; % 总下层函数评估次数
ins.record = record;  % 记录每次迭代的精英个体






function POP = parallelLLS(POP,CMA,BI,perct)  
% parallelLLS 并行下层搜索函数
% 输入:
%   POP - 当前种群
%   CMA - CMA-ES 参数结构体
%   BI - 优化问题相关参数结构体
%   perct - 当前迭代的进度比例
% 输出:
%   POP - 更新后的种群
    global cmaes_wins temp_wins;
    global cmaes_wins_record temp_wins_record;

    % 重塑上层决策变量矩阵
     UX = reshape([POP.UX],BI.u_dim,[]); % 将决策变量数组变成列向量后形成矩阵


      % 计算调整比例 p
     p = 1/(1+exp(-6*perct));


     % 动态调整上层变量的上下界
     BI.ux_ub = BI.u_ub*(1-p)+(p)*max(UX,[],2)';
     BI.ux_lb = BI.u_lb*(1-p)+(p)*min(UX,[],2)';



     Num = 2; % 分组数量， 分两组
     Lpopsize = Num*CMA.lambda;  % 下层种群大小
     Maxgen = ceil(BI.LmaxFEs/Num);  % 最大生成次数



    ImprIter = 5;   % 平稳判断迭代次数
    %initialize



    % 初始化子种群结构
    sub    = get_structure('subclass');
    
%     case 'subclass' 
%         str = struct(...
%             ...%%%%%%%%%%%%%%%%%%%%%%%%%%%公有参数
%             'center',[],...               %个体权重向量反映了每个个体的相对重要性或者影响力
%             'pdist',[],...                %center是子区域的中心
%             'inter',[],...                %inter是子种群
%             'num_ind',[],...              % num_ind子种群的规模  
%             'RF',0,...                   %   
%             'temp',[],...                   %   
%             'idealmin',[],...             % 
%             'weight',[]);                 % 每个子种群的权重




    subpop = repmat(sub, [1,Lpopsize/Num]);%复制10个结构体
    NNbs   = 3;   % 邻居数量
    


     % 初始化子种群
    for i = 1 : Lpopsize/Num
        subpop(i).center = generateweight(POP(i).UX,BI);    % 生成权重中心
        subpop(i).sigma = CMA.sigma;    % 设置初始步长
        [POP(i).LF,POP(i).LC] = LL_evaluate(POP(i).UX,POP(i).LX,BI.fn);        % 评估下层函数

        subpop(i).elite = POP(i); % 初始化精英个体



        %---------------------------------------------------------------
        Temp = POP(i);
        Temp.LX = CMA.xmean(BI.u_dim+1:end); % 设置临时下层变量
        Temp.X = [Temp.UX,Temp.LX];
        %---------------------------------------------------------------

        [Temp.LF,Temp.LC] = LL_evaluate(Temp.UX,Temp.LX,BI.fn);         % 评估临时下层函数 



        if lowerLevelComparator(Temp,subpop(i).elite) % 先比较适应度 

%   isNoWorseThan - 如果 P 的适应度不差于 Q，则为 true，否则为 false,Q为空返回true
%比较 P 和 Q 的适应度，p>q返回true
            
             % 记录 Temp 的胜利
            temp_wins = temp_wins + 1;


            subpop(i).elite = Temp;    % 更新精英个体
            POP(i) = Temp; 
            C1 = CMA.C(BI.u_dim+1:end,BI.u_dim+1:end);      % 更新协方差矩阵      
            subpop(i).sigma = CMA.sigma/2; % 调整步长
            subpop(i).xmean = CMA.xmean(BI.u_dim+1:end);
        else



            Temp_C = (POP(i).LX - Temp.LX)'*(POP(i).LX - Temp.LX);   % 计算协方差修正
            C1 = CMA.C(BI.u_dim+1:end,BI.u_dim+1:end)*(1-perct)+perct*Temp_C;
            C1 = triu(C1) + triu(C1,1)';  % 对称化协方差矩阵
            subpop(i).sigma = CMA.sigma; % 保持步长不变
            subpop(i).xmean = CMA.xmean(BI.u_dim+1:end)*(1-perct)+perct*POP(i).LX;   % 更新均值


            % 记录 CMA-ES 生成个体的胜利
            cmaes_wins = cmaes_wins + 1;
        end       



        [B1,D1] = eig(C1);   % 计算特征值分解
        D1 = sqrt(diag(D1))'; 
        subpop(i).B = B1; % 设置特征向量矩阵
        subpop(i).D = D1;  % 设置特征值根
        subpop(i).RF = 0;  % 初始化停止标志
        subpop(i).inter = [];  % 初始化中间个体
    end

    % 记录每轮的胜利次数
    cmaes_wins_record = [cmaes_wins_record, cmaes_wins];
    temp_wins_record = [temp_wins_record, temp_wins];


    Q = [];   % 初始化临时种群



    for i = 1 : Lpopsize/Num
         [subpop(i).snbs,subpop(i).neighbours ]= findneighbours(i,subpop,NNbs); % 查找邻居
        for j = 1 : 2*Num
            subQ(j) = get_structure('popclass'); % 创建新个体结构
            subQ(j).UX = subpop(i).elite.UX; % 继承上层决策变量
            Temp_LX = subpop(i).xmean;        

            subQ(j).LX = Temp_LX + subpop(i).sigma*randn(1,BI.l_dim).* subpop(i).D * subpop(i).B'; % 生成新下层变量
            flag = subQ(j).LX > BI.l_ub;  % 检查上界


            subQ(j).LX(flag) = (Temp_LX(flag) + BI.l_ub(flag))/2;  % 修正上界超出部分
            flag = subQ(j).LX < BI.l_lb; % 检查下界
            subQ(j).LX(flag) = (Temp_LX(flag) + BI.l_lb(flag))/2;  % 修正下界超出部分


            subQ(j).X = [subQ(j).UX,subQ(j).LX];  % 合并决策变量


            [subQ(j).LF,subQ(j).LC] = LL_evaluate(subQ(j).UX,subQ(j).LX,BI.fn);      % 评估下层函数


        end 
        subpop(i).inter = subQ;  % 保存中间个体
            Q = [Q,subQ];   % 添加到临时种群

            clear subQ; % 清除临时变量
    end
       Q = [Q,POP]; % 合并当前种群


       W = generateweight([Q.UX],BI);    % 生成权重


       team = group(subpop, W);  % 分组


       theta = 5;  % 选择参数


       subpop = select(subpop,Q,team,BI,Num,theta,perct);     % 选择操作  


       LPOP = [subpop.inter];   % 获取中间个体

       Record = zeros(Lpopsize/Num,Maxgen);  % 初始化记录矩阵



       for i = 1 : Lpopsize/Num
           Temp = subpop(i).temp;         
           if sum(abs(Temp.LX-subpop(i).elite.LX))>0
               Temp.UX = POP(i).UX;
               Temp.X = [Temp.UX,Temp.LX];
               [Temp.LF,Temp.LC] = LL_evaluate(POP(i).UX,Temp.LX,BI.fn); % 重新评估
               if lowerLevelComparator(Temp,subpop(i).elite)
                   subpop(i).elite = Temp; % 更新精英
               end
           end
           Record(i,1) = subpop(i).elite.LF; % 记录初始下层函数值
       end


          % 迭代生成新的下层种群
   
    for gen = 2:Maxgen      
        Q = []; % 重置临时种群


        for i = 1 : Lpopsize/Num
            Inter_temp = [subpop(i).inter];
            Num_temp = length(Inter_temp);
            if Num_temp > 0
                for j = 1 : Num_temp    
                    subQ(j) = get_structure('popclass'); % 创建新个体结构



                    X= CM(j,i,subpop,Num,BI);  %generate new offspring by SBX     % 通过交叉产生新个体


                    if all([subpop([i,subpop(i).neighbours]).RF] ==1)&&j==1 
                       subQ(j).UX = subpop(i).elite.UX;    % 保持上层变量不变
                    else                  
                       subQ(j).UX = X(1:BI.u_dim);       % 更新上层变量
                    end   
                    subQ(j).LX = X(BI.u_dim+1:end);  % 更新下层变量
                    subQ(j).X = [subQ(j).UX,subQ(j).LX]; % 合并决策变量   
                    [subQ(j).LF,subQ(j).LC] = LL_evaluate(subQ(j).UX,subQ(j).LX,BI.fn);  % 评估下层函数
                end


                Q = [Q,subQ];  % 添加到临时种群
                clear subQ;   % 清除临时变量
            end           
        end                 
%        fprintf('gen=%d\n', gen)     
      Q = [LPOP,Q];  % 合并中间个体
      W =  generateweight([Q.UX],BI);      % 生成权重
      team = group(subpop, W);  % 分组
%       theta = 5 + min(100,exp(10*perct))+ exp(10*gen/Maxgen);
      subpop = select(subpop,Q,team,BI,Num,theta,perct);       % 选择操作
      LPOP = [subpop.inter]; % 获取中间个体



      for i = 1 : Lpopsize/Num    
         Temp = subpop(i).temp;
          if sum(abs(Temp.LX-subpop(i).elite.LX))>0              
              Temp.UX = POP(i).UX;
              Temp.X = [Temp.UX,Temp.LX];
              [Temp.LF,Temp.LC] = LL_evaluate(POP(i).UX,Temp.LX,BI.fn);  % 重新评估
              if lowerLevelComparator(Temp,subpop(i).elite)
                  subpop(i).elite = Temp;  % 更新精英
              end           
          end     
         Record(i,gen) = subpop(i).elite.LF;  % 记录下层函数值


            % 检查是否达到平稳状态
         if ((gen>ImprIter&&abs(Record(i,gen -ImprIter+1)-Record(i,gen))/(max(1,(abs(Record(i,1))+abs(Record(i,gen))))) < 1e-4)...
                 ||(gen > ImprIter && abs(Record(i,gen -ImprIter+1)-Record(i,gen)) < 10*BI.l_ftol))
             subpop(i).RF = 1;  % 设置停止标志
         else
             subpop(i).RF = 0; % 继续优化
         end
      end

      % 如果所有子种群达到停止条件，退出迭代
%          X = reshape([LPOP.X],BI.dim,[]);
        if all([subpop.RF]>0)
             break;
        end
   end     


     % 更新主种群的下层变量和评估结果
     for i = 1 : Lpopsize/Num
         POP(i).LX = subpop(i).elite.LX;        % 更新下层变量
         POP(i).LF = subpop(i).elite.LF;        % 更新下层函数值
         POP(i).LC = subpop(i).elite.LC;       % 更新下层附加信息
         POP(i).RF = subpop(i).RF>0;             %更新停止标志
     end

function [SNbs,Nbs] = findneighbours(i,subpop,NNbs)
% findneighbours 查找指定个体的邻居
% 输入:
%   i     - 当前个体的索引
%   subpop - 子种群结构体数组
%   NNbs  - 邻居数量
% 输出:
%   SNbs - 与当前个体相同中心的个体索引
%   Nbs  - 最近邻的个体索引




W = [subpop.center];  % 获取所有子种群个体的中心权重
dst = W'*subpop(i).center;  % 计算当前个体中心与所有其他个体中心的点积（相似度度量）
[value,indx]  = sort(-dst);  % 对相似度进行降序排序，得到排序值和对应的索引
Num = length(value);   % 总个体数量
T = length(find(value==value(1))); % 找到与当前个体相同中心的个体数量
SNbs = indx(1:T)';    % 获取这些相同中心个体的索引



if T+NNbs<Num
    Nbs = indx(T+1:T+NNbs)';  % 获取最近的 NNbs 个邻居的索引
else
    Nbs = indx(T+1:Num)';   % 如果剩余个体不足，则获取所有剩余个体的索引
end 
  
  
function cfit_ = combineConstraintWithFitness(fit_,c) 
% combineConstraintWithFitness 将约束条件与适应度值结合
% 输入:
%   fit_ - 个体的适应度值数组
%   c    - 个体的约束违反值数组
% 输出:
%   cfit_ - 结合约束后的适应度值数组




if all(c==0)
	cfit_ = fit_;   % 如果所有个体都满足约束，适应度值不变
else
    isFeasible = c == 0; % 判断哪些个体满足约束
    if any(isFeasible)
        cfit_(isFeasible)  = fit_(isFeasible);  % 满足约束的个体适应度保持不变
        cfit_(~isFeasible) = min(fit_(isFeasible)) - c(~isFeasible);  % 违反约束的个体适应度降低
    else
        cfit_ = - c; % 如果所有个体都违反约束，适应度值为约束违反程度的负值
    end 
end



function POP = assignLowerFitness(POP)
% assignLowerFitness 为种群分配下层适应度
% 输入:
%   POP - 种群结构体数组，包含下层适应度和约束
% 输出:
%   POP - 更新后的种群，包含分配后的适应度值




fit_ = combineConstraintWithFitness([POP.LF],[POP.LC]); % 结合约束条件计算适应度
for i = 1 : length(POP)
    POP(i).fit = fit_(i); % 将计算后的适应度值赋值给每个个体
end

function POP = assignUpperFitness(POP,BI)
% assignUpperFitness 为种群分配上层适应度
% 输入:
%   POP - 种群结构体数组，包含上层适应度和约束
%   BI  - 优化问题相关参数结构体
% 输出:
%   POP - 更新后的种群，包含分配后的适应度值




CV = [POP.UC]; % 获取所有个体的上层约束违反值
if ~BI.isLowerLevelConstraintsIncludedInUpperLevel
    CV = CV + [POP.LC];  % 如果下层约束包含在上层中，合并下层约束
end
% CV = CV / max(CV);  % （注释掉）可选的约束归一化
fit_ = combineConstraintWithFitness([POP.UF],CV);   % 结合约束条件计算适应度


for i = 1 : length(POP)
	POP(i).fit = fit_(i);   % 将计算后的适应度值赋值给每个个体
	% ????? how to handle RF
	% according to latest experimental results, incorporating RF here seems to be meaningless
      % （注释内容：根据最新的实验结果，将 RF（是否停止标志）纳入适应度似乎没有意义）

end

function Q = refine(P,CMA,BI)
% refine 精炼个体，通过下层搜索优化
% 输入:
%   P    - 当前个体
%   CMA  - CMA-ES 参数结构体
%   BI   - 优化问题相关参数结构体
% 输出:
%   Q    - 精炼后的个体


Q = P;  % 初始化精炼后的个体为当前个体
[Q.LX,Q.LF,Q.LC,Q.RF] = lowerLevelSearch(Q.UX,CMA,BI);  % 执行下层搜索优化
if lowerLevelComparator(Q,P)
    Q.RF = max(Q.RF,P.RF);   % 更新停止标志
    [Q.UF,Q.UC] = UL_evaluate(Q.UX,Q.LX,BI.fn); % 重新评估上层函数
else
	Q = P; % 如果精炼后的个体不优于原个体，则保持原个体不变
end

function isNoWorseThan = upperLevelComparator(P,Q,BI)
% upperLevelComparator 比较两个个体的上层适应度
% 输入:
%   P - 个体 P
%   Q - 个体 Q
%   BI - 优化问题相关参数结构体
% 输出:
%   isNoWorseThan - 如果 P 的适应度不差于 Q，则为 true，否则为 false



if isempty(Q)  
    isNoWorseThan = true;  % 如果 Q 为空，则 P 被视为不差
else
    tmp = assignUpperFitness([P Q],BI);  % 为 P 和 Q 分配上层适应度
    isNoWorseThan = tmp(1).fit >= tmp(2).fit;  % 比较 P 和 Q 的适应度
end

function isNoWorseThan = lowerLevelComparator(P,Q)
% lowerLevelComparator 比较两个个体的下层适应度
% 输入:
%   P - 个体 P
%   Q - 个体 Q
% 输出:
%   isNoWorseThan - 如果 P 的适应度不差于 Q，则为 true，否则为 false





if isempty(Q)
    isNoWorseThan = true;  % 如果 Q 为空，则 P 被视为不差
else
    tmp = assignLowerFitness([P Q]);   % 为 P 和 Q 分配下层适应度
    isNoWorseThan = tmp(1).fit >= tmp(2).fit;  % 比较 P 和 Q 的适应度
end 

function [bestLX,bestLF,bestLC,bestRF] = lowerLevelSearch(xu,CMA,BI)
% lowerLevelSearch 执行下层搜索优化
% 输入:
%   xu  - 上层决策变量
%   CMA - CMA-ES 参数结构体
%   BI  - 优化问题相关参数结构体
% 输出:
%   bestLX - 最佳下层决策变量
%   bestLF - 最佳下层适应度值
%   bestLC - 最佳下层约束违反值
%   bestRF - 是否达到停止条件





sigma0 = 1;  % 初始步长
LCMA.xmean = CMA.xmean(BI.u_dim+1:end); % 初始化下层均值
LCMA.sigma = sigma0;  % 设置下层步长


LCMA.C = CMA.C(BI.u_dim+1:end,BI.u_dim+1:end) * CMA.sigma^2;  % 初始化下层协方差矩阵


LCMA.pc = CMA.pc(BI.u_dim+1:end) * CMA.sigma; % 初始化下层进化路径 pc
LCMA.ps = zeros(1,BI.l_dim);  % 初始化下层进化路径 ps
lambda = 4+floor(3*log(BI.l_dim)); % 设置种群大小 lambda


mu = floor(lambda/2);  % 选择父代数量 mu
weights = log(mu+1/2)-log(1:mu);  % 计算权重



weights = weights/sum(weights);   % 归一化权重
mueff=sum(weights)^2/sum(weights.^2);  % 有效选择数量


cc = (4+mueff/BI.l_dim) / (BI.l_dim+4 + 2*mueff/BI.l_dim);  % 进化路径更新常数 cc



cs = (mueff+2) / (BI.l_dim+mueff+5); % 协方差更新常数 cs



c1 = 2 / ((BI.l_dim+1.3)^2+mueff);  % 学习率 c1



cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((BI.l_dim+2)^2+mueff));  % 学习率 cmu


damps = 1 + 2*max(0, sqrt((mueff-1)/(BI.l_dim+1))-1) + cs; % 步长调整参数 damps


chiN=BI.l_dim^0.5*(1-1/(4*BI.l_dim)+1/(21*BI.l_dim^2));  % 理想收敛路径长度


[LCMA.B,LCMA.D] = eig(LCMA.C);  % 计算协方差矩阵的特征分解



LCMA.D = sqrt(diag(LCMA.D))';   % 获取特征值的平方根


LCMA.invsqrtC = LCMA.B * diag(LCMA.D.^-1) * LCMA.B'; % 计算协方差矩阵的逆平方根


cy = sqrt(BI.l_dim)+2*BI.l_dim/(BI.l_dim+2);   % 缩放因子 cy
bestIndv = []; % 初始化最佳个体
bestRF = false;  % 初始化停止标志



maxIter = ceil(BI.LmaxFEs / lambda);    % 计算最大迭代次数

ImprIter = ceil(BI.LmaxImprFEs / lambda);  % 计算判断平稳状态的迭代次数


record = zeros(1,maxIter);  % 初始化记录数组
for iter = 1 : maxIter
	for i = 1 : lambda

         % 生成新个体的下层决策变量
	    Q(i).LX = LCMA.xmean + LCMA.sigma * randn(1,BI.l_dim) .* LCMA.D * LCMA.B';
	    flag = Q(i).LX > BI.l_ub; % 检查上层下层变量是否超过上界


	    Q(i).LX(flag) = (LCMA.xmean(flag) + BI.l_ub(flag))/2;
	    flag = Q(i).LX < BI.l_lb;  % 修正超出上界的部分


	    Q(i).LX(flag) = (LCMA.xmean(flag) + BI.l_lb(flag))/2;  % 检查下层变量是否低于下界
	    [Q(i).LF,Q(i).LC] = LL_evaluate(xu,Q(i).LX,BI.fn); % 修正低于下界的部分
    end

	Q = assignLowerFitness(Q);  % 评估下层函数% 为生成的个体分配下层适应度
 
	[~, Xindex] = sort(-[Q.fit]);     % 按适应度降序排序，获取排序索引
	xold = LCMA.xmean;    % 保存旧的均值 


	X = cell2mat(arrayfun(@(q)q.LX,Q','UniformOutput',false));  % 获取前 mu 个个体的下层决策变量

	Y = bsxfun(@minus,X(Xindex(1:mu),:),xold) / LCMA.sigma;  % 计算偏差

	Y = bsxfun(@times,Y,min(1,cy./sqrt(sum((Y*LCMA.invsqrtC').^2,2))));  % 缩放偏差

	delta_xmean = weights * Y; % 计算均值的变化

	LCMA.xmean = LCMA.xmean + delta_xmean * LCMA.sigma;  % 更新均值

	C_mu = Y' * diag(weights) * Y;  % 计算协方差矩阵的更新部分

	LCMA.ps = (1-cs)*LCMA.ps + sqrt(cs*(2-cs)*mueff) * delta_xmean * LCMA.invsqrtC; % 更新进化路径 ps

	LCMA.pc = (1-cc)*LCMA.pc + sqrt(cc*(2-cc)*mueff) * delta_xmean;  % 更新进化路径 pc

	LCMA.C = (1-c1-cmu) * LCMA.C + c1 * (LCMA.pc'*LCMA.pc) + cmu * C_mu;   % 更新协方差矩阵


	delta_sigma = (cs/damps)*(norm(LCMA.ps)/chiN - 1); % 计算步长的变化量


	LCMA.sigma = LCMA.sigma * exp(min(CMA.delta_sigma_max,delta_sigma));  % 更新步长，限制其变化范围

	LCMA.C = triu(LCMA.C) + triu(LCMA.C,1)';  % 对称化协方差矩阵

	[LCMA.B,LCMA.D] = eig(LCMA.C);  % 重新计算特征分解

	LCMA.D = sqrt(diag(LCMA.D))';  % 获取特征值的平方根

    LCMA = repairCMA(LCMA);  % 修复 CMA-ES 参数（自定义函数
	LCMA.invsqrtC = LCMA.B * diag(LCMA.D.^-1) * LCMA.B';   % 重新计算协方差矩阵的逆平方根

	% elite-preservation right????

     % 精英个体保留
	if lowerLevelComparator(Q(Xindex(1)),bestIndv)
		bestIndv = Q(Xindex(1));   % 更新最佳个体
	end
    record(iter) = bestIndv.LF;    % 记录当前最佳下层适应度


    % 检查是否满足终止条件
    if (iter>ImprIter && abs(record(iter)-record(iter-ImprIter+1))/(abs(record(1))+abs(record(iter))) < 1e-4) ...
		|| (iter > ImprIter && abs(record(iter) - record(iter-ImprIter+1)) < 10*BI.l_ftol) ...
        || LCMA.sigma / sigma0 < 1e-2 ...
        || LCMA.sigma / sigma0 > 1e2
    	bestRF = true;  % 设置停止标志
        break;  % 退出迭代
    end
end
bestLX = bestIndv.LX; % 获取最佳下层决策变量
bestLF = bestIndv.LF;  % 获取最佳下层适应度
bestLC = bestIndv.LC; % 获取最佳下层约束违反值


function CMA = initCMAES(BI)
% initCMAES 初始化 CMA-ES 参数
% 输入:
%   BI - 优化问题相关参数结构体
% 输出:
%   CMA - 初始化后的 CMA-ES 参数结构体


CMA.lambda = 4+floor(3*log(BI.dim));  % 计算种群大小 lambda


% CMA.sigma = min(BI.xrange(2,:)-BI.xrange(1,:))/2; % （注释掉）另一种初始化步长的方法


CMA.sigma = 0.3*median(BI.xrange(2,:)-BI.xrange(1,:));    % 初始化步长为变量范围中位数的 0.3 倍


CMA.mu = floor(CMA.lambda/2);   % 计算父代数量 mu


CMA.weights = log(CMA.mu+1/2)-log(1:CMA.mu);   % 计算权重

CMA.weights = CMA.weights/sum(CMA.weights);  % 归一化权重

CMA.mueff=sum(CMA.weights)^2/sum(CMA.weights.^2);  % 有效选择数量

CMA.cc = (4+CMA.mueff/BI.dim) / (BI.dim+4 + 2*CMA.mueff/BI.dim); % 进化路径更新常数 cc

CMA.cs = (CMA.mueff+2) / (BI.dim+CMA.mueff+5);  % 协方差更新常数 cs

CMA.c1 = 2 / ((BI.dim+1.3)^2+CMA.mueff);  % 学习率 c1

CMA.cmu = min(1-CMA.c1, 2 * (CMA.mueff-2+1/CMA.mueff) / ((BI.dim+2)^2+CMA.mueff)); % 学习率 cmu

CMA.damps = 1 + 2*max(0, sqrt((CMA.mueff-1)/(BI.dim+1))-1) + CMA.cs; % 步长调整参数 damps

CMA.chiN=BI.dim^0.5*(1-1/(4*BI.dim)+1/(21*BI.dim^2)); % 理想收敛路径长度

CMA.pc = zeros(1,BI.dim);   % 初始化进化路径 pc

CMA.ps = zeros(1,BI.dim);    % 初始化进化路径 ps
 

CMA.B = eye(BI.dim); % 初始化特征向量矩阵 B 为单位矩阵

CMA.D = ones(1,BI.dim);  % 初始化特征值向量 D 为全 1

CMA.C = CMA.B * diag(CMA.D.^2) * CMA.B'; % 初始化协方差矩阵 C

CMA.invsqrtC = CMA.B * diag(CMA.D.^-1) * CMA.B';   % 计算协方差矩阵的逆平方根

CMA.xmean = (BI.xrange(2,:)-BI.xrange(1,:)).*rand(1,BI.dim)+BI.xrange(1,:);  % 初始化均值为变量范围内的随机值

CMA.cy = sqrt(BI.dim)+2*BI.dim/(BI.dim+2); % 缩放因子 cy


CMA.delta_sigma_max = 1;  % 最大步长变化限制

function CMA = updateCMAES(CMA,POP,BI)
% updateCMAES 更新 CMA-ES 参数
% 输入:
%   CMA - 当前 CMA-ES 参数结构体
%   POP - 当前种群结构体数组，包含适应度值
%   BI  - 优化问题相关参数结构体
% 输出:
%   CMA - 更新后的 CMA-ES 参数结构体



[~, Xindex] = sort(-[POP.fit]);  % 按适应度降序排序，获取排序索引


xold = CMA.xmean;  % 保存旧的均值


X = cell2mat(arrayfun(@(p)[p.UX p.LX],POP(Xindex(1:CMA.mu))','UniformOutput',false));   % 获取前 mu 个个体的决策变量


Y = bsxfun(@minus,X,xold)/CMA.sigma; % 计算偏差

Y = bsxfun(@times,Y,min(1,CMA.cy./sqrt(sum((Y*CMA.invsqrtC').^2,2))));   % 缩放偏差

delta_xmean = CMA.weights * Y;   % 计算均值的变化

CMA.xmean = CMA.xmean + delta_xmean * CMA.sigma; % 更新均值

C_mu = Y' * diag(CMA.weights) * Y; % 计算协方差矩阵的更新部分

CMA.ps = (1-CMA.cs)*CMA.ps + sqrt(CMA.cs*(2-CMA.cs)*CMA.mueff) * delta_xmean * CMA.invsqrtC;  % 更新进化路径 ps

CMA.pc = (1-CMA.cc)*CMA.pc + sqrt(CMA.cc*(2-CMA.cc)*CMA.mueff) * delta_xmean; % 更新进化路径 pc

CMA.C = (1-CMA.c1-CMA.cmu) * CMA.C + CMA.c1 * (CMA.pc'*CMA.pc) + CMA.cmu * C_mu;  % 更新协方差矩阵

delta_sigma = (CMA.cs/CMA.damps)*(norm(CMA.ps)/CMA.chiN - 1);   % 计算步长的变化量

CMA.sigma = CMA.sigma * exp(min(delta_sigma,CMA.delta_sigma_max));   % 更新步长，限制其变化范围

CMA.C = triu(CMA.C) + triu(CMA.C,1)';    % 对称化协方差矩阵

[CMA.B,CMA.D] = eig(CMA.C);   % 重新计算特征分解

CMA.D = sqrt(diag(CMA.D))';  % 获取特征值的平方根

CMA = repairCMA(CMA);  % 修复 CMA-ES 参数（自定义函数）

CMA.invsqrtC = CMA.B * diag(CMA.D.^-1) * CMA.B';   % 重新计算协方差矩阵的逆平方根 

function [F,C] = LL_evaluate(xu,xl,fn)
% LL_evaluate 评估下层函数
% 输入:
%   xu - 上层决策变量
%   xl - 下层决策变量
%   fn - 函数句柄或函数名称
% 输出:
%   F - 下层函数值
%   C - 下层约束违反值



[F,~,C] = llTestProblem(xl,fn,xu); % 调用下层测试问题函数，获取函数值和约束

C = sum(max(0,C)); % 计算约束违反程度（累加所有正值）



function [F,C] = UL_evaluate(UPop,LPOP,fn)
% UL_evaluate 评估上层函数
% 输入:
%   UPop - 上层种群的上层决策变量数组
%   LPOP - 上层种群对应的下层种群的下层决策变量数组
%   fn   - 函数句柄或函数名称
% 输出:
%   F - 上层函数值
%   C - 上层约束违反值

[F,~,C] = ulTestProblem(UPop, LPOP, fn); % 调用上层测试问题函数，获取函数值和约束

C = sum(max(0,C));  % 计算约束违反程度（累加所有正值）


function Model = repairCMA(Model)
% repairCMA 修复和调整 CMA-ES 模型参数，确保协方差矩阵和步长的数值稳定性
% 输入:
%   Model - CMA-ES 参数结构体，包含字段 D（特征值向量）、C（协方差矩阵）、sigma（步长）、pc 和 ps（进化路径）
% 输出:
%   Model - 修复和调整后的 CMA-ES 参数结构体


dim = length(Model.D); % 获取模型的维度，即特征值向量 D 的长度
% limit condition of C to 1e14
  % 限制协方差矩阵 C 的条件数不超过 1e14

if any(Model.D<=0) % 检查是否存在非正的特征值
    Model.D(Model.D<0) = 0; % 将所有非正的特征值设为 0
    tmp = max(Model.D)/1e7;  % 计算一个调整因子 tmp，为当前最大特征值的 1e-7 倍
    Model.C = Model.C + tmp * eye(dim);  % 将 tmp 倍的单位矩阵加到协方差矩阵 C 上，提升特征值的最小值    
    Model.D = Model.D + tmp * ones(1,dim);  % 将 tmp 加到所有特征值 D 上，确保所有特征值为正
end


  % 确保协方差矩阵 C 的条件数（最大特征值与最小特征值之比）不超过 1e7
if max(Model.D) > 1e7 * min(Model.D)   % 检查最大特征值是否超过最小特征值的 1e7 倍
    tmp = max(Model.D)/1e7 - min(Model.D); % 计算需要增加的值 tmp，以缩小特征值的比例
    Model.C = Model.C + tmp * eye(dim);  % 将 tmp 倍的单位矩阵加到协方差矩阵 C 上，提升最小特征值
    Model.D = Model.D + tmp * ones(1,dim);  % 将 tmp 加到所有特征值 D 上，确保条件数不超过 1e7
end
% rescale sigma

% 重新调整步长 sigma，防止步长过大导致搜索不稳定
if Model.sigma > 1e7 * max(Model.D) % 检查步长 sigma 是否超过最大特征值的 1e7 倍
    fac = Model.sigma / max(Model.D); % 计算调整因子 fac，使得 sigma 不再过大
    Model.sigma = Model.sigma / fac; % 缩小步长 sigma
    Model.D = Model.D * fac; % 缩放特征值 D，以匹配步长的调整
    Model.pc = Model.pc * fac;  % 缩放进化路径 pc
    Model.C = Model.C * fac^2; % 缩放协方差矩阵 C，以保持协方差矩阵与特征值的一致性
end