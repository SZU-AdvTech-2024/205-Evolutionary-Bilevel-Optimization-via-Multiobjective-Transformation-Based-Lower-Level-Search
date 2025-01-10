function [ulPop] = initialization(ulDim, ulDimMin, ulDimMax, llDimMin, llDimMax, ulPopSize, testProblemName)
    % 该函数用于初始化上层种群 ulPop，同时考虑约束条件。
    % 输入：
    %   ulDim: 上层维度的大小
    %   ulDimMin: 上层维度的最小边界
    %   ulDimMax: 上层维度的最大边界
    %   llDimMin: 下层维度的最小边界
    %   llDimMax: 下层维度的最大边界
    %   ulPopSize: 上层种群的大小
    %   testProblemName: 要求解的测试问题的名称
    % 输出：
    %   ulPop: 初始化后的上层种群



    significantDistance = 1e-6; % 定义显著距离，用于判断个体之间的距离是否足够大
    data = createData(ulDimMin, ulDimMax, llDimMin, llDimMax, testProblemName, []); % 创建问题相关数据
    if data.totalConstraints>0
        display('Ensuring feasibility for the constraints') % 如果存在约束，输出提示信息
        constrModel = createModelConstraints(data);  % 创建约束模型
    end
    
    for i=1:ulPopSize % 确保新生成的个体与已有个体之间的距离足够大
        if data.totalConstraints>0 && i<=ceil(ulPopSize/2)
            % 确保前半部分种群满足约束条件

            %Ensuring relaxed feasibility
            weights(1) = rand();  % 随机生成第一个权重
            weights(2) = 1-rand();  % 随机生成第二个权重，保证两个权重之和为 1
            
            ulPop(i,:) = ensureFeasibility(ulDimMin, ulDimMax, llDimMin, llDimMax, data, constrModel, weights, []);
            for j=1:i-1    % 确保新生成的个体与已有个体之间的距离足够大
                if computeDistance(ulPop(i,:), ulPop(j,:))<=significantDistance % 如果距离不够大，对个体进行变异
                    ulPop(i,:) = mutation(ulPop(i,:),1);
                    ulPop(i,:) = checkLimitsReflection(ulPop(i,:),ulDimMin,ulDimMax); % 检查边界并反射修正
                    break;
                end
            end
        
        
        
        
        else
            %Normal initialization % 正常初始化剩余种群个体
            ulPop(i,:) = ulDimMin + rand(1, ulDim).*(ulDimMax-ulDimMin);  % 在边界内随机生成个体
        end


    end

function data = createData(ulDimMin, ulDimMax, llDimMin, llDimMax, testProblemName, member)

     % 创建问题数据，包括上下层维度和种群初始化
    % 输入：
    %   ulDimMin, ulDimMax: 上层变量的最小和最大边界
    %   llDimMin, llDimMax: 下层变量的最小和最大边界
    %   testProblemName: 要求解的测试问题的名称
    %   member: 如果为空，则创建整个边界范围内的种群；如果不为空，则创建在该成员附近的种群

    %If the input parameter 'member' is empty then population is created in the entire
    %box constraint set. If the input parameter 'member' is not empty then
    %population is created around the 'member'

    dimMin = [ulDimMin llDimMin]; % 合并上下层最小边界
    dimMax = [ulDimMax llDimMax]; % 合并上下层最大边界
    ulDim = size(ulDimMin,2);  % 上层维度
    llDim = size(llDimMin,2); % 下层维度
    dim = ulDim + llDim;  % 总维度
    popSize = (dim+1)*(dim+2)/2+5*dim;  % 种群大小的计算方式     ?????___-----2024-11-16是否可改进


    if ~isempty(member) %___-----2024-11-16是否可改进
        % 如果给定了成员，则在该成员附近生成种群
        for i=1:popSize-1
            memberVariance = (dimMax - dimMin)/5;
            pop(i,:) = member + (sqrt(memberVariance)).*randn(size(member)); % 生成在成员附近的随机个体 假设 member 是一个长度为 5 的向量，randn(size(member)) 会生成一个长度为 5 的随机数向量 （均值为 0，方差为 1）的随机数。
        end
        pop(popSize,:) = member; % 最后一个个体为给定的成员本身
        pop=checkLimitsReflection(pop, dimMin, dimMax); % 检查种群是否在边界内，反射修正



        % 如果没有给定成员，则在整个边界范围内初始化种群
    else
        pop=initialize(popSize,dim,dimMin,dimMax); % 生成在成员附近的随机个体
        member = pop(popSize,:);  % 指定最后一个个体为成员
    end
    

     % 计算每个个体的上下层目标函数值和约束
    for i=1:popSize
        [ulFunctionVals(i,:), ulEquality, ulInequality]=ulTestProblem(pop(i,1:ulDim), pop(i,ulDim+1:end), testProblemName); %调用上层测试问题函数 ulTestProblem 计算第 i 个个体的上层目标函数值和约束条件。
        if ~isempty(ulEquality)  %如果上层等式约束不为空，将其存储在 ulEqualityConstrVals 中的第 i 行。
            ulEqualityConstrVals(i,:) = ulEquality;
        else
            ulEqualityConstrVals = []; %如果为空，则将 ulEqualityConstrVals 设为空，表示没有等式约束。
        end
        if ~isempty(ulInequality)
            ulInequalityConstrVals(i,:) = ulInequality;
        else
            ulInequalityConstrVals = [];
        end
        [llFunctionVals(i,:), llEquality, llInequality]=llTestProblem(pop(i,ulDim+1:end), testProblemName, pop(i,1:ulDim));
        if ~isempty(llEquality)
            llEqualityConstrVals(i,:) = llEquality;
        else
            llEqualityConstrVals = [];
        end
        if ~isempty(llInequality)
            llInequalityConstrVals(i,:) = llInequality;
        else
            llInequalityConstrVals = [];
        end
    end
    

    % 组织数据结构以存储所有约束和目标函数信息
    data.totalConstraints = size(ulEqualityConstrVals,2)+size(ulInequalityConstrVals,2)+size(llEqualityConstrVals,2)+size(llInequalityConstrVals,2);
    %size(ulEqualityConstrVals, 2)：返回上层等式约束的数量（列数）。
    %size(ulInequalityConstrVals, 2)：返回上层不等式约束的数量（列数）。
    %size(llEqualityConstrVals, 2)：返回下层等式约束的数量（列数）。
    %size(llInequalityConstrVals, 2)：返回下层不等式约束的数量（列数）。
    %所有这些数量相加得到 totalConstraints，即整个问题中的约束数量。



    data.ul.pop = pop(:,1:ulDim);  %将上层变量的种群存储在 data 结构体中
    data.ll.pop = pop(:,ulDim+1:end); %将下层变量的种群存储在 data 结构体中。
    
    data.ul.functionVals = ulFunctionVals;% 存储上层目标函数值
    data.ul.equalityConstrVals = ulEqualityConstrVals; % 存储上层等式约束值。
    data.ul.inequalityConstrVals = ulInequalityConstrVals; %存储上层不等式约束值。
     
    data.ll.functionVals = llFunctionVals;  % 存储下层目标函数值。
    data.ll.equalityConstrVals = llEqualityConstrVals;  %存储下层等式约束值。
    data.ll.inequalityConstrVals = llInequalityConstrVals; %存储下层不等式约束值。
    
function constrModel = createModelConstraints(data)


    % 根据数据创建约束模型，生成约束的二次逼近
    % 输入：
    %   data: 包含上下层种群和约束信息的数据结构
    % 输出：
    %   constrModel: 约束模型，包含等式和不等式约束的二次逼近



    pop = [data.ul.pop data.ll.pop]; %pop 是包含所有上层和下层种群的矩阵，通过将 data.ul.pop 和 data.ll.pop 拼接在一起得到。
    
    equalityConstrVals = [data.ul.equalityConstrVals, data.ll.equalityConstrVals];  % 将上层和下层的等式约束值拼接在一起，形成完整的等式约束集合。
    inequalityConstrVals = [data.ul.inequalityConstrVals, data.ll.inequalityConstrVals]; %将上层和下层的不等式约束值拼接在一起，形成完整的不等式约束集合。

    numEqualityConstr = size(equalityConstrVals,2); % 计算等式约束的数量，即 equalityConstrVals 中的列数。
    numInequalityConstr = size(inequalityConstrVals,2); % 计算不等式约束的数量，即 inequalityConstrVals 中的列数。
    
    if numEqualityConstr ~=0 % 对等式约束进行二次逼近： % 如果存在等式约束（numEqualityConstr 不为 0），则对每一个等式约束进行处理。
        for i=1:numEqualityConstr
            %将生成的二次逼近模型存储在 approx.equalityConstr{i} 中。
            approx.equalityConstr{i} = quadApprox(equalityConstrVals(:,i), pop); %调用 quadApprox 函数对第 i 个等式约束进行二次逼近，输入为对应的约束值（equalityConstrVals(:, i)）和种群变量（pop）
        
        
        
%             %----------------------------------------------------------------拟合过程演示代码--------------------------
%         
%             % 动态显示拟合过程
%             figure;
%             values = equalityConstrVals(:, i);
%             model = approx.equalityConstr{i};
%             scatter3(pop(:, 1), pop(:, 2), values, 'filled');
%             hold on;
%             [X1, X2] = meshgrid(linspace(min(pop(:, 1)), max(pop(:, 1)), 20), ...
%                                 linspace(min(pop(:, 2)), max(pop(:, 2)), 20));
%             Z = model.constant + model.linear(1) * X1 + model.linear(2) * X2 + ...
%                 model.sqmatrix(1, 1) * X1.^2 + model.sqmatrix(2, 2) * X2.^2 + ...
%                 2 * model.sqmatrix(1, 2) * X1 .* X2;
%             surf(X1, X2, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
%             title(['Equality Constraint ', num2str(i), ' - Dynamic Approximation']);
%             xlabel('x_1'); ylabel('x_2'); zlabel('Constraint Value');
%             hold off;
%             pause(1); % 暂停以创建动态效果
%             %----------------------------------------------------------------拟合过程演示代码--------------------------
%         
        end
    else
        approx.equalityConstr = []; % 如果不存在等式约束，则将 approx.equalityConstr 设为空数组。
    end

    if numInequalityConstr ~= 0 %如果存在不等式约束（numInequalityConstr 不为 0），则对每一个不等式约束进行处理。

        for i=1:numInequalityConstr
            approx.inequalityConstr{i} = quadApprox(inequalityConstrVals(:,i), pop);
        end
    else
        approx.inequalityConstr = [];
    end
    
    constrModel = approx;
    
function [ulMember, llMember] = ensureFeasibility(ulDimMin, ulDimMax, llDimMin, llDimMax, data, constrModel, weights, member)

    % 确保生成的个体满足约束条件
    % 输入：
    %   ulDimMin, ulDimMax: 上层变量的最小和最大边界
    %   llDimMin, llDimMax: 下层变量的最小和最大边界
    %   data: 包含上下层种群和约束信息的数据结构
    %   constrModel: 约束模型
    %   weights: 上下层目标函数的权重
    %   member: 起始点个体
    % 输出：
    %   ulMember, llMember: 满足约束条件的上层和下层成员

    %If input parameter 'member' is provided then the math program is
    %solved with the 'member' as starting point otherwise any random point
    %from data is chosen

    epsilonZero = 1e-6; 
    
    dimMin = [ulDimMin llDimMin];
    dimMax = [ulDimMax llDimMax];
    ulDim = size(ulDimMin,2);
    llDim = size(llDimMin,2);
    dim = ulDim + llDim;
        
    %Copying constraint model to approx
    approx = constrModel;
    
    %If weights are not provided then function to be optimized is taken as
    %0, which means that only feasibility wrt constraints is needed
    %otherwise a new function with upper and lower level weights is
    %constructed
    if isempty(weights)
        functionVals = zeros(popSize,1);
    else
        functionVals = weights(1)*data.ul.functionVals + weights(2)*data.ll.functionVals; %根据权重计算上层和下层目标函数的加权和。
    end
    
    pop = [data.ul.pop data.ll.pop];
    popSize = size(pop,1);
    if isempty(member)  % 如果没有指定 member，则从种群中随机选择一个个体作为起始点。
        r = ceil(rand*popSize);  %r 是用于从种群 pop 中随机选择一个个体 ceil向上取整
        member = pop(r,:);
    else
        functionVals = weights(1)*data.ul.functionVals + weights(2)*data.ll.functionVals; %如果指定了 member，则使用这个个体。
    end
        
    %Adding function model to approx
    approx.function = quadApprox(functionVals, pop); %使用 quadApprox 函数拟合目标函数值，并添加到约束模型 approx 中。

    options = optimset('Algorithm','interior-point','Display','off','TolX',1e-10,'TolFun',1e-10); %设置优化选项，包括使用内点算法、禁用显示以及设置收敛精度。
    
    if isLinear(approx,epsilonZero) %判断约束是否为线性。如果所有二次项系数接近于零，则可以认为是线性约束。 
        if ~isempty(approx.equalityConstr)%如果约束是线性的，使用 linprog 求解。
            for i=1:length(approx.equalityConstr)
                A_equality(i,:) = approx.equalityConstr{i}.linear;
                b_equality(i) = -approx.equalityConstr{i}.constant;%将线性约束的系数提取出来，构造等式约束矩阵 A_equality 和 b_equality
            end
        else
            A_equality = [];
            b_equality = [];
        end
        if ~isempty(approx.inequalityConstr)
            for i=1:length(approx.inequalityConstr)
                A_inequality(i,:) = approx.inequalityConstr{i}.linear;
                b_inequality(i) = -approx.inequalityConstr{i}.constant;
            end
        else
            A_inequality = [];
            b_inequality = [];
        end
        
        optionsLinprog = optimset('Display','off'); %使用 linprog 函数进行线性规划，求得最优个体 eliteIndiv。
        [eliteIndiv] = linprog(-approx.function.linear,A_inequality,b_inequality',A_equality,b_equality',dimMin,dimMax,member,optionsLinprog);
        eliteIndiv = eliteIndiv';
    else
        [eliteIndiv] = fmincon(@(x) -approximatedFunction(x,approx.function),member,[],[],[],[],dimMin,dimMax,@(x) approximatedConstraints(x,approx.equalityConstr,approx.inequalityConstr),options);
        %使用 approximatedConstraints 作为约束函数。
    end
    ulMember = eliteIndiv(1:ulDim);
    llMember = eliteIndiv(ulDim+1:end); %最终的最优个体被分为上层成员 ulMember 和下层成员 llMember
    
    
function approxFunctionValue = approximatedFunction(pop, parameters) %如果约束不是线性的，使用 fmincon 求解。

    approxFunctionValue = parameters.constant + pop*parameters.linear + pop*parameters.sqmatrix*pop';

function [c, ceq] = approximatedConstraints(pop, parametersEqualityConstr, parametersInequalityConstr)

    if ~isempty(parametersEqualityConstr)
        for i=1:length(parametersEqualityConstr)
            ceq(i) = parametersEqualityConstr{i}.constant + pop*parametersEqualityConstr{i}.linear + pop*parametersEqualityConstr{i}.sqmatrix*pop';
        end
    else
        ceq = [];
    end
    
    if ~isempty(parametersInequalityConstr)
        for i=1:length(parametersInequalityConstr)
            c(i) = parametersInequalityConstr{i}.constant + pop*parametersInequalityConstr{i}.linear + pop*parametersInequalityConstr{i}.sqmatrix*pop';
        end
    else
        c = []; %Negative value suggests that there is no active inequality
    end
    
function bool = isLinear(approximateFunctionParameters,epsilonZero)

    bool = true;
    if sum(sum(abs(approximateFunctionParameters.function.sqmatrix)>epsilonZero))>0
        bool = false;
    end
    for i=1:length(approximateFunctionParameters.equalityConstr)
        if sum(sum(abs(approximateFunctionParameters.equalityConstr{i}.sqmatrix)>epsilonZero))>0
            bool = false;
        end
    end
    for i=1:length(approximateFunctionParameters.inequalityConstr)
        if sum(sum(abs(approximateFunctionParameters.inequalityConstr{i}.sqmatrix)>epsilonZero))>0
            bool = false;
        end
    end
    
    
function pop=initialize(popSize,dim,dimMin,dimMax)

    for i=1:popSize
        pop(i,:) = dimMin + rand(1, dim).*(dimMax-dimMin);
    end
    
function d = computeDistance(matrix1, matrix2)
 
    %Computes pairwise distance between rows of matrix1 and matrix2
    sz1 = size(matrix1, 1);
    sz2 = size(matrix2, 1);
    
    for i = 1:sz1
        for j = 1:sz2
            d(i,j) = sqrt(sum((matrix1(i,:)-matrix2(j,:)).^2));
        end
    end
    
    
function members = mutation(members, probMutation)
 
    numOffsprings=size(members,1);
    ulDim=size(members,2);
    mum=20;
    for i=1:numOffsprings
        r = rand(1,ulDim);
        index = r<0.5;
        delta(index) = (2*r(index)).^(1/(mum+1)) - 1;
        index = r>=0.5;
        delta(index) = 1 - (2*(1 - r(index))).^(1/(mum+1));
 
        % Generate the corresponding child element.
        r = rand(1,ulDim);
        if ~all(r>=probMutation)
            members(i,r<probMutation) = members(i,r<probMutation) + delta(r<probMutation);
        end
    end
    
function members=checkLimitsReflection(members, dimMin, dimMax)
    %This function reflects an infeasible point into the variable bounds. If the
    %point lies far away, it assigns it a random position in the bounds.
    numOfMembers = size(members,1);
    dimMinMatrix = dimMin(ones(1,numOfMembers),:);
    dimMaxMatrix = dimMax(ones(1,numOfMembers),:);
    i = 0;
    while sum(sum(members<dimMinMatrix)) || sum(sum(members>dimMaxMatrix))
        I = members<dimMinMatrix-(dimMaxMatrix-dimMinMatrix);
        J = members>dimMaxMatrix+(dimMaxMatrix-dimMinMatrix);
        randI = rand(size(I));
        randJ = rand(size(J));
        members(I) = dimMinMatrix(I) + randI(I).*(dimMaxMatrix(I)-dimMinMatrix(I));
        members(J) = dimMinMatrix(J) + randJ(J).*(dimMaxMatrix(J)-dimMinMatrix(J));
        members(members<dimMinMatrix)=members(members<dimMinMatrix) + 2*(dimMinMatrix(members<dimMinMatrix)-members(members<dimMinMatrix));
        members(members>dimMaxMatrix)=members(members>dimMaxMatrix) + 2*(dimMaxMatrix(members>dimMaxMatrix)-members(members>dimMaxMatrix));
    end
    
