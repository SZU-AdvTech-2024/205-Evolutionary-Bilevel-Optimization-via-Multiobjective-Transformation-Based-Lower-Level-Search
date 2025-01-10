% 完整的测试代码，用于调试 createModelConstraints 函数并演示拟合过程

% 假设上层和下层变量的数量以及种群大小
ulDim = 2;
llDim = 2;
ulPopSize = 3;
ulDimMin = [0, 0];
ulDimMax = [5, 5];
llDimMin = [0, 0];
llDimMax = [5, 5];

% 创建测试数据
data.ul.pop = [1.0, 2.0; 1.5, 2.5; 2.0, 3.0];
data.ll.pop = [3.0, 4.0; 3.5, 4.5; 4.0, 5.0];

data.ul.equalityConstrVals = [0.1; -0.2; 0.0];
data.ul.inequalityConstrVals = [0.3; 0.1; 0.2];

data.ll.equalityConstrVals = [-0.1; 0.2; -0.1];
data.ll.inequalityConstrVals = [0.4; 0.3; 0.5];

% 调用 createModelConstraints 函数创建约束模型
constrModel = createModelConstraints(data);

% 打印约束模型以供调试
disp('Constraint Model:');
disp(constrModel);

% 画图以演示拟合过程
figure;
for i = 1:size(data.ul.equalityConstrVals, 2)
    % 获取原始数据
    values = data.ul.equalityConstrVals(:, i);
    
    % 生成拟合模型
    model = quadApprox(values, [data.ul.pop data.ll.pop]);
    
    % 绘制原始数据点
    subplot(1, size(data.ul.equalityConstrVals, 2), i);
    scatter3(data.ul.pop(:, 1), data.ul.pop(:, 2), values, 'filled');
    hold on;
    
    % 绘制拟合曲面
    [X1, X2] = meshgrid(linspace(min(data.ul.pop(:, 1)), max(data.ul.pop(:, 1)), 20), ...
                        linspace(min(data.ul.pop(:, 2)), max(data.ul.pop(:, 2)), 20));
    
    % 计算拟合曲面上的值
    Z = model.constant + model.linear(1) * X1 + model.linear(2) * X2 + ...
        model.sqmatrix(1, 1) * X1.^2 + model.sqmatrix(2, 2) * X2.^2 + ...
        2 * model.sqmatrix(1, 2) * X1 .* X2;
    
    % 绘制曲面
    surf(X1, X2, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    title(['Equality Constraint ', num2str(i)]);
    xlabel('x_1'); ylabel('x_2'); zlabel('Constraint Value');
    hold off;
end

% 调整画图布局
sgtitle('Quadratic Approximation of Equality Constraints');

% 动态演示拟合过程
figure;
for t = 1:100
    % 获取动态数据（模拟噪声或扰动）
    noise = 0.05 * randn(size(data.ul.equalityConstrVals));
    values = data.ul.equalityConstrVals + noise;
    
    % 生成拟合模型
    model = quadApprox(values, [data.ul.pop data.ll.pop]);
    
    % 清除当前图像
    clf;
    
    % 绘制原始数据点
    scatter3(data.ul.pop(:, 1), data.ul.pop(:, 2), values, 'filled');
    hold on;
    
    % 绘制拟合曲面
    [X1, X2] = meshgrid(linspace(min(data.ul.pop(:, 1)), max(data.ul.pop(:, 1)), 20), ...
                        linspace(min(data.ul.pop(:, 2)), max(data.ul.pop(:, 2)), 20));
    
    % 计算拟合曲面上的值
    Z = model.constant + model.linear(1) * X1 + model.linear(2) * X2 + ...
        model.sqmatrix(1, 1) * X1.^2 + model.sqmatrix(2, 2) * X2.^2 + ...
        2 * model.sqmatrix(1, 2) * X1 .* X2;
    
    % 绘制曲面
    surf(X1, X2, Z, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
    title(['Dynamic Quadratic Approximation (Iteration ', num2str(t), ')']);
    xlabel('x_1'); ylabel('x_2'); zlabel('Constraint Value');
    hold off;
    
    % 暂停以创建动态效果
    pause(0.1);
end

% 定义 quadApprox 函数，用于生成二次逼近模型
function model = quadApprox(values, pop)
    % 确定种群大小和维度
    [numSamples, numFeatures] = size(pop);
    
    % 构建设计矩阵 X
    X = ones(numSamples, 1 + numFeatures + numFeatures * (numFeatures + 1) / 2);
    X(:, 2:1+numFeatures) = pop; % 线性项

    % 构建二次项
    index = 1 + numFeatures;
    for i = 1:numFeatures
        for j = i:numFeatures
            X(:, index) = pop(:, i) .* pop(:, j);
            index = index + 1;
        end
    end

    % 使用最小二乘法拟合 values
    coef = X \ values;

    % 将拟合得到的系数分解到模型中
    model.constant = coef(1);
    model.linear = coef(2:1+numFeatures);
    
    % 构建二次项矩阵 Q
    model.sqmatrix = zeros(numFeatures, numFeatures);
    index = 1 + numFeatures;
    for i = 1:numFeatures
        for j = i:numFeatures
            model.sqmatrix(i, j) = coef(index);
            model.sqmatrix(j, i) = coef(index); % 保证 Q 是对称的
            index = index + 1;
        end
    end
end

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
end
