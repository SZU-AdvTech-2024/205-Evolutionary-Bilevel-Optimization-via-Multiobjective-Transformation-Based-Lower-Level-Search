function [upperMember lowerMember upperFunctionValue lowerFunctionValue]=getOptimalSolutionSMD(upperDimensions,lowerDimensions,testProblemName)
    
     % 主函数，用于获取 SMD 问题的最优解
    % 输入:
    %   upperDimensions: 上层变量的维度数量
    %   lowerDimensions: 下层变量的维度数量
    %   testProblemName: 要解决的 SMD 测试问题的名称
    % 输出:
    %   upperMember: 上层变量的最优值
    %   lowerMember: 下层变量的最优值
    %   upperFunctionValue: 上层目标函数的值
    %   lowerFunctionValue: 下层目标函数的值



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Function call here
     % 调用对应测试问题函数的函数句柄
    fhandle = str2func(testProblemName);  %其实就是拿到函数接口
    [upperMember lowerMember] = fhandle(upperDimensions,lowerDimensions);% 调用本文件下的具体的 SMD 函数 如function [upperMember lowerMember] = smd12(upperDimensions,lowerDimensions)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

     % 计算下层目标函数和约束
    [lowerFunctionValue lowerEqualityConstrVals lowerInequalityConstrVals]=llTestProblem(lowerMember, testProblemName, upperMember);
    % 计算上层目标函数和约束
    [upperFunctionValue upperEqualityConstrVals upperInequalityConstrVals]=ulTestProblem(upperMember, lowerMember, testProblemName);


    % 不同 SMD 问题的函数定义
function [upperMember lowerMember] = smd1(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2); %上层和下层共享的变量数量
    p = upperDimensions - r; % 仅上层的变量数量
    q = lowerDimensions - r; % 仅下层的变量数量

    upperMember = zeros(1,upperDimensions);% 初始化上层变量
    lowerMember = zeros(1,lowerDimensions);% 初始化下层变量
    
function [upperMember lowerMember] = smd2(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = [zeros(1,q) ones(1,r)];
    
function [upperMember lowerMember] = smd3(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = zeros(1,lowerDimensions);
    
function [upperMember lowerMember] = smd4(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = zeros(1,lowerDimensions);
    
function [upperMember lowerMember] = smd5(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = [ones(1,q) zeros(1,r)];
    
function [upperMember lowerMember] = smd6(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = floor((lowerDimensions - r)/2 - eps);
    s = ceil((lowerDimensions - r)/2 + eps);

    upperMember = zeros(1,upperDimensions);
    lowerMember = zeros(1,lowerDimensions);

function [upperMember lowerMember] = smd7(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = [zeros(1,q) ones(1,r)];

function [upperMember lowerMember] = smd8(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = [ones(1,q) zeros(1,r)];

function [upperMember lowerMember] = smd9(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = zeros(1,lowerDimensions);

function [upperMember lowerMember] = smd10(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = 1/sqrt(p+r-1)*ones(1,upperDimensions);
    lowerMember = [1/sqrt(q-1)*ones(1,q) atan(1/sqrt(p+r-1)*ones(1,r))];

function [upperMember lowerMember] = smd11(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = zeros(1,upperDimensions);
    lowerMember = [zeros(1,q) exp(-1/sqrt(r))*ones(1,r)];

function [upperMember lowerMember] = smd12(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = 1/sqrt(p+r-1)*ones(1,upperDimensions);  %数乘向量获得最优解
    lowerMember = [1/sqrt(q-1)*ones(1,q) atan(1/sqrt(p+r-1)-1/sqrt(r))*ones(1,r)];
    
function [upperMember lowerMember] = smd13(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = lowerDimensions - r;

    upperMember = [ones(1,p) zeros(1,r)];
    lowerMember = [zeros(1,p) ones(1,r)];

function [upperMember lowerMember] = smd14(upperDimensions,lowerDimensions)

    r = floor(upperDimensions/2);
    p = upperDimensions - r;
    q = floor((lowerDimensions - r)/2 - eps);
    s = ceil((lowerDimensions - r)/2 + eps);

    upperMember = [ones(1,p) zeros(1,r)];
    lowerMember = zeros(1,lowerDimensions);