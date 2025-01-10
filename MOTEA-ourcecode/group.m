%group
function team = group(pop, normal_obj)

 % 获取目标值矩阵的维度，Dim表示目标函数的维数，num_p表示个体数目
    [Dim,num_p]       = size(normal_obj); 

    % 提取种群中所有子种群的中心点，center是一个Dim x num_class的矩阵
    center            = [pop.center];

       % 获取子种群的数量，即中心点的数量
    num_class         = size(center,2);

     % 创建一个空的cell数组，用于存储分组后的个体
    team              = cell(num_class,1);     

     % 计算每个个体到所有子种群中心的距离
    dis               = center'*normal_obj;  

     % 对每个个体计算其与子种群中心的距离，并找到距离最小的中心
    [minval,minindex] = max(dis,[],1);    

     % 将每个个体分配到距离最小的子种群组中
    for i = 1:num_p
        team{minindex(i)}  = [team{minindex(i)},i]; % 根据minindex将个体索引添加到对应的组
    end   
end
