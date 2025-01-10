function subpop = select(subpop,LPOP,team,BI,Num,theta,perct)
 num_class = length(team); % 获取子种群数量（即子区域数）
for i=1:num_class % 遍历每个子种群
    num_p1            = length([team{i}]); % 当前子种群的个体数
    num_p             = num_p1 + length([team{subpop(i).neighbours}]);  % 当前子种群的个体数加上与邻居子种群合并后的个体数

    if num_p<Num % 如果子种群总个体数小于 Num（最大个体数）
           subpop(i).inter    = selection(LPOP,subpop(i),BI,Num,1e5,perct); % 从LPOP中选择Num个个体
           subpop(i).temp = selection(subpop(i).inter,subpop(i),BI,1,theta,perct);   % 从inter子种群中再选择1个个体作为临时个体
    else % 如果子种群总个体数大于等于Num
        if num_p1<Num+1  % 如果当前子种群的个体数小于Num+1
            selind    = LPOP([team{[i,subpop(i).neighbours]}]);   % 从当前子种群和邻居子种群中选择个体
            subpop(i).inter  = selection(selind,subpop(i),BI,Num,theta,perct); % 从合并的个体集合中选择Num个个体
            subpop(i).temp = selection(subpop(i).inter,subpop(i),BI,1,theta,perct);         % 从inter子种群中选择1个个体作为临时个体 
        else % 如果当前子种群的个体数大于等于Num+1
            selind    = LPOP([team{i}]); % 只从当前子种群选择个体
            subpop(i).inter = selection(selind,subpop(i),BI,Num,theta,perct);  % 从当前子种群中选择Num个个体
            subpop(i).temp = selection(subpop(i).inter,subpop(i),BI,1,theta,perct);   % 从inter子种群中选择1个个体作为临时个体
        end
    end
end
end