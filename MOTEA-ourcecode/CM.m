function cld_x = CM(SQ,i,subpop,Num,BI)  
         % 获取当前个体 i 的邻居集和已选择的邻居集
        Nbs = setdiff(subpop(i).neighbours,subpop(i).snbs);

        % 计算群体中心点的数量
        N   = length([subpop.center]);      

            % 计算当前个体 i 的可选候选集 WQ（不包含已选择邻居和邻居）
        WQ  = setdiff([1:N],[subpop(i).snbs,subpop(i).neighbours]);     

         % 获取当前个体 i 的交互信息
        Inter = [subpop(i).inter];


        % 从交互信息中选择一个父母个体 parent1
        parent1 = Inter(SQ).X;    


        % 获取变量的上下限
        x_min = [BI.ux_lb,BI.l_lb]';
        x_max = [BI.ux_ub,BI.l_ub]';


         % 获取父母个体的维度
        d = length(parent1);   
        %*****************************


         % 随机选择父母的策略
        if (rand < 0.8 || isempty(WQ))  

            % 如果随机数小于 0.8，或者没有可用的候选集 WQ，则选择当前个体的邻居或交互集中的父母
%             r = unidrnd(length(Inter));
            r = mod(SQ,2)+1; % 根据 SQ 的奇偶性选择另一个父母的交互信息
            parent2 = Inter(r).X;


            % 如果 parent1 和 parent2 非常相似（差异小于 1e-7），且邻居集 Nbs 不为空，则选择一个新的父母
            if max(abs(parent1-parent2))<1e-7 &&~isempty(Nbs) 

                % 随机选择一个邻居并从该邻居的交互信息中选择 parent2
                r2 = floor(rand*(size(Nbs,2)))+1;



                Inter = [subpop(Nbs(r2)).inter];


                r = unidrnd(length(Inter));  % 随机选择交互信息中的一个个体


                parent2 = Inter(r).X;

                 % 如果 parent1 和 parent2 依然相似，则选择邻居的最优解（elite）
                if max(abs(parent1-parent2))<1e-7    
                   parent2 = subpop(Nbs(r2)).elite.X;
                end
            end
        else

            % 如果随机数大于 0.8，选择 WQ 中的候选个体
            r3 = floor(rand*(size(WQ,2)))+1;


            Inter = [subpop(WQ(r3)).inter];  % 从 WQ 中选择一个个体


            r = unidrnd(length(Inter));  % 随机选择交互信息中的一个个体
            parent2 = Inter(r).X;


             % 如果 parent1 和 parent2 非常相似，则选择 WQ 的最优解（elite）
             if max(abs(parent1-parent2))<1e-7   
                 parent2 = subpop(WQ(r3)).elite.X;
             end
        end

         % 如果 parent1 和 parent2 差异非常小，则选择当前个体的最优解（elite）作为 parent1
         if max(abs(parent1-parent2))<min(1e-7,max(abs(subpop(i).elite.X-parent2)))  
          parent1 = subpop(i).elite.X; 
         end           

         % 调用 SBXCM 函数进行模拟二进制交叉（SBX），生成新的子代 cld_x
         cld_x = SBXCM(parent1,parent2,x_min,x_max,BI);    
end

