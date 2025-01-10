function SPOP = selection(TPOP,subpop,BI,Num,theta,perct)   
      LC   = [TPOP.LC]; % 获取种群中所有个体的可行性标志，LC为每个个体的约束违规程度
      W    = subpop.center;   % 获取当前子种群的中心权重
      weights = generateweight([TPOP.UX],BI);  % 调用 generateweight 函数计算所有个体的权重，传入 TPOP 中的UX和BI
      LF = (-[TPOP.LF]+subpop.elite.LF);  % 计算种群的局部适应度，LF为每个个体的局部适应度
      LF = LF./max(1,abs(subpop.elite.LF));  % 对局部适应度进行归一化
      isFeasible = LC == 0; % 判断每个个体是否满足约束条件，LC==0表示个体可行

       % 对于可行的个体，保持其LF不变；对于不可行的个体，调整其LF
      if any(isFeasible)
          LF(isFeasible)  = LF(isFeasible); % 可行个体的LF值不变
          LF(~isFeasible) = max(LF(isFeasible)) + LC(~isFeasible); % 不可行个体的LF值设为最大可行个体LF加上它的约束违规程度
      end


       % 如果LF的最大值减去最小值小于1，则进行调整；否则，按非线性方式转换LF值
     if max(LF)-min(0,min(LF))<1
         if min(LF)<0
             LF = LF -min(LF); % 如果最小LF小于0，则将LF平移至非负范围
         end
     else         
          LF = 1./(1+(1+3*perct)*exp(-LF)); % 如果LF的范围较大，则应用非线性变换
          if min(LF) < 1/(1+(1+3*perct))        
              LF = abs(LF - min([LF])); % 如果最小LF小于某个阈值，则将LF的最小值归零
          else
              LF = abs(LF - 1/(1+(1+3*perct)));  % 否则将LF的最小值标准化为一个较小的值
          end  
     end

     % 计算适应度的两个部分 d1 和 d2，用于选择个体
      d1 = (1+LF.^(1/2)).*(W'*weights); % d1为加权适应度
      d2 = (1+LF.^(1/2)).*sqrt(abs(1-(W'*weights).^2)); % d2为另一种加权适应度，用于平衡权重和局部适应度的偏差


       % 如果Num为1，使用 d1 和 d2 的平方根之和作为适应度
    if Num==1       
       fit_ = d1+theta*d2.^(1/2);       % 适应度计算公式，d1和d2根据theta调整    
    else

         % 如果d1小于1，则使用d1和d2的平方
        if all(d1<1)
            fit_ = d1+theta*d2.^(2); % 如果d1全都小于1，则d2的平方作为加权项
        else  
            fit_ = d1+theta*d2.^(2); % 否则，还是使用d2的平方作为加权项
        end
    end

      % 对所有个体按照适应度值进行排序，选择最好的Num个个体
        [val,index] = sort(fit_); % 按照适应度值升序排序
        SPOP = TPOP(index(1:Num)); % 选择适应度最好的Num个个体

%       if all(LC==0)
%           [val,index] = sort(fit_);
%           SPOP = TPOP(index(1:Num));      
%       else
%           isFeasible = LC == 0;
%           if any(isFeasible)
%               cfit_(isFeasible)  = fit_(isFeasible);
%               cfit_(~isFeasible) = max(fit_(isFeasible)) + LC(~isFeasible);
%               [val,index] = sort(cfit_);
%               SPOP = TPOP(index(1:Num));
%           else
%               cfit_ =  LC;
%               [val,index] = sort(cfit_);
%               SPOP = TPOP(index(1:Num));
%           end
%       end
end
