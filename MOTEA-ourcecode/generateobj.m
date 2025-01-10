function objv=generateobj(LPOP,subpop,BI)
num = length(LPOP); % 获取 LPOP 数组的长度，即个体的数量
Temp = generateweight([LPOP.UX],BI);    % 调用 generateweight 函数，生成一个权重矩阵 Temp
LF = -[LPOP.LF]; % 获取个体的 LF 值，并取其负值
% if -subpop.elite.LF < min(LF)
%     disp('error');
% end


 % 如果 subpop.elite 存在，则调整 LF
if ~isempty(subpop.elite)
    LF = LF - min([-subpop.elite.LF,LF]);   % 调整 LF，使其基于最优个体（elite）
else
    LF = LF - min([LF]);  % 否则，直接调整 LF
end


 % 对每个个体进行目标值计算
for j = 1 : num
%     indiv = LPOP(j);
%     x = [(indiv.UX-BI.u_lb)./(BI.u_ub-BI.u_lb),0];
%     M = BI.u_dim+1;
%     cofficentMatrix = tril(ones(M,M),-1);
%     X =cos( pi/2*repmat(x(1:M),M,1));
%         X = X.*cofficentMatrix;
%         for i = 1:M-1
%             X(i,i) = sin(pi/2*x(i));
%         end
%         X(M,M) = 1;
%         X = X+(tril(ones(M,M),-1))';
%         y = cumprod(X,2);
% %     X = repmat(x(1:M),M,1);
% %     X = X.*cofficentMatrix;
% %     for i = 1:M-1
% %         X(i,i) = 1-x(i);
% %     end
% %     X(M,M) = 1;
% %     X = X+(tril(ones(M,M),-1))';
% %     y = cumprod(X,2);
%     objv(:,j) = (2-1/(1+log(1+abs(LF(j)))))*y(:,end);
% %     objv(:,j) = (1+1/(1+exp(-LF(j))))*y(:,end);
% objv(:,j) = (2-1/(1+log(1+abs(LF(j)))))*W(:,j);

   % 通过公式计算目标值
objv(:,j) = (1/2+1/(1+exp(-LF(j))))*Temp(:,j);
end

end