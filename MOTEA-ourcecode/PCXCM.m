function cld_x = PCXCM(parent1,parent2,parent3,g,d)%根据父代个体的基因进行交叉操作，生成一个子代个体

cld_x = parent1 + 0.1*(parent1-g)+ d/sum(abs(parent1-g))*(parent2-parent3)/2;