% run MOTEA----source code------

addpath(genpath(pwd));
maxRuns = 19;
BI_list = [];
% BI_list = [BI_list;getBLOPinfo('SMD',2:2,20)];
%BI_list = [BI_list;getBLOPinfo('TP',1:8)];
BI_list = [BI_list;getBLOPinfo('SMD',12:12,10)];
% BI_list = [BI_list;getBLOPinfo('SMD',1:12,20)]; 
% BI_list = [BI_list;getBLOPinfo('GoldMini ng')];
% BI_list = [BI_list;getBLOPinfo('DecisionMaking')];


algName = 'MOTEA';
y1=[];y2=[];y3=[];y4=[];y5=[];
 i = 1;
for BI = BI_list'
    if BI.dim == 3
        % for gold mining problem
        BI.UmaxFEs = 1500;
        BI.UmaxImprFEs = 70;
        BI.LmaxFEs = 150;
        BI.LmaxImprFEs = 15;
    elseif BI.dim == 5 || BI.dim == 6 || strcmp(BI.fn(1:2),'tp')
        % for generic benchmark test problem or decision making problem 
        BI.UmaxFEs = 2500;
        BI.UmaxImprFEs = 350;
        BI.LmaxFEs = 250;
        BI.LmaxImprFEs = 25;
    elseif BI.dim == 10
        BI.UmaxFEs = 3500;
        BI.UmaxImprFEs = 500;
        BI.LmaxFEs = 350;
        BI.LmaxImprFEs = 35;
    elseif BI.dim == 20
        BI.UmaxFEs = 5000;
        BI.UmaxImprFEs = 750;
        BI.LmaxFEs = 500;
        BI.LmaxImprFEs = 50;
    else
        error('unknown dimensionality');
    end
    
    BI.u_N = 50;
    BI.l_N = 50;


    x1=[];x2=[];x3=[];x4=[];x5=[];
  seed = 68;randn('state',seed);rand('state',seed); 
  % 	parfor runNo = 1:maxRuns
    for runNo = 1:maxRuns
		tic;
        ins = MOTEA(BI);
		ins.runPath = pwd;
		ins.runNo = runNo;
		ins.BI = BI;
        ins.BI.fn = strvcat(ins.BI.fn);
		ins.alg = algName;
		ins.runTime = toc;
        fprintf('%s %s #%d [%g,%g]\n', ins.alg, ins.BI.fn, ins.runNo, ins.UF, ins.LF);    
    end	   
end  
   