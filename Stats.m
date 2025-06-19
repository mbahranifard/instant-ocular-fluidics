%% outflow facility
x4 = [22.12, 34.32, 31.99, 21.64, 50.47];
% [S4,h4,stats4] = signrank([22.12, 34.32, 31.99])
pd4 = fitdist(transpose(x4),'normal');
% ci4 = paramci(pd4);
% ci4(:,1) - mean(x4')
[h4, p4] = ttest(x4,0,'Alpha',0.05)
% [hk4,pk4] = kstest(x4)
[hk4,pk4] = swtest(x4)

x5 = [20.51, 20.86, 24.1, 21.3, 35.2];
pd5 = fitdist(transpose(x5),'normal');
% [S5,h5,stats5] = signrank(x5)
[h5, p5] = ttest(x5,0,'Alpha',0.05)
[hk5,pk5] = swtest(x5)

x6 = [79.8, 79.15, 79.1, 86.23, 62.3];
pd6 = fitdist(transpose(x6),'normal');
[h6, p6] = ttest(x6,0,'Alpha',0.05)
[hk6,pk6] = swtest(x6)

x7 = [53.38, 58.99, 70.32, 73.37, 53.5 ];
pd7 = fitdist(transpose(x7),'normal');
[h7, p7] = ttest(x7,0,'Alpha',0.05)
[hk7,pk7] = swtest(x7)

%% Compliance
% Ctot = [67.7, 73.1, 59.1, 49.9, 52.5, 48.3, 35.9, 34.4, 71.75, 67.01, 67, 63.15] ;
% pdtot = fitdist(transpose(Ctot),'normal')
% [hc, pc] = ttest(Ctot,0,'Alpha',0.05);
% [hktot,pktot] = swtest(Ctot);
% [SC,hC,statsC] = signrank(Ctot)
% % [SC,hC,statsC] = ranksum(Ctot)

clear all
 Comp = readtable ('CompPlot.csv');
 Trans = Comp(string(table2cell(Comp(:,1)))=='Transient',2:3);
 Stead = Comp(string(table2cell(Comp(:,1)))=='Steady',2:3);
 
 
 for i=1:4
   StepStr = {'Step A','Step B','Step C','Step D'} ;
   StepCompT = table2array(Trans(string(table2cell(Trans(:,1)))==StepStr(i),2));
   StepCompS = table2array(Stead(string(table2cell(Stead(:,1)))==StepStr(i),2));
   [~,PT(i)] = swtest(StepCompT);
   [~,PS(i)] = swtest(StepCompS);
   [~,Ptot(i)]=ttest2(StepCompT,StepCompS)  
 end
 PT
 PS
 Ptot
 
 
 
 
 
 
 
 
 
