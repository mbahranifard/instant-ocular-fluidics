function[fig1,fig2,fig3,fig4,fig5,fig6,fig7,fig8,fig15,fig16,C10,CIOP,sC10,sCIOP,Pf] = DirectComp (stepP,stepQ,stepPa,steptmin,actoff,Cc,FSys,Yr,B,PrC,Cj,Cq,KeepStepsUp,gracep,PLOT)

NStep = size(stepP,2);

%% for test
FSys=3.5;
%% %%%%%%%%%%%%%%%%%%%%%
%% overall figure

%%
fig51 = figure(51)
toverall = (steptmin(:,4:7)-steptmin(1,4))*60;
% toverall = steptmin;
% for i=2:length(KeepStepsUp)
for i=4:7
    if i~=4:7
        ax1=subplot(2,1,1)
        hold all
        plot(toverall(:,i) ,stepQ(:,i),'color',[1 1 1]*.3,'LineWidth',2)

        ax2=subplot(2,1,2)
        hold on
        plot(toverall(:,i), stepP(:,i),'color',[1 1 1]*.3,'LineWidth',2)
    else
%         ax1=subplot(2,1,1)
%         hold all
%         plot(toverall(:,i) ,stepQ(:,i),'color','b','LineWidth',2)
% 
%         ax2=subplot(2,1,2)
%         hold on
%         plot(toverall(:,i), stepP(:,i),'color','b','LineWidth',2)
        ax1=subplot(2,1,1)
        hold all
        plot(toverall(:,i-3) ,stepQ(:,i),'color','b','LineWidth',2)

        ax2=subplot(2,1,2)
        hold on
        plot(toverall(:,i-3), stepP(:,i),'color','b','LineWidth',2)
    end
end
%%

% for target = KeepStepsUp
% for target = 5
for target = 8
    
if target == 1      %this is because the StepVars function doens't give data for step 1. Needs to be fixed     
    continue
end

kinterval = [actoff(target)+10:actoff(target)+210]; %% for Liz's data

% kinterval = [actoff(target):170];
% if  numel(kinterval) <15    %% un-comment for your own data
%     continue
% end

i = find(KeepStepsUp == target);    
P = rmmissing(stepP (:,target));
t = rmmissing(steptmin (:,target));
Q = rmmissing(stepQ (:,target));
Pa = rmmissing(stepPa (:,target));
%% Filter P and Q
dt = diff(t(1:2));
Pnofilt = P;
iSG = round(0.005/dt); %it was 0.005
P = sgolayfilt(P,1,iSG*2+1,[],1);

Qnofilt = Q;
iSG = round(0.0025/dt);
% Q = sgolayfilt(Q,1,iSG*2+1,[],1);
gaussFilter = gausswin(2*iSG + 1)';
gaussFilter = gaussFilter / sum(gaussFilter); % normalize
Q = conv(Q, gaussFilter, 'same');

%% Q calculated by Cq
% Qref = Cq .* (Pa-Pnofilt);
% Qrefnofilt = Qref;
% iSG = round(0.005/dt);
% Qref = sgolayfilt(Qref,1,iSG*2+1,[],1);

%% P first derivative
derP = diff(P)/dt;
derPnofilt = derP;
iSG = round(0.0025/dt);
derP = sgolayfilt(derP,1,iSG*2+1,[],1);

[derP3,derP4,der2P3,der2P5] = ThreePointDiff(P,dt);
derP3nofilt = derP3; derP4nofilt = derP4; der2P3nofilt = der2P3; der2P5nofilt = der2P5; 
derP3 = sgolayfilt(derP3,1,iSG*2+1,[],1); derP4 = sgolayfilt(derP4,1,iSG*2+1,[],1); der2P3 = sgolayfilt(der2P3,1,iSG*2+1,[],1); der2P5 = sgolayfilt(der2P5,1,iSG*2+1,[],1);

% gaussFilter = gausswin(2*iSG + 1)';
% gaussFilter = gaussFilter / sum(gaussFilter); % normalize
% derP = conv(derP, gaussFilter, 'same');

%%
% figure;
% hold on
% plot(t,derP3)
% [~,index] = findpeaks(derP3,'MinPeakDistance',length(t)-10);
% plot(t,derP3)
% plot(t(index),derP3(index),'og','MarkerSize',8)

%% P second derivative
der2P = diff(derP)/dt;
der2Pnofilt = der2P;
iSG = round(0.0025/dt);
der2P = sgolayfilt(der2P,1,iSG*2+1,[],1);
% gaussFilter = gausswin(2*iSG + 1)';
% gaussFilter = gaussFilter / sum(gaussFilter); % normalize
% der2P = conv(der2P, gaussFilter, 'same');

%% flow rate derivative
derQ = diff(Q)/dt;
derQnofilt = derQ;
iSG = round(0.0025/dt);
derQ = sgolayfilt(derQ,1,iSG*2+1,[],1);

[derQ3,derQ4,~,~] = ThreePointDiff(Q,dt);
derQ3nofilt = derQ3; derQ4nofilt = derQ4;
derQ3 = sgolayfilt(derQ3,1,iSG*2+1,[],1); derQ4 = sgolayfilt(derQ4,1,iSG*2+1,[],1); 

% gaussFilter = gausswin(2*iSG + 1)';
% gaussFilter = gaussFilter / sum(gaussFilter); % normalize
% derQ = conv(derQ, gaussFilter, 'same');

%% k2 and k1
k2 = Q(1:end-1)./derP - FSys;
k1 = P(1:end-1)./derP;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THIS IS A TEST %%%%%%%%%%
% k2 = Qref(1:end-1)./derP - FSys;
% k1 = P(1:end-1)./derP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k1nofilt = k1;
iSG = round(0.005/dt); % was 0.01
k1 = sgolayfilt(k1,1,iSG*2+1,[],1);

k2nofilt = k2;
iSG = round(0.005/dt); % was 0.01
k2 = sgolayfilt(k2,1,iSG*2+1,[],1);

[cf,rstat] = fit(k1(kinterval),k2(kinterval),'poly1'); 
kfit = coeffvalues(cf);
C10 (target) = kfit (1); 
cf_confint = confint(cf);
a = kfit(1);
% b = cf_coeff(2);
a_uncert = (cf_confint(2,1) - a); %% 95% confidence interval
sC10 (target) = a_uncert;
% b_uncert = (cf_confint(2,2) - cf_confint(1,2))/2;
% rcoeff = corrcoef(k1(kinterval), k2(kinterval));
rval =['R^2= ',num2str(rstat.rsquare)];
if sign(max(k1(kinterval))) == -1 
k1val = [min(k1(kinterval))*1.5 max(k1(kinterval))*0.5];
elseif sign(min(k1(kinterval))) == 1
k1val = [min(k1(kinterval))*0.5 max(k1(kinterval))*1.5];
else
k1val = [min(k1(kinterval))*1.5 max(k1(kinterval))*1.5];
end
k2val = feval(cf,k1val);


%% IOP calculations 
% IOP = P - Q/Cc;
IOP = P(1:end-1) - Q(1:end-1)/Cc + FSys/Cc*derP;
% IOP = P - Q/Cc + (FSys/Cc)*derP3;

%diagnositc
% figure(130)
% hold on
% plot(t,P,'.-')
% yyaxis right
% 
% plot(t,Q,'.-')
% yyaxis left

IOPnofilt = IOP;
iSG = round(0.005/dt); % was 0.01
IOP = sgolayfilt(IOP,1,iSG*2+1,[],1);

derIOP = derP(1:end-1) - derQ(1:end-1)/Cc + (FSys/Cc)*der2P;
% derIOP = derP3 - derQ3/Cc + (FSys/Cc)*[der2P; 0; 0];
derIOPnofilt = derIOP;
iSG = round(0.001/dt); % was 0.01
derIOP = sgolayfilt(derIOP,1,iSG*2+1,[],1);
% gaussFilter = gausswin(2*iSG + 1)';
% gaussFilter = gaussFilter / sum(gaussFilter); % normalize
% derIOP = conv(derIOP, gaussFilter, 'same');

% k1N = P(1:end-1)./derIOP;
 k1N = IOP(1:end-1)./derIOP;
%   k1N = IOP./derIOP;
k2N = (Q(2:end-1)-FSys*derP(1:end-1))./derIOP;
% k2N = (Q-FSys*derP3)./derIOP;

k1Nnofilt = k1N;
iSG = round(0.001/dt); % was 0.01
k1N = sgolayfilt(k1N,1,iSG*2+1,[],1);

% sigma = 10; % pick sigma value for the gaussian
% gaussFilter = gausswin(2*iSG + 1)';
% gaussFilter = gaussFilter / sum(gaussFilter); % normalize
% k1N = conv(k1N, gaussFilter, 'same');

k2Nnofilt = k2N;
iSG = round(0.001/dt); % was 0.01
k2N = sgolayfilt(k2N,1,iSG*2+1,[],1);
% gaussFilter = gausswin(2*iSG + 1)';
% gaussFilter = gaussFilter / sum(gaussFilter); % normalize
% k2N = conv(k2N, gaussFilter, 'same');

% kinterval = [actoff(target):170];
IOPfit = polyfit(k1N(kinterval),k2N(kinterval),1);

[cf,rstat3] = fit(k1N(kinterval),k2N(kinterval),'poly1'); 
IOPfit= coeffvalues(cf);
CIOP (target) = IOPfit (1); 
cf_confint = confint(cf);
a = IOPfit(1);
% b = cf_coeff(2);
a_uncert = (cf_confint(2,1) - a);
sCIOP (target) = a_uncert;
% b_uncert = (cf_confint(2,2) - cf_confint(1,2))/2;
rval3 =['R^2= ',num2str(rstat3.rsquare)];
if sign(max(k1N(kinterval))) == -1 
k1Nval = [min(k1N(kinterval))*1.5 max(k1N(kinterval))*0.5];
elseif sign(min(k1N(kinterval))) == 1
k1Nval = [min(k1N(kinterval))*0.5 max(k1N(kinterval))*1.5];
else
k1Nval = [min(k1N(kinterval))*1.5 max(k1N(kinterval))*1.5];
end
k2Nval = feval(cf,k1Nval);

%% Instant Compliance and outflow facility
% FEye10 = k2 - k1 * C10 (target);
% FEyeIOP = k2N - k1N * CIOP(target);
% FEye10 = k2 - k1 .* 11.3 .* ((P(1:end-1)/8).^(-1.1));
% FEye13 = k2 - k1 * C13 (target);
Pf (target) = mean(P(kinterval));

InstOut = (k2N-IOPfit(2))./k1N;
InstComp = k2N-IOPfit(1)*k1N;


%% calculate C and Phi for different intervals (Incremental)
shifter = 50;
shifterfull = 3000;
InterLength = 100;
for L=0:(shifterfull-InterLength)/shifter    
    NewInterval = [kinterval(1)+L*shifter:kinterval(1)+L*shifter+InterLength];
%     IOPfit{1} = polyfit(k1N(NewInterval),k2N(NewInterval),1);
    [cf,rstat4] = fit(k1N(NewInterval),k2N(NewInterval),'poly1'); 
    IOPfitinc= coeffvalues(cf);
    Rsq(L+1) = rstat4.rsquare;
%     [p]=coeffvalues(cf);
    CI=confint(cf)-[IOPfitinc; IOPfitinc];
    CI=CI(2,:);
    FacTrace(L+1) = IOPfitinc(1);
    MEC(L+1)=CI(1);
    CompTrace(L+1)=IOPfitinc(2);
    MEF(L+1) = CI(2);
    Wherab(L+1) = NewInterval(1); 
%     MEB=CI(2);
end





%% Estimation of IOP derivative difference from P derivative
% estder = derQ(2:end) - der2P * FSys;    
% estder = derQ(2:end);

%% Figures %%
    
    %% overall schematic
    fig51=figure(51)
    hold all
%     axes(ax1)
subplot(2,1,1)
    plot((t(kinterval)-steptmin(1,4))*60,Q(kinterval),'color','g','LineWidth',2);
    hold on
%     axes(ax2)
subplot(2,1,2)

    plot((t(kinterval)-steptmin(1,4))*60,P(kinterval),'color','g','LineWidth',2);



    %% Q and P plots
    if PLOT
fig20=figure(120);
set(fig20,'position',[100 100 1100 580],'color','w')
sgtitle('Q (\color{blue}Filter \color{red}NoFilter)');subplot(3,4,i);hold all;box on;grid on 
% plot (Q,'color','b');
hold on 
% plot(Qnofilt,'color','r');
% plot(actoff(target),Q(actoff(target)),'marker', 'O','color','g');
% plot(kinterval,Q(kinterval),'g','linewidth', 1.2);
plot((t-min(t))*60,Q,'color','b','LineWidth',2);
set(gca,'xlim',[0.7*min(kinterval) 1.2*max(kinterval)]) 
xlabel('datapoints'); ylabel('Q (nl/min)');
title(['Step ' num2str(target)]);
end
    if PLOT
fig21=figure(121);
set(fig21,'position',[100 100 1100 580],'color','w')
% set(fig21,'position',[100 100 450 200],'color','w')
sgtitle('P (\color{blue}Filter \color{red}NoFilter)');subplot(3,4,i);hold all;box on;grid on 
% plot (P,'color','b');
plot ((t-min(t))*60,P,'color','b','LineWidth',2);
hold on 
% plot(Pnofilt,'color','r');
% plot(actoff(target),P(actoff(target)),'marker', 'O','color','g');
plot((t(kinterval)-min(t))*60,P(kinterval),'g','linewidth', 2.5);
% set(gca,'xlim',[0.7*min(kinterval) 1.2*max(kinterval)]) 
set(gca,'xlim',[-10 200])
% xlabel('datapoints'); ylabel('P (mmHg)');

ylabel('$Pressure$ (mmHg)','Interpreter','latex','FontSize',12);
xlabel('$Time$ (Sec)','Interpreter','latex','FontSize',12);
title(['Step ' num2str(target)]);
end



    %% k1 and k2 
% if PLOT
% fig1 = figure(101);
% set(fig1,'position',[100 100 1100 580],'color','w')
% sgtitle('k2 vs k1'); 
% subplot(3,4,i);hold all;box on;grid on
% scatter(k1(1:actoff(target)),k2(1:actoff(target)),8,'b','filled'); 
% hold on
% scatter(k1(actoff(target):actoff(target)+300),k2(actoff(target):actoff(target)+300),8,'b','filled'); 
% scatter(k1(kinterval),k2(kinterval),8,'g','filled');
% plot(k1(actoff(target)),k2(actoff(target)),'marker','O','color','g');  % plot actuator stop point
% text(0.05,0.9, rval,'Units','normalized','fontsize',8,'HorizontalAlignment','left','color','k')
% plot(k1val,k2val,'color','k');   % plot fitter line
% kfitdisp = sprintf('y = %.1f x + (%.1f)', kfit(1), kfit(2));
% text(0.05,0.80, kfitdisp,'Units','normalized','fontsize',8,'HorizontalAlignment','left','color','k')
% 
% set(gca,'xlim',[k1val(1) k1val(2)]) 
% ylabel('k2 (nl/mmHg)'); xlabel('k1 (min)');
% title(['Step ' num2str(target)]);
% 
% end
% 
% if PLOT
% fig5=figure(105);
% set(fig5,'position',[100 100 1100 580],'color','w')
% sgtitle('k1 (\color{blue}Filter \color{red}NoFilter)');subplot(3,4,i);hold all;box on;grid on 
% % plot ((1:actoff(target)+gracep),k1(1:actoff(target)+gracep),'color','b');
% plot (k1,'color','b');
% hold on 
% plot(k1nofilt,'color','r');
% plot(actoff(target),k1(actoff(target)),'marker', 'O','color','g');
% plot(kinterval,k1(kinterval),'g','linewidth', 1.2);
% set(gca,'xlim',[0.7*min(kinterval) 1.2*max(kinterval)]) 
% xlabel('datapoints'); ylabel('k1');
% title(['Step ' num2str(target)]);
% end

% if PLOT
% fig6=figure(106);
% set(fig6,'position',[100 100 1100 580],'color','w')
% sgtitle('k2 (\color{blue}Filter \color{red}NoFilter)');
% subplot(3,4,i);hold all;box on;grid on 
% plot (k2,'color','b');
% hold on 
% plot(k2nofilt,'color','r');
% plot(actoff(target),k2(actoff(target)),'marker', 'O','color','g');
% plot(kinterval,k2(kinterval),'g','linewidth', 1.2);
% set(gca,'xlim',[0.7*min(kinterval) 1.2*max(kinterval)]) 
% xlabel('datapoints'); ylabel('k2');
% title(['Step ' num2str(target)]);
% end

%% IOP k2N vs k1N plot
if PLOT
%     close all
    
fig2 = figure(102);
% set(fig2,'position',[100 100 450 400],'color','w')
set(fig2,'position',[100 100 1100 580],'color','w')
sgtitle('k2N vs k1N'); 
subp =  subplot(3,4,i);

hold all;box on;grid off 
marker_size_scatter = 20;
scatter(k1N(kinterval(1):kinterval(end)+100),k2N(kinterval(1):kinterval(end)+100),...
    marker_size_scatter,...
    (t(kinterval(1):kinterval(end)+100)-min(t(kinterval(1):kinterval(end)+100)))*60,'filled')
% scatter(k1N(actoff(target):actoff(target)+400),k2N(actoff(target):actoff(target)+400),...
%     marker_size_scatter,...
%     (t(actoff(target):actoff(target)+400)-min(t(actoff(target):actoff(target)+400)))*60,'filled')
% scatter(k1N(actoff(target):actoff(target)+500),k2N(actoff(target):actoff(target)+500),10,'b','filled');
hold on
% scatter(k1N(kinterval),k2N(kinterval),8,'g','filled');
% plot(k1N(actoff(target)),k2N(actoff(target)),'marker','O','color','r');
%  rval3
text(0.05,0.9, sprintf('R^2 = %.4f',rstat3.rsquare),'Units','normalized',...
    'FontName','Times','FontSize',9,...
    'HorizontalAlignment','left','color','k')

%%overlay the fit; the entire line region as dashed, and kinterval is the
%%region used for curve fitting
plot(k1Nval,k2Nval,'-','color',[0,0,0,.8],'LineWidth',2);
% % plot([k1Nval(1),max(k1N(kinterval))],[k2Nval(1),max(k2N(kinterval))],'-','color',[0,0,0,1],'LineWidth',2)
% plot(k1N(kinterval),k2N(kinterval),'.k','MarkerSize',10)

% kfitdisp = ['y = ' num2str(kfit(1),1) '* x + ' num2str(kfit(2),1)];
kfitdisp = sprintf('y = %.1f x + %.1f', IOPfit(1), IOPfit(2));
text(0.05,0.80, kfitdisp,'Units','normalized',...
    'FontName','Times','FontSize',9,...
    'HorizontalAlignment','left','color','k','Interpreter','latex')

ylabel('$k_2$ (nl/mmHg)','Interpreter','latex','FontSize',12);
xlabel('$k_1$ (min)','Interpreter','latex','FontSize',12);
xlim([k1Nval(1)-10 k1Nval(2)+30]);

colormap coolwarm
title(['Step ' num2str(target)]);
% h = colorbar;
% ylabel(h,'Time after pressure change (Sec)','Interpreter','latex','FontSize',12)

% set(subp,'ylim',[-300 7000],'xlim',[-10 210],'FontName','Times','FontSize',12,'interpreter','latex');
% set(subp,'FontName','Times','FontSize',12,'interpreter','latex');



end

%%

if PLOT
fig131 = figure(131);
% set(fig131,'position',[100 100 1100 580],'color','w')
set(fig21,'position',[100 100 450 200],'color','w')
% sgtitle('Increment Outflow'); 
% subplot(3,4,i);
hold all;box on;grid off 
% errorbar(Wherab,FacTrace,MEC,'o','markeredgecolor','none',... % going up
%                     'color','r','markerfacecolor','b','linestyle','none','markersize',8); 
erh = errorbar((t(Wherab)-min(t))*60,FacTrace,MEC,'o','markeredgecolor','none',... % going up
                    'color','r','markerfacecolor','b','linestyle','none','markersize',8); 
erh.LineWidth = 1;                

yline(Cj(target),'-.','linewidth',2,'color','k');


% RsqStr = repelem({'R-Sq='} , numel(Rsq));                               
% formatSpec = "%.2f";
% stroutflo = compose(formatSpec,Rsq) ;                           
% text(Wherab,FacTrace+3,stroutflo,'VerticalAlignment','bottom','FontSize',8,'Rotation',90)               
% set(gca,'xlim',[Wherab(1)-20 Wherab(end)+20]);
% ylabel('C_eye (nl/min/mmHg)'); xlabel('Data Point');
incylim=get(gca,'YLim');
set(gca,'ylim',[incylim(1) incylim(2)+1],'FontSize',12)
ylabel('$Instantaneous\ Facility$ (nl/min/mmHg)','Interpreter','latex','FontSize',14);
xlabel('$Time$ (Sec)','Interpreter','latex','FontSize',14);
% title(['Step ' num2str(target)]);
end

if PLOT
fig132 = figure(132);
set(fig132,'position',[100 100 1100 580],'color','w')
% set(fig21,'position',[100 100 450 200],'color','w')

sgtitle('Increment Compliance'); 
subplot(3,4,i);
hold all;box on;grid off 
errorbar((t(Wherab)-min(t))*60,transpose(CompTrace),MEF,'o','markeredgecolor','none',... % going up
                    'color','r','markerfacecolor','b','linestyle','none','markersize',8);  
% [cfcomp,~] = fit((t(Wherab)-min(t))*60,transpose(CompTrace),'poly5'); 
% Xcomp = [0:1:max((t(Wherab)-min(t))*60)];
% Ycomp = feval(cfcomp,Xcomp);
% plot(Xcomp,Ycomp,'color','g');
                
                
% RsqStr = repelem({'R-Sq='} , numel(Rsq));                               
% formatSpec = "%.2f";
% stroutflo = compose(formatSpec,Rsq) ;                           
% text(Wherab,FacTrace+3,stroutflo,'VerticalAlignment','bottom','FontSize',8,'Rotation',90)   
yline(92,'-.','linewidth',2,'color','k');
incylim=get(gca,'YLim');
set(gca,'ylim',[incylim(1) incylim(2)+20])
% set(gca,'xlim',[Wherab(1)-20 Wherab(end)+20]);
% ylabel('\phi_eye (nl/min/mmHg)'); xlabel('Data Point');
ylabel('$Instant\ Compliance$ (nl/mmHg)','Interpreter','latex','FontSize',12);
xlabel('$Time$ (Sec)','Interpreter','latex','FontSize',12);
title(['Step ' num2str(target)]);
end


if PLOT
fig110=figure(110);
set(fig110,'position',[100 100 1100 580],'color','w')
sgtitle('k1N (\color{blue}Filter \color{red}NoFilter)');subplot(3,4,i);hold all;box on;grid on 
% plot ((1:actoff(target)+gracep),k1(1:actoff(target)+gracep),'color','b');
plot (k1N,'color','b');
hold on 
plot(k1Nnofilt,'color','r');
plot(actoff(target),k1N(actoff(target)),'marker', 'O','color','g');
plot(kinterval,k1N(kinterval),'g','linewidth', 1.2);
set(gca,'xlim',[0.7*min(kinterval) 1.2*max(kinterval)]) 
xlabel('datapoints'); ylabel('k1N');
title(['Step ' num2str(target)]);
end

if PLOT
fig111=figure(111);
set(fig111,'position',[100 100 1100 580],'color','w')
sgtitle('k2N (\color{blue}Filter \color{red}NoFilter)');
subplot(3,4,i);hold all;box on;grid on 
plot (k2N,'color','b');
hold on 
plot(k2Nnofilt,'color','r');
plot(actoff(target),k2N(actoff(target)),'marker', 'O','color','g');
plot(kinterval,k2N(kinterval),'g','linewidth', 1.2);
set(gca,'xlim',[0.7*min(kinterval) 1.2*max(kinterval)]) 
xlabel('datapoints'); ylabel('k2N');
title(['Step ' num2str(target)]);
end

%% IOP and derivative of IOP plots
if PLOT
fig112=figure(112);
set(fig112,'position',[100 100 1100 580],'color','w')
sgtitle('IOP (\color{blue}Filter \color{red}NoFilter)');subplot(3,4,i);hold all;box on;grid on 
% plot ((1:actoff(target)+gracep),k1(1:actoff(target)+gracep),'color','b');
plot (IOP,'color','b');
hold on 
plot(IOPnofilt,'color','r');
plot(actoff(target),IOP(actoff(target)),'marker', 'O','color','g');
plot(kinterval,IOP(kinterval),'g','linewidth', 1.2);
set(gca,'xlim',[0.7*min(kinterval) 1.2*max(kinterval)]) 
xlabel('datapoints'); ylabel('IOP (mmHg)');
title(['Step ' num2str(target)]);
end


if PLOT
fig113=figure(113);
set(fig113,'position',[100 100 1100 580],'color','w')
sgtitle('IOP'' (\color{blue}Filter \color{red}NoFilter)');subplot(3,4,i);hold all;box on;grid on 
% plot ((1:actoff(target)+gracep),k1(1:actoff(target)+gracep),'color','b');
plot (derIOP,'color','b');
hold on 
plot(derIOPnofilt,'color','r');
plot(actoff(target),derIOP(actoff(target)),'marker', 'O','color','g');
plot(kinterval,derIOP(kinterval),'g','linewidth', 1.2);
set(gca,'xlim',[0.7*min(kinterval) 1.2*max(kinterval)]) 
xlabel('datapoints'); ylabel('dIOP/dt (mmHg/min)');
title(['Step ' num2str(target)]);
end
%% IOP method validation plots
if PLOT    
fig7=figure(107);
set(fig7,'position',[100 100 1100 580],'color','w')
sgtitle('\color{blue}Q \color{black}and \color{red}(P'' x \phi_S)');
subplot(3,4,i);hold all;box on;grid on 
plot(Q,'b');
hold on 
plot(derP*FSys,'r');
plot(kinterval,Q(kinterval),'color','g');
xlabel('datapoints'); ylabel('nl/min');
title(['Step ' num2str(target)]);
end

if PLOT    
fig8=figure(108);
set(fig8,'position',[100 100 1100 580],'color','w')
sgtitle('\color{blue}Q'' \color{black}and \color{red}(P'''' x \phi_S)');
subplot(3,4,i);hold all;box on;grid on 
plot(derQ,'b');
hold on 
plot(der2P*FSys,'r');
plot(kinterval,derQ(kinterval),'color','g');
xlabel('datapoints'); ylabel('nl/min^2');
title(['Step ' num2str(target)]);
end

%% Instant Compliance and outflow
if PLOT
fig3 = figure(103);
set(fig3,'position',[100 100 1100 580],'color','w')
sgtitle('Compliance (\color{blue}Eq.10 \color{red}IOPMethod)'); 
subplot(3,4,i);hold all;box on;grid on 
plot(InstComp,'color','b');
% hold on
% plot(InstComp,'color','r');
plot(actoff(target),InstComp(actoff(target)),'marker', 'O','color','b');
scatter(kinterval,InstComp(kinterval),8,'g','filled');

% plot([actoff(target)+1:actoff(target)+6],FEye(actoff(target)+1:actoff(target)+6),'color','g','linewidth',1.2)
ylabel('\phi (nl/min)'); xlabel('datapoints');
title(['Step ' num2str(target)]);
xlim([1 350]);
end

if PLOT
fig130 = figure(130);
set(fig130,'position',[100 100 1100 580],'color','w')
sgtitle('Instant Outflow'); 
subplot(3,4,i);
hold all;box on;grid off 
plot(InstOut,'color','b');
% scatter(actoff(target):actoff(target)+500,InstOut(actoff(target):actoff(target)+500),8,'b','filled');
% hold on
scatter(kinterval,InstOut(kinterval),8,'g','filled');
plot(actoff(target),InstOut(actoff(target)),'marker','O','color','r');
set(gca,'xlim',[actoff(target) actoff(target)+500]);
ylabel('C_e_y_e (nl/min/mmHg)'); xlabel('Data Point');
title(['Step ' num2str(target)]);
% colormap coolwarm; colorbar
end

%% derivative of P plot
if PLOT
fig4=figure(104);
set(fig4,'position',[100 100 1100 580],'color','w')
sgtitle('P'' (\color{blue}Filter \color{red}NoFilter)');subplot(3,4,i);hold all;box on;grid on 
% plot ((1:actoff(target)+gracep),k1(1:actoff(target)+gracep),'color','b');
plot (derP,'color','b');
hold on 
plot(derPnofilt,'color','r');
plot(actoff(target),derP(actoff(target)),'marker', 'O','color','g');
plot(kinterval,derP(kinterval),'g','linewidth', 1.2);
set(gca,'xlim',[0.7*min(kinterval) 1.2*max(kinterval)]) 
xlabel('datapoints'); ylabel('dP/dt (mmHg/min)');
title(['Step ' num2str(target)]);
end

%% derivative of Q plot
if PLOT
fig114=figure(114);
set(fig114,'position',[100 100 1100 580],'color','w')
sgtitle('Q'' (\color{blue}Filter \color{red}NoFilter)');subplot(3,4,i);hold all;box on;grid on 
% plot ((1:actoff(target)+gracep),k1(1:actoff(target)+gracep),'color','b');
plot (derQ,'color','b');
hold on 
plot(derQnofilt,'color','r');
plot(actoff(target),derQ(actoff(target)),'marker', 'O','color','g');
plot(kinterval,derQ(kinterval),'g','linewidth', 1.2);
set(gca,'xlim',[0.7*min(kinterval) 1.2*max(kinterval)]) 
xlabel('datapoints'); ylabel('dP/dt (mmHg/min)');
title(['Step ' num2str(target)]);
end
%% Second derivative of P plot
if PLOT
fig4=figure(109);
set(fig4,'position',[100 100 1100 580],'color','w')
sgtitle('der2P (\color{blue}Filter \color{red}NoFilter)');subplot(3,4,i);hold all;box on;grid on 
% plot ((1:actoff(target)+gracep),k1(1:actoff(target)+gracep),'color','b');
plot (der2P,'color','b');
hold on 
plot(der2Pnofilt,'color','r');
plot(actoff(target),der2P(actoff(target)),'marker', 'O','color','g');
plot(kinterval,der2P(kinterval),'g','linewidth', 1.2);
set(gca,'xlim',[0.7*min(kinterval) 1.2*max(kinterval)]) 
xlabel('datapoints'); ylabel('dP/dt2 (mmHg/min)');
title(['Step ' num2str(target)]);
end

%% Q and P plots
% if PLOT
% fig15 = figure(115),
% set(fig15,'position',[100 100 1100 580],'color','w')
% sgtitle('P (\color{blue}Filter \color{red}NoFilter)');subplot(3,4,i);hold all;box on;grid on 
% plot (P,'color','b');
% hold on 
% plot(Pnofilt,'color','r');
% plot(actoff(target),P(actoff(target)),'marker', 'O','color','g');
% plot(kinterval,P(kinterval),'g','linewidth', 1.2);
% set(gca,'xlim',[0.7*min(kinterval) 1.2*max(kinterval)]) 
% xlabel('datapoints'); ylabel('P (mmHg)');
% title(['Step ' num2str(target)]);
% end

% if PLOT
% fig16 = figure(116),
% set(fig16,'position',[100 100 1100 580],'color','w')
% % sgtitle('Q (\color{blue}Filter \color{red}NoFilter)');subplot(3,4,i);hold all;box on;grid on 
% sgtitle('Q (\color{red}flow sensor \color{black}C_q(P_a - P))');subplot(3,4,i);hold all;box on;grid on
% % plot (Q,'color','b');
% hold on 
% plot(t-min(t),Q,'color','r');
% % plot(Qref,'color','k')
% plot(t-min(t),Qref,'color','k')
% % plot(actoff(target),Q(actoff(target)),'marker', 'O','color','g');
% % plot(kinterval,Q(kinterval),'g','linewidth', 1.2);
% % set(gca,'xlim',[0.7*min(kinterval) 1.2*max(kinterval)]) 
% % xlabel('datapoints'); ylabel('Q (nl/min)');
% xlabel('time (min)'); ylabel('Q (nl/min)');
% title(['Step ' num2str(target)]);
% end




% if PLOT
% fig10=figure(110);
% set(fig10,'position',[100 100 1100 580],'color','w')
% sgtitle('k1 time derivative');
% subplot(3,4,i);hold all;box on;grid on 
% plot ((2:actoff(target)+gracep+1),derk1(1:actoff(target)+gracep),'color','b');
% hold on 
% plot((2:actoff(target)+gracep+1),dk1nofilt(1:actoff(target)+gracep),'color','r');
% plot((actoff(target)+1),dk1nofilt(actoff(target)),'marker', 'O','color','g');
% xlabel('datapoints'); ylabel('dk1/dt ');
% title(['Step ' num2str(target)]);
% end

% if PLOT
% fig11=figure(111);
% set(fig11,'position',[100 100 1100 580],'color','w')
% sgtitle('k2 time derivative');
% subplot(3,4,i);hold all;box on;grid on 
% plot ((2:actoff(target)+gracep+1),derk2(1:actoff(target)+gracep),'color','b');
% hold on 
% % plot((2:actoff(target)+gracep+1)*60,dk2nofilt(1:actoff(target)+gracep)/dt,'color','r');
% plot((2:actoff(target)+gracep+1),dk2nofilt(1:actoff(target)+gracep),'color','r');
% plot((actoff(target)+1),dk2nofilt(actoff(target)),'marker', 'O','color','g');
% xlabel('datapoints'); ylabel('dk2/dt (nl/mmHg/min)');
% title(['Step ' num2str(target)]);
% end

if ~PLOT
    fig1 = [];
    fig2 = [];
    fig3 = [];
    fig4 = [];
    fig5 = [];
    fig6 = [];
    fig7 = [];
    fig8 = [];
%     fig9 = [];
%     fig10 = [];
%     fig11 = [];                
%     fig12 = [];
%     fig13 = [];
    fig15 = [];
    fig16 = [];
end


end


