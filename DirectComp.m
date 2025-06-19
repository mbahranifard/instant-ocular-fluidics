function[FEyej,Pf,fig1,fig2,fig4,fig6,fig7,fig8,slope] = DirectComp (stepP,stepQ,stepPa,steptmin,actoff,Cc,FSys,Yr,B,PrC,KeepStepsUp,gracep,PLOT)


%%
if PLOT
fig1 = figure(101);
set(fig1,'position',[100 100 2000 1200],'color','w')
sgtitle('Pressure');
subplot(3,4,i);hold all;box on;grid on 
plot (P,'color','b');
hold on 
plot(Pnofilt,'color','r');
plot(actoff(target),Pnofilt(actoff(target)),'marker', 'O','color','g');
ylabel ('P (mmHg)'); xlabel('DataPoint');
title(['Step ' num2str(target)]);
end

%% first derivative
pest = mean(diff (P(1:actoff(target))));
% dt = diff(t(1:2));
Pest = [P(1)-pest;P];
dP = Pest(2:end)-Pest(1:end-1);

% dP = [dP(1)-pest;dP];
% d2P =  dP(2:end)-dP(1:end-1);
% figure(102);
% plot(t(1:actoff(target)+500),dP(1:actoff(target)+500))
% hold on 
% plot(t(actoff(target)),dP(actoff(target)),'marker', 'O','color','r');

dPnofilt = dP;
iSG = round(0.01/dt);
dP = sgolayfilt(dP,1,iSG*2+1,[],1);

if PLOT
fig2=figure(102);
set(fig2,'position',[100 100 2000 1200],'color','w')
sgtitle('Pressure time derivative');
subplot(3,4,i);hold all;box on;grid on 
plot (dP(1:actoff(target)+gracep)/dt,'color','b');
hold on 
plot(dPnofilt(1:actoff(target)+gracep)/dt,'color','r');
plot(actoff(target),dPnofilt(actoff(target))/dt,'marker', 'O','color','g');
xlabel('datapoint'); ylabel('dP/dt (mmHg/min)');
title(['Step ' num2str(target)]);
end
%% second derivative
% pest2  = mean(diff (dP(1:actoff(target))));
% dPest = [dP(1)-pest;dP];
% d2P =  dPest(2:end)-dPest(1:end-1);
% % figure(104);
% % plot(t(1:actoff(target)+500),d2P(1:actoff(target)+500))
% % hold on 
% % plot(t(actoff(target)),d2P(actoff(target)),'marker', 'O','color','r');
% 
% d2Pnofilt = d2P;
% iSG = round(0.01/dt);
% d2P = sgolayfilt(d2P,1,iSG*2+1,[],1);
% 
% fig3 = figure(103);
% set(fig3,'position',[100 100 2000 1200],'color','w')
% title('Pressure second derivative');
% subplot(3,4,i);hold all;box on;grid on 
% 
% plot (d2P(1:actoff(target)+gracep),'color','b');
% hold on 
% plot(d2Pnofilt(1:actoff(target)+gracep),'color','r');
% plot(actoff(target),d2Pnofilt(actoff(target)),'marker', 'O','color','g');
% title(['Step ' num2str(target)]);
% xlabel('datapoint'); ylabel('P'' (mmHg/s^2)');
%%

%% IOP

% Panofilt = Pa;
% iSG = round(0.01/dt);
% Pa = sgolayfilt(Pa,1,iSG*2+1,[],1);


Qnofilt = Q;
iSG = round(0.01/dt);
Q = sgolayfilt(Q,1,iSG*2+1,[],1);

IOP = P - (Q/Cc) + FSys/Cc * dP/dt;    %%%% if noisy use Cq(Pa-P)
IOPnofilt = IOP;
iSG = round(0.01/dt);
IOP = sgolayfilt(IOP,1,iSG*2+1,[],1);

if PLOT
fig4 = figure(104);
set(fig4,'position',[100 100 2000 1200],'color','w')
sgtitle ('IOP vs filter');
subplot(3,4,i);hold all;box on;grid on 

plot (IOP(1:actoff(target)+gracep),'color','b');
hold on 
plot(IOPnofilt(1:actoff(target)+gracep),'color','r');
plot(actoff(target),IOPnofilt(actoff(target)),'marker', 'O','color','g');
xlabel('datapoint'); ylabel('IOP (mmHg)');
title(['Step ' num2str(target)]);
end
% fig5 = figure(105);
% set(fig5,'position',[100 100 2000 1200],'color','w')
% title('IOP  vs P');
% subplot(3,4,i);hold all;box on;grid on 
% 
% plot(IOP,'color','b');
% hold on
% plot(P,'color','r');
% ylabel('P (mmHg)'); xlabel('datapoint');
% title(['Step ' num2str(target)]);

%% first derivative of Pa
paest = mean(diff (Pa(1:actoff(target))));
Patest = [Pa(1)-pest;Pa];
dPa = Patest(2:end)-Patest(1:end-1);

if PLOT
fig6 = figure(106);
set(fig6,'position',[100 100 2000 1200],'color','w')
sgtitle('Pa derivative');
subplot(3,4,i);hold all;box on;grid on 
plot(1:actoff(target)+gracep,dPa(1:actoff(target)+gracep)/dt);
hold on 
plot(actoff(target),dPa(actoff(target))/dt,'marker', 'O','color','r');
ylabel('dPa/dt (mmHg/min)'); xlabel('datapoint');
title(['Step ' num2str(target)]);
end
% 
% dIOP = dP/dt-(Cq/Cc)*(dPa/dt - dP/dt)-(FSys/Cc)*(d2P/dt^2);
%% IOP Derivative
IOPest = mean(diff (IOP(1:actoff(target))));
IOPtest = [IOP(1)-pest;IOP];
dIOP = IOPtest(2:end)-IOPtest(1:end-1);

dIOPnofilt = dIOP;
iSG = round(0.01/dt);
dIOP = sgolayfilt(dIOP,1,iSG*2+1,[],1);

if PLOT
fig7 = figure(107);
set(fig7,'position',[100 100 2000 1200],'color','w')
sgtitle('IOP derivative');
subplot(3,4,i);hold all;box on;grid on 
plot (dIOP(1:actoff(target)+gracep)/dt,'color','b');
hold on 
plot(dIOPnofilt(1:actoff(target)+gracep)/dt,'color','r');
plot(actoff(target),dIOPnofilt(actoff(target))/dt,'marker', 'O','color','g');
ylabel('dIOP/dt (mmHg/min)'); xlabel('datapoint');
title(['Step ' num2str(target)]);
end
%%


FEye = (Cc*P - (exp(Yr)*(IOP/PrC).^B + Cc) .*IOP)./(dIOP/dt);
% close all
if PLOT
fig8 = figure(108);
set(fig8,'position',[100 100 2000 1200],'color','w')
sgtitle('Compliance'); 
subplot(3,4,i);hold all;box on;grid on 
plot(FEye);
hold on
plot(actoff(target),FEye(actoff(target)),'marker', 'O','color','r');
plot([actoff(target)+1:actoff(target)+6],FEye(actoff(target)+1:actoff(target)+6),'color','g','linewidth',1.2)
ylabel('\phi (nl/min)'); xlabel('datapoints');
title(['Step ' num2str(target)]);
xlim([130 180]);
end

mFEye = mean(FEye(actoff(target)+1:actoff(target)+6));
Pf (target) = mean(P(actoff(target)+1:actoff(target)+6));
risefit=polyfit(t(actoff(target)+1:actoff(target)+6),FEye(actoff(target)+1:actoff(target)+6),1);
slope(target) = risefit(1);

FEyej (target) = mFEye;

if ~PLOT
    fig1 = [];
    fig2 = [];
    fig3 = [];
    fig4 = [];
    fig5 = [];
    fig6 = [];
    fig7 = [];
    fig8 = [];
end


end

end
