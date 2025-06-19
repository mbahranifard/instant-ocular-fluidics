function test_IOP()
% close all
clc

%time vector
t = linspace(0,5,400);

%pressure
p = 20 - 10*(1-exp(-t/1));

%syringe system compliance
phi_s = 0;

%diagnostic figure
figure('Units','pixels','Position',[200,200,1000,600])

Q0 = linspace(2000,8000,6);

for i=1:length(Q0)
    %flow
    % Q = 6000-2000*1./(t/5+.1);
    Q = -(2000*1./(t/5+.1) - Q0(i));
    
    subplot(2,2,1)
    hold on
    plot(t,p,'-')
    ylabel('Pressure')
    xlabel('Time')

    yyaxis right
    hold on
    grid on

    plot(t,Q,'-')
    ylabel('Flow')


    % fit variables
    k1 = p./(gradient(p)./gradient(t));
    k2 = Q./(gradient(p)./gradient(t))- phi_s;


    subplot(2,2,2)
    hold on
    plot(k1,k2,'-.','LineWidth',2)
    hold on
    ylabel('k2')
    xlabel('k1')
    leg{i} = sprintf('Q_0 = %.2f',Q0(i));
    
    
    subplot(2,2,3)
%     hold on
    semilogy(t,k1,'-.','LineWidth',2)
    hold on
    ylabel('k1')
    xlabel('t')
    
    subplot(2,2,4)
%     hold on
    semilogy(t,k2,'-.','LineWidth',2)
    hold on
    ylabel('k2')
    xlabel('t')
end
subplot(2,2,2)
hold on
legend(leg,'location','best')
end