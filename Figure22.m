
%% DR potential comparison and customer targeting
% @Copyright by NingQi-Tsinghua University-10/04/2020
% contact-qn18@mails.tsinghua.edu.cn

load DRpotentialresult_1.mat

for i=1:size(Cfinal2,1)
    
    C = Cfinal2(i,1);
    R = RTfinal(i,1);

    Tout= 40;
    Tran = 0.5;%tduration
    Tin1 = Tinmean(i,1);
    dt=2;
    Tin2 = Tin1 + dt;
    M = exp(Tran/(R*C));
    Preduce1(i,1) = ((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) <= 0).*((Tout-Tin1)./R)...
        +((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) > 0).*((M.*dt)./(R.*(M-1)))*onoffprobability(i,13);
    
end

[P1,I1]=sort(Preduce1,'descend');

for i=1:size(Cfinal2,1)
    
    C = Cfinal2(i,1);
    R = RTfinal(i,1);

    Tout= 35;
    Tran = 0.5;%tduration
    Tin1 = Tinmean(i,1);
    dt=2;
    Tin2 = Tin1 + dt;
    M = exp(Tran/(R*C));
    Preduce2(i,1) = ((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) <= 0).*((Tout-Tin1)./R)...
        +((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) > 0).*((M.*dt)./(R.*(M-1)))*onoffprobability(i,18);
    
end

[P2,I2]=sort(Preduce2,'descend');

for i=1:size(Cfinal2,1)
    
    C = Cfinal2(i,1);
    R = RTfinal(i,1);

    Tout= 40;
    Tran = 1.5;%tduration
    Tin1 = Tinmean(i,1);
    dt=2;
    Tin2 = Tin1 + dt;
    M = exp(Tran/(R*C));
    Preduce3(i,1) = ((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) <= 0).*((Tout-Tin1)./R)...
        +((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) > 0).*((M.*dt)./(R.*(M-1)))*onoffprobability(i,13);
    
end

[P3,I3]=sort(Preduce3,'descend');

for i=1:size(Cfinal2,1)
    
    C = Cfinal2(i,1);
    R = RTfinal(i,1);

    Tout= 35;
    Tran = 1.5;%tduration
    Tin1 = Tinmean(i,1);
    dt=2;
    Tin2 = Tin1 + dt;
    M = exp(Tran/(R*C));
    Preduce4(i,1) = ((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) <= 0).*((Tout-Tin1)./R)...
        +((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) > 0).*((M.*dt)./(R.*(M-1)))*onoffprobability(i,18);
    
end

[P4,I4]=sort(Preduce4,'descend');


%% targeting VS random

for i=1:119
    
a1=randperm(i);
b1=Preduce1(a1');
P1random(i,1)=sum(b1);
P1strategy(i,1)=sum(P1(1:i));

end

for i=1:119
    
a2=randperm(i);
b2=Preduce2(a2');
P2random(i,1)=sum(b2);
P2strategy(i,1)=sum(P2(1:i));

end

for i=1:119
    
a3=randperm(i);
b3=Preduce3(a3');
P3random(i,1)=sum(b3);
P3strategy(i,1)=sum(P3(1:i));

end

for i=1:119
    
a4=randperm(i);
b4=Preduce4(a4');
P4random(i,1)=sum(b4);
P4strategy(i,1)=sum(P4(1:i));

end

figure(1)
set(gcf,'unit','centimeters','position',[0,0,8,6])
subplot(2,2,1)
plot(P1strategy,'b-','LineWidth',1)
hold on
plot(P1random,'r-.','LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',8)

hold on
set(gca, 'Position', [ 0.1, 0.6, 0.35, 0.35 ]);
xlim([0 120])
ylim([0 120])

subplot(2,2,2)
plot(P2strategy,'b-','LineWidth',1)
hold on
plot(P2random,'r-.','LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',8)

set(gca, 'Position', [ 0.6, 0.6, 0.35, 0.35 ]);
xlim([0 120])
ylim([0 120])
subplot(2,2,3)
plot(P3strategy,'b-','LineWidth',1)
hold on
plot(P3random,'r-.','LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',8)

xlim([0 120])
ylim([0 50])
set(gca, 'Position', [ 0.1, 0.1, 0.35, 0.35 ]);

subplot(2,2,4)
plot(P4strategy,'b-','LineWidth',1)
hold on
plot(P4random,'r-.','LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',8)

set(gca, 'Position', [ 0.6, 0.1, 0.35, 0.35 ]);
xlim([0 120])
ylim([0 50])