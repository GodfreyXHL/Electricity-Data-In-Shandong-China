
%% DR potential comparison with different time scale
% @Copyright by NingQi-Tsinghua University-10/04/2020
% contact-qn18@mails.tsinghua.edu.cn

load DRpotentialresult_1.mat

%◊Ó Ê   25.5  Ωœ Ê  27.5   ≤ª Ê   30
dtnew=[5;7;9.5];
for i=1:3
    
   C = 558.65 %Cfinal2(10,1);
   R = 0.0012 %RTfinal(10,1);

    Tout= 35;
    Tran = 0.1;%tduration
    Tin1 = 20.6318;
    dt=dtnew(i,1);
    Tin2 = Tin1 + dt;
    M = exp(Tran/(R*C));
    Preduce1(i,1) = ((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) <= 0).*((Tout-Tin1)./R)...
        +((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) > 0).*((M.*dt)./(R.*(M-1)))*onoffprobability(i,18);
    
end

[P1,I1]=sort(Preduce1,'ascend');

for i=1:3
    
    Tran = 0.5;%tduration
   C = 558.65 %Cfinal2(10,1);
   R = 0.0012 %RTfinal(10,1);

    Tout= 35;
    
    Tin1 = 20.6318;
    dt=dtnew(i,1);
    Tin2 = Tin1 + dt;
    M = exp(Tran/(R*C));
    Preduce2(i,1) = ((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) <= 0).*((Tout-Tin1)./R)...
        +((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) > 0).*((M.*dt)./(R.*(M-1)))*onoffprobability(i,18);
    
end

[P2,I2]=sort(Preduce2,'ascend');

for i=1:3
    
 
    Tran = 1;%tduration
   
    
    C = 558.65 %Cfinal2(10,1);
    R = 0.0012 %RTfinal(10,1);

    Tout= 35;
    
    Tin1 = 20.6318;
    
    dt=dtnew(i,1);
    Tin2 = Tin1 + dt;
    M = exp(Tran/(R*C));
    Preduce3(i,1) = ((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) <= 0).*((Tout-Tin1)./R)...
        +((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) > 0).*((M.*dt)./(R.*(M-1)))*onoffprobability(i,18);
    
end

[P3,I3]=sort(Preduce3,'ascend');

for i=1:3
   

    Tout= 35;
    Tran = 1.5;%tduration
    C = 558.65 %Cfinal2(10,1);
    R = 0.0012 %RTfinal(10,1);

    Tout= 35;
    
    Tin1 = 20.6318;
    
    dt=dtnew(i,1);
    Tin2 = Tin1 + dt;
    M = exp(Tran/(R*C));
    Preduce4(i,1) = ((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) <= 0).*((Tout-Tin1)./R)...
        +((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) > 0).*((M.*dt)./(R.*(M-1)))*onoffprobability(i,18);
    
end

[P4,I4]=sort(Preduce4,'ascend');


for i=1:3
    
    
    Tout= 35;
    Tran = 2;%tduration
    C = 558.65 %Cfinal2(10,1);
    R = 0.0012 %RTfinal(10,1);

    Tout= 35;
    
    Tin1 = 20.6318;
    
    dt=dtnew(i,1);
    Tin2 = Tin1 + dt;
    M = exp(Tran/(R*C));
    Preduce5(i,1) = ((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) <= 0).*((Tout-Tin1)./R)...
        +((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) > 0).*((M.*dt)./(R.*(M-1)))*onoffprobability(i,18);
    
end

[P5,I5]=sort(Preduce5,'ascend');

P=[P1(1:3) P2(1:3) P3(1:3) P4(1:3) P5(1:3) ];
I=[I1(1:3) I2(1:3) I3(1:3) I4(1:3) I5(1:3) ];

figure(1)
set(gcf,'unit','centimeters','position',[0,0,8,6])
bar(P');
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman} Duration (h)')
ylabel('\fontsize{8}\fontname{Times new roman}Peak clipping Potential (kW)')
hold on
set(gca,'xticklabel',{'0.1','0.5','1','1.5','2'})
legend('Customer comfort index = 0.3','Customer comfort index = 0.7','Customer comfort index = 1.3','fontsize',8)
ylim([0 17000])

