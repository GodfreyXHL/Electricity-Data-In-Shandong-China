% The robustnessevaluation.m is used for robustness evaluation of DR  potential
% @Copyright by NingQi-Tsinghua University-10/04/2020
% contact-qn18@mails.tsinghua.edu.cn

% load data from parameter estimation
clear all;
load DRpotentialresult_1.mat
C = 558.65 %Cfinal2(10,1);
R = 0.0012 %RTfinal(10,1);


tout= 30:0.1:40;% Tout
tran = 0:0.1:2;%tduration
[Tran,Tout] = meshgrid(tran,tout);

%Tin1 = Tinmean(10,1)*ones(length(Tout(:,1)),length(Tout(1,:))); %Tinmean(10,1) 室内温度
Tin1 = 20.6318*ones(length(Tout(:,1)),length(Tout(1,:))); %Tinmean(10,1) 室内温度

dt=2; %dT
Tin2 = Tin1 + dt;
M = exp(Tran/(R*C));
Preduce = ((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) <= 0).*((Tout-Tin1)./R)...
    +((Tran-(R*C*log((Tout-Tin1)./(Tout-Tin2)))) > 0).*((M.*dt)./(R.*(M-1)));

figure(1)
set(gcf,'unit','centimeters','position',[0,0,8,6])
meshz(Tran,Tout,Preduce);
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}DR Duration (h)')
ylabel('\fontsize{8}\fontname{Times new roman}Ambient Temperature')
hold on
plot3(0.1*ones(1,size(tout,2)),tout,Preduce(:,2)','r-','LineWidth',4);
zlabel('\fontsize{8}\fontname{Times new roman}DR Potential (kW)')
view([90 90 90])
zlim([0 15000])

Trange=20:0.05:30;
E=6.112*exp(17.67*Trange./(Trange+243.5));
RH1=0.5;
E1=0.5*E(1,1);
anew=0.272;
bnew=0.245;
cnew=-7.2;
P=anew*Trange+bnew*E1+cnew-2.4;
Pmv=abs(P);

for i=1:size(Pmv,2)
    if Pmv(1,i)<0.2
        Pmvnew(1,i)=0;
    end
    if Pmv(1,i)>=0.2
        Pmvnew(1,i)=Pmv(1,i)-0.2;
    end
end

figure
plot(Pmv);
