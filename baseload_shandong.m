%% baseload_garment.m is used for load decomposition of garment
clear all
timescale=24;
%load Data.mat
load data1.mat
%

%% baseload month-April
figure(1)
set(gcf,'unit','centimeters','position',[0,0,16,12])
plot(P4');
hold on
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Time (15min)')
ylabel('\fontsize{10.5}\fontname{Times new roman}Load Power (kW)')
Pcluster=P4;

%% adaptive clustering

for i = 1:5
    [cidx4,cmeans4,sumd4,D4] = kmeans(Pcluster,i,'dist','sqeuclidean');
   
    indicate(i) = sum(sumd4); %SSE
    [DB(i),CH(i),KL(i),Han(i),st(i)] = valid_internal_deviation(Pcluster,cidx4,0);
    Si(i)=mean(silhouette(Pcluster,cidx4));
    %     'Davies-Bouldin (DB)', 'Calinski-Harabasz (CH)', 'Krzanowski-Lai
    %     (KL)','Hartigan',Silhouette
end

figure(2)
set(gcf,'unit','centimeters','position',[0,0,16,12])
subplot(2,3,1)
plot(indicate,'b-','LineWidth',1) 
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Cluster')
ylabel('\fontsize{10.5}\fontname{Times new roman}SSE')
hold on
plot(find(abs(diff(diff(indicate)))==max(abs(diff(diff(indicate)))))+1,indicate(find(abs(diff(diff(indicate)))==max(abs(diff(diff(indicate)))))+1),'r*')
subplot(2,3,2)
plot(DB,'b-','LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Cluster')
ylabel('\fontsize{10.5}\fontname{Times new roman}Davies-Bouldin')
hold on
plot(find(DB==min(DB)),DB(find(DB==min(DB))),'r*')
subplot(2,3,3)
plot(CH,'b-','LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Cluster')
ylabel('\fontsize{10.5}\fontname{Times new roman}Calinski-Harabasz')
hold on
plot(find(CH==max(CH)),CH(find(CH==max(CH))),'r*')
subplot(2,3,4)
plot(KL,'b-','LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Cluster')
ylabel('\fontsize{10.5}\fontname{Times new roman}Krzanowski-Lai')
hold on
plot(find(KL==max(KL)),KL(find(KL==max(KL))),'r*')
subplot(2,3,5)
plot(Han,'b-','LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Cluster')
ylabel('\fontsize{10.5}\fontname{Times new roman}Hartigan')
hold on
plot(find(Han==min(Han)),Han(find(Han==min(Han))),'r*')
subplot(2,3,6)
plot(Si,'b-','LineWidth',1)
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Cluster')
ylabel('\fontsize{10.5}\fontname{Times new roman}Silhouette')
hold on
plot(find(Si==max(Si)),Si(find(Si==max(Si))),'r*')
a4=mode([find(abs(diff(diff(indicate)))==max(abs(diff(diff(indicate)))))+1 find(DB==min(DB)) find(CH==max(CH)) find(KL==max(KL)) find(Han==min(Han)) find(Si==max(Si))]); %% best k

[cidx4,cmeans4,sumd4,D4] = kmeans(Pcluster,2,'dist','sqeuclidean');

[Y4,I4]=sort(mean(cmeans4'));

%% final clustering of April
figure(3)
set(gcf,'unit','centimeters','position',[0,0,8,6])

plot(Pcluster(find(cidx4==I4(1)),:)','b-','LineWidth',1)
times4(1,1)=length(find(cidx4==I4(1)));
hold on

plot(Pcluster(find(cidx4==I4(2)),:)','r-','LineWidth',1)
times4(1,2)=length(find(cidx4==I4(2)));
hold on

set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Time (15min)')
ylabel('\fontsize{8}\fontname{Times new roman}Power (kW)')
% two cluster six days for work

P4baseloadweekday=Pcluster(find(cidx4==I4(2)),:);
P4baseloadweekend=Pcluster(find(cidx4==I4(1)),:);
temp4weekday=T4(find(cidx4==I4(2)),:);
temp4weekend=T4(find(cidx4==I4(1)),:);

% obtain the quasi-baseload for April
figure(4)
set(gcf,'unit','centimeters','position',[0,0,16,12])
plot([temp4weekday;temp4weekend],[P4baseloadweekday;P4baseloadweekend],'o','MarkerSize',2);
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Outdoor Temperature (℉)')
ylabel('\fontsize{10.5}\fontname{Times new roman}Load Power (kW)')

%% clustering for march (k=4,CH)
Pcluster=P3;
[cidx3,cmeans3,sumd3,D3] = kmeans(Pcluster,4,'dist','sqeuclidean');
[Y3,I3]=sort(mean(cmeans3'));

P3baseloadweekday=Pcluster(find(cidx3==I3(3)),:);
P3baseloadweekend=Pcluster(find(cidx3==I3(1)),:);
P3otherweekday=Pcluster(find(cidx3==I3(4)),:);
P3otherweekend=Pcluster(find(cidx3==I3(2)),:);

temp3weekday=T3(find(cidx3==I3(3)),:);
temp3weekend=T3(find(cidx3==I3(1)),:);
temp3otherweekday=T3(find(cidx3==I3(4)),:);
temp3otherweekend=T3(find(cidx3==I3(2)),:);

figure(5)
set(gcf,'unit','centimeters','position',[0,0,8,6])
plot(P3baseloadweekday','k-.','LineWidth',1)
hold on
plot(P3baseloadweekend','b-.','LineWidth',1)
hold on
plot(P3otherweekday','r-.','LineWidth',1)
hold on
plot(P3otherweekend','y-.','LineWidth',1)
hold on
plot(cmeans3(I3(3),:),'k-','LineWidth',4)
hold on
plot(cmeans3(I3(1),:),'b-','LineWidth',4)
hold on
plot(cmeans3(I3(4),:),'r-','LineWidth',4)
hold on
plot(cmeans3(I3(2),:),'y-','LineWidth',4)
hold on
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Time (15min)')
ylabel('\fontsize{8}\fontname{Times new roman}Power (kW)')
% obtain the quasi-baseload 

figure(6)
set(gcf,'unit','centimeters','position',[0,0,16,12])
plot([temp3weekday;temp3weekend],[P3baseloadweekday;P3baseloadweekend],'o','MarkerSize',2);
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Outdoor Temperature (℉)')
ylabel('\fontsize{10.5}\fontname{Times new roman}Load Power (kW)')

%% clustering for May (k=4,CH)

Pcluster=P5;
[cidx5,cmeans5,sumd5,D5] = kmeans(Pcluster,4,'dist','sqeuclidean');
[Y5,I5]=sort(mean(cmeans5'));

P5baseloadweekday=Pcluster(find(cidx5==I5(3)),:);
P5baseloadweekend=Pcluster(find(cidx5==I5(1)),:);
P5otherweekday=Pcluster(find(cidx5==I5(4)),:);
P5otherweekend=Pcluster(find(cidx5==I5(2)),:);
temp5weekday=T5(find(cidx5==I5(3)),:);
temp5weekend=T5(find(cidx5==I5(1)),:);
temp5otherweekday=T5(find(cidx5==I5(4)),:);
temp5otherweekend=T5(find(cidx5==I5(2)),:);

figure(7)
set(gcf,'unit','centimeters','position',[0,0,8,6])
plot(P5baseloadweekday','k-.','LineWidth',1)
hold on
plot(P5baseloadweekend','b-.','LineWidth',1)
hold on
plot(P5otherweekday','r-.','LineWidth',1)
hold on
plot(P5otherweekend','y-.','LineWidth',1)
hold on
plot(cmeans5(I5(3),:),'k-','LineWidth',4)
hold on
plot(cmeans5(I5(1),:),'b-','LineWidth',4)
hold on
plot(cmeans5(I5(4),:),'r-','LineWidth',4)
hold on
plot(cmeans5(I5(2),:),'y-','LineWidth',4)
hold on
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Time (15min)')
ylabel('\fontsize{8}\fontname{Times new roman}Power (kW)')
% obtain the quasi-baseload 

figure(8)
set(gcf,'unit','centimeters','position',[0,0,16,12])
plot([temp5weekday;temp5weekend],[P5baseloadweekday;P5baseloadweekend],'o','MarkerSize',2);
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Outdoor Temperature (℉)')
ylabel('\fontsize{10.5}\fontname{Times new roman}Load Power (kW)')

%% clustering for September (k=4,CH)

Pcluster=P9;
[cidx9,cmeans9,sumd9,D9] = kmeans(Pcluster,4,'dist','sqeuclidean');
[Y9,I9]=sort(mean(cmeans9'));

P9baseloadweekday=Pcluster(find(cidx9==I9(3)),:);
P9baseloadweekend=Pcluster(find(cidx9==I9(1)),:);
P9otherweekday=Pcluster(find(cidx9==I9(4)),:);
P9otherweekend=Pcluster(find(cidx9==I9(2)),:);

temp9weekday=T9(find(cidx9==I9(3)),:);
temp9weekend=T9(find(cidx9==I9(1)),:);
temp9otherweekday=T9(find(cidx9==I9(4)),:);
temp9otherweekend=T9(find(cidx9==I9(2)),:);

figure(9)
set(gcf,'unit','centimeters','position',[0,0,8,6])
plot(P9baseloadweekday','k-.','LineWidth',1)
hold on
plot(P9baseloadweekend','b-.','LineWidth',1)
hold on
plot(P9otherweekday','r-.','LineWidth',1)
hold on
plot(P9otherweekend','y-.','LineWidth',1)
hold on
plot(cmeans9(I9(3),:),'k-','LineWidth',4)
hold on
plot(cmeans9(I9(1),:),'b-','LineWidth',4)
hold on
plot(cmeans9(I9(4),:),'r-','LineWidth',4)
hold on
plot(cmeans9(I9(2),:),'y-','LineWidth',4)
hold on
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Time (15min)')
ylabel('\fontsize{8}\fontname{Times new roman}Power (kW)')
% obtain the quasi-baseload 

figure(10)
set(gcf,'unit','centimeters','position',[0,0,16,12])
plot([temp9weekday;temp9weekend],[P9baseloadweekday;P9baseloadweekend],'o','MarkerSize',2);
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Outdoor Temperature (℉)')
ylabel('\fontsize{10.5}\fontname{Times new roman}Load Power (kW)')

%% clustering for October (k=4,CH)

Pcluster=P10;
[cidx10,cmeans10,sumd10,D10] = kmeans(Pcluster,4,'dist','sqeuclidean');
[Y10,I10]=sort(mean(cmeans10'));

P10baseloadweekday=Pcluster(find(cidx10==I10(3)),:);
P10baseloadweekend=Pcluster(find(cidx10==I10(1)),:);
P10otherweekday=Pcluster(find(cidx10==I10(4)),:);
P10otherweekend=Pcluster(find(cidx10==I10(2)),:);

temp10weekday=T10(find(cidx10==I10(3)),:);
temp10weekend=T10(find(cidx10==I10(1)),:);
temp10otherweekday=T10(find(cidx10==I10(4)),:);
temp10otherweekend=T10(find(cidx10==I10(2)),:);

figure(11)
set(gcf,'unit','centimeters','position',[0,0,8,6])
plot(P10baseloadweekday','k-.','LineWidth',1)
hold on
plot(P10baseloadweekend','b-.','LineWidth',1)
hold on
plot(P10otherweekday','r-.','LineWidth',1)
hold on
plot(P10otherweekend','y-.','LineWidth',1)
hold on
plot(cmeans10(I10(3),:),'k-','LineWidth',4)
hold on
plot(cmeans10(I10(1),:),'b-','LineWidth',4)
hold on
plot(cmeans10(I10(4),:),'r-','LineWidth',4)
hold on
plot(cmeans10(I10(2),:),'y-','LineWidth',4)
hold on
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Time (15min)')
ylabel('\fontsize{8}\fontname{Times new roman}Power (kW)')
% obtain the quasi-baseload 

figure(12)
set(gcf,'unit','centimeters','position',[0,0,16,12])
plot([temp10weekday;temp10weekend],[P10baseloadweekday;P10baseloadweekend],'o','MarkerSize',2);
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Outdoor Temperature (℉)')
ylabel('\fontsize{10.5}\fontname{Times new roman}Load Power (kW)')

%% clustering for November (k=4,CH)

Pcluster=P11;
[cidx11,cmeans11,sumd11,D11] = kmeans(Pcluster,4,'dist','sqeuclidean');
[Y11,I11]=sort(mean(cmeans11'));

P11baseloadweekday=Pcluster(find(cidx11==I11(3)),:);
P11baseloadweekend=Pcluster(find(cidx11==I11(1)),:);
P11otherweekday=Pcluster(find(cidx11==I11(4)),:);
P11otherweekend=Pcluster(find(cidx11==I11(2)),:);

temp11weekday=T11(find(cidx11==I11(3)),:);
temp11weekend=T11(find(cidx11==I11(1)),:);
temp11otherweekday=T11(find(cidx11==I11(4)),:);
temp11otherweekend=T11(find(cidx11==I11(2)),:);

figure(13)
set(gcf,'unit','centimeters','position',[0,0,8,6])
plot(P11baseloadweekday','k-.','LineWidth',1)
hold on
plot(P11baseloadweekend','b-.','LineWidth',1)
hold on
plot(P11otherweekday','r-.','LineWidth',1)
hold on
plot(P11otherweekend','y-.','LineWidth',1)
hold on
plot(cmeans11(I11(3),:),'k-','LineWidth',4)
hold on
plot(cmeans11(I11(1),:),'b-','LineWidth',4)
hold on
plot(cmeans11(I11(4),:),'r-','LineWidth',4)
hold on
plot(cmeans11(I11(2),:),'y-','LineWidth',4)
hold on
set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Time (15min)')
ylabel('\fontsize{8}\fontname{Times new roman}Power (kW)')
% obtain the quasi-baseload 

figure(14)
set(gcf,'unit','centimeters','position',[0,0,8,6])
plot([temp11weekday;temp11weekend],[P11baseloadweekday;P11baseloadweekend],'o','MarkerSize',2);
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Outdoor Temperature (℉)')
ylabel('\fontsize{10.5}\fontname{Times new roman}Load Power (kW)')

%% Statistics baseload set
Pbaseloadweekday=[P3baseloadweekday;P4baseloadweekday;P5baseloadweekday;P9baseloadweekday;P10baseloadweekday;P11baseloadweekday];
Pbaseloadweekend=[P3baseloadweekend;P4baseloadweekend;P5baseloadweekend;P9baseloadweekend;P10baseloadweekend;P11baseloadweekend];
Potherweekday=[P3otherweekday;P5otherweekday;P9otherweekday;P10otherweekday;P11otherweekday];
Potherweekend=[P3otherweekend;P5otherweekend;P9otherweekend;P10otherweekend;P11otherweekend];

figure(15)
set(gcf,'unit','centimeters','position',[0,0,8,6])


plot(Potherweekday','r:','LineWidth',1)
hold on
plot(Potherweekend','y:','LineWidth',1)
hold on
plot(mean(Potherweekday),'r-','LineWidth',2)
hold on
plot(mean(Potherweekend),'y-','LineWidth',2)

plot(Pbaseloadweekday','k:','LineWidth',1)
hold on
plot(Pbaseloadweekend','b:','LineWidth',1)
hold on
plot(mean(Pbaseloadweekday),'k-','LineWidth',2)
hold on
plot(mean(Pbaseloadweekend),'b-','LineWidth',2)

set(gca,'FontName','Times New Roman','FontSize',8)
xlabel('\fontsize{8}\fontname{Times new roman}Time (15min)')
ylabel('\fontsize{8}\fontname{Times new roman}Load Power (kW)')


%% distribution test
%lognormal distribution
k = 2;%kσ

for i = 1:size(Pbaseloadweekday,2)
    parmhat = lognfit(Pbaseloadweekday(:,i),0.05);
    mu = parmhat(1);
    sigma = parmhat(2);
    lowerweekday(i) = exp(mu-k*sigma);
    higherweekday(i) = exp(mu+k*sigma);
    meanweekday(i)=exp(mu);
end

k=2;
for i = 1:size(Pbaseloadweekend,2)
    parmhat = lognfit(Pbaseloadweekend(:,i),0.05);
    mu = parmhat(1);
    sigma = parmhat(2);
    lowerweekend(i) = exp(mu-k*sigma);
    higherweekend(i) = exp(mu+k*sigma);
    meanweekend(i)=exp(mu);
end


lognmeanweekend=meanweekend;
lognmeanweekday=meanweekday;


%% Generalized extremum distribution

k = 2;%kσ
     
for i = 1:size(Pbaseloadweekday,2)
    phat = gevfit(Pbaseloadweekday(:,i),0.05);
 
    mu =phat(3);
    sigma = phat(2);
    lowerweekday(i) = mu-k*sigma;
    higherweekday(i) = mu+k*sigma;
    meanweekday(i)=mu;
end

k=2;
for i = 1:size(Pbaseloadweekend,2)
    phat = gevfit(Pbaseloadweekend(:,i),0.05);
     mu =phat(3);
    sigma = phat(2);
    lowerweekend(i) = mu-k*sigma;
    higherweekend(i) = mu+k*sigma;
    meanweekend(i)=mu;
end

pmeanweekend=meanweekend;
pmeanweekday=meanweekday;

%% weekday|weekend

weekend8_judge=[5;12;19;26];
weekend7_judge=[1;8;15;22;29];
weekend6_judge=[3;10;16;17;18;24];

%% ACLs calculation

air8_es=P8;
for i=1:size(P8,1)
    
air8_es(i,:)=P8(i,:)-lognmeanweekday;
end

for i=1:size(weekend8_judge,1)
    air8_es(weekend8_judge(i),:)=P8(weekend8_judge(i),:)-lognmeanweekend;
end

air7_es=P7;
for i=1:size(P7,1)
    
air7_es(i,:)=P7(i,:)-lognmeanweekday*(100)./100;
end

for i=1:size(weekend7_judge,1)
    air7_es(weekend7_judge(i),:)=P7(weekend7_judge(i),:)-lognmeanweekend;
end

air6_es=P6;
for i=1:size(P6,1)
    
air6_es(i,:)=P6(i,:)-lognmeanweekday*(100)./100;
end

for i=1:size(weekend6_judge,1)
    air6_es(weekend6_judge(i),:)=P6(weekend6_judge(i),:)-lognmeanweekend;
end

air6_es=max(air6_es,0);
air7_es=max(air7_es,0);
air8_es=max(air8_es,0);

airsummer_es=[air6_es;air7_es;air8_es];
tempsummer=[T6;T7;T8];

figure(16)
set(gcf,'unit','centimeters','position',[0,0,16,12])
plot(tempsummer,airsummer_es,'bo','MarkerSize',2);
hold on
set(gca,'FontName','Times New Roman','FontSize',10.5)
xlabel('\fontsize{10.5}\fontname{Times new roman}Outdoor Temperature (℉)')
ylabel('\fontsize{10.5}\fontname{Times new roman}ACL (kW)')

save baseloadgarment.mat
airsummer1=zeros(92*24,1);
tempsummer1=zeros(92*24,1);
t=1;
for i=1:92
    for j=1:24
        airsummer1(t,1)=airsummer_es(i,j);
        tempsummer1(t,1)=tempsummer(i,j);
        t=t+1;
    end
end
air7new=zeros(744,1)
t=1;
for i=1:31
    for j=1:24
        air7new(t,1)=air7_es(i,j);
        t=t+1;
    end
end
    
    
figure(17)
set(gcf,'Units','centimeters','Position',[6 6 14.5 9]);%设置图片大小
subplot(2,1,1)
a=plot(airsummer1);
%axis([0,24,30000,90000]);
title({'Load Data:';'shoulder season.'},'fontname','Times New Roman'); % 两行的内容用分号隔开，再用大括号括起来。
xlabel('\fontsize{8}\fontname{Times new roman}Time (h)')
ylabel('\fontsize{8}\fontname{Times new roman}Power (kW)')

hold on 
subplot(2,1,2)
b=plot(tempsummer1);
%axis([0,24,30000,90000]);
title({'Load Data:';'summer season.'},'fontname','Times New Roman');
xlabel('\fontsize{8}\fontname{Times new roman}Time (h)')
ylabel('\fontsize{8}\fontname{Times new roman}Power (kW)')




figure(18)
set(gcf,'Units','centimeters','Position',[6 6 16 9]);%设置图片大小
subplot(2,1,1)
a=plot(air7new(400:502,1));
%axis([0,24,30000,90000]);
%title({'Load Data:';'shoulder season.'},'fontname','Times New Roman'); % 两行的内容用分号隔开，再用大括号括起来。
set(gca,'FontName','Times New Roman','FontSize',8)
title('Air conditioning Load','fontname','Times New Roman','fontsize',8);
xlabel('\fontsize{8}\fontname{Times new roman}Time (h)')
ylabel('\fontsize{8}\fontname{Times new roman}Power (kW)')

hold on 
subplot(2,1,2)
b=plot(temprature7(400:502,1),'r');
%axis([0,24,30000,90000]);
set(gca,'FontName','Times New Roman','FontSize',8)
title('Temperature','fontname','Times New Roman','fontsize',8);
xlabel('\fontsize{8}\fontname{Times new roman}Time (h)')
ylabel('\fontsize{8}\fontname{Times new roman}Temperature (℃)')


figure(19)
set(gcf,'unit','centimeters','position',[3,3,16,6])
subplot(1,2,1)
plot(T3,P3,'o','MarkerSize',2);
axis([10,30,45000,65000]);
set(gca,'FontName','Times New Roman','FontSize',8)
title({'Original Data:';'Power-Temperature'},'fontname','Times New Roman','fontsize',8);
xlabel('\fontsize{8}\fontname{Times new roman}Outdoor Temperature (℃)')
ylabel('\fontsize{8}\fontname{Times new roman}Load Power (kW)')

hold on
subplot(1,2,2)
plot([temp3weekday;temp3weekend],[P3baseloadweekday;P3baseloadweekend],'o','MarkerSize',2);
axis([10,30,45000,65000]);
set(gca,'FontName','Times New Roman','FontSize',8)
title({'Identified Data:';'Power-Temperature'},'fontname','Times New Roman','fontsize',8);
xlabel('\fontsize{8}\fontname{Times new roman}Outdoor Temperature (℃)')
ylabel('\fontsize{8}\fontname{Times new roman}Load Power (kW)')