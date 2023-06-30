% Created 2023-07-01
% Authors: Jun Steed Huang et al
% Analysis the RdRp sequence for COVID-19
  clc;
  clear;
  clf;
% read data from Excelsheet
  data=xlsread('NewProteaseAnalysisV4.xls', 3, 'A3:CL999');
% save data
% load('data.mat')

nameVirus={'CoV1','CoV2','Bat','BA5','XBB','Omicron','Murine','Mers','OC43','F229E','NL63','HKU1','RaTG13','Pangolin','MurineHep','PorcRes','PorcGas','Canine'};
[numAcid,numVirus]=size(data);
numVirus = numVirus/5;

for i=1:numVirus
    % extract the data    
    colStart = i*5-4;
    colEnd = i*5-4+3;
    oneVirus = getVirusFourLabel(data(:,colStart:colEnd));
    
    % b=[nameVirus{1,i},num2str(i)];
    nameV=[nameVirus{1,i}];
    eval([nameV,'=oneVirus']);
    % dynamic rename data
    clear oneVirus;
end

% Run equivalent inhibitor distance from 0.1 to 10 nm with 1pm step 
i = 1;
clear record
for z = 0.1:0.01:10
    recordZ(i)= z;
    record1(i)= log(Langlands(CoV1,z));
    record2(i)= log(Langlands(CoV2,z));
    record3(i)= log(Langlands(Bat,z));
    record4(i)= log(Langlands(BA5,z));
    record5(i)= log(Langlands(XBB,z));
    record6(i)= log(Langlands(Omicron,z));
    record7(i)= log(Langlands(Murine,z));
    record8(i)= log(Langlands(Mers,z));
    record9(i)= log(Langlands(OC43,z));
    record10(i)= log(Langlands(F229E,z));
    record11(i)= log(Langlands(NL63,z));
    record12(i)= log(Langlands(HKU1,z));
    record13(i)= log(Langlands(RaTG13,z));
    record14(i)= log(Langlands(Pangolin,z));
    record15(i)= log(Langlands(MurineHep,z));
    record16(i)= log(Langlands(PorcRes,z));
    record17(i)= log(Langlands(PorcGas,z));
    record18(i)= log(Langlands(Canine,z));
    % next virus 
    i=i+1;
end

% Calulate the distance among each others
record=(record1+record2+record3+record4+record5++record6+record7+record8+record9+record10+record11+record12+record13+record14++record15+record16+record17+record18)/18;
semi(1)=(mean((record1-record).^(1/18)))^18;
semi(2)=(mean((record2-record).^(1/18)))^18;
semi(3)=(mean((record3-record).^(1/18)))^18;
semi(4)=(mean((record4-record).^(1/18)))^18;
semi(5)=(mean((record5-record).^(1/18)))^18;
semi(6)=(mean((record6-record).^(1/18)))^18;
semi(7)=(mean((record7-record).^(1/18)))^18;
semi(8)=(mean((record8-record).^(1/18)))^18;
semi(9)=(mean((record9-record).^(1/18)))^18;
semi(10)=(mean((record10-record).^(1/18)))^18;
semi(11)=(mean((record11-record).^(1/18)))^18;
semi(12)=(mean((record12-record).^(1/18)))^18;
semi(13)=(mean((record13-record).^(1/18)))^18;
semi(14)=(mean((record14-record).^(1/18)))^18;
semi(15)=(mean((record15-record).^(1/18)))^18;
semi(16)=(mean((record16-record).^(1/18)))^18;
semi(17)=(mean((record17-record).^(1/18)))^18;
semi(18)=(mean((record18-record).^(1/18)))^18;

% The real of spectrum
figure(1)
loglog(recordZ,real(record1),recordZ,real(record2),'--',recordZ,real(record3),recordZ,real(record4),'--',recordZ,real(record5),recordZ,real(record6),'--',recordZ,real(record7),recordZ,real(record8),'--',recordZ,real(record9),recordZ,real(record10),'--',recordZ,real(record11),recordZ,real(record12),'--',recordZ,real(record13),recordZ,real(record14),'--',recordZ,real(record15),'-.',recordZ,real(record16),'-.',recordZ,real(record17),'-.',recordZ,real(record18),'-.','LineWidth',1);
legend(nameVirus);
title('Real of Equivalent Mass Charge vs Inhibitor Distance');
xlabel('inhibitor distance');
ylabel('mass charge ratio');
set(gcf,'Position',[0 50 900 400]);
% set(gca,'xtick',0.1:100:10) 
set(gca,'position',[0.04,0.14,0.94,0.8] );
set(gca,'FontSize',12);
print('-dtiff','-r600','Fig 1 real effects of different Z of the Langlands program on viruses'); 

% The imag of spectrum
figure(2)
loglog(recordZ,imag(record1),recordZ,imag(record2),'--',recordZ,imag(record3),recordZ,imag(record4),'--',recordZ,imag(record5),recordZ,imag(record6),'--',recordZ,imag(record7),recordZ,imag(record8),'--',recordZ,imag(record9),recordZ,imag(record10),'--',recordZ,imag(record11),recordZ,imag(record12),'--',recordZ,imag(record13),recordZ,imag(record14),'--',recordZ,imag(record15),'-.',recordZ,imag(record16),'-.',recordZ,imag(record17),'-.',recordZ,imag(record18),'-.','LineWidth',1);
legend(nameVirus);
title('Imag of Equivalent Mass Charge vs Inhibitor Distance');
xlabel('inhibitor distance');
ylabel('mass charge ratio');
set(gcf,'Position',[50 100 900 400]);
% set(gca,'xtick',0.1:100:10)
set(gca,'position',[0.04,0.14,0.94,0.8] );
set(gca,'FontSize',12);
print('-dtiff','-r600','Fig 2 imag effects of different Z of the Langlands program on viruses'); 

% The kin distance in between
figure(3)
plot(log(semi(1)),'+');
% pause;
hold on;
for j = 2:5
plot(log(semi(j)),'+');
% pause;
end
for j = 6:10
plot(log(semi(j)),'d');
% pause;
end
for j = 11:15
plot(log(semi(j)),'s');
% pause;
end
for j = 16:18
plot(log(semi(j)),'o');
% pause;
end
hold off;
legend(nameVirus);
title('Semi Kin Distance');
xlabel('conservative momentum');
ylabel('divertive momentum');
set(gcf,'Position',[100 150 900 400]);
% set(gca,'xtick',0.1:100:10)
set(gca,'position',[0.04,0.14,0.94,0.8] );
set(gca,'FontSize',12);
print('-dtiff','-r600','Fig 3 kin distance effects of different Z of the Langlands program on viruses'); 