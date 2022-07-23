%% 
%parameters
% %peptide that binds Bfl-1 with <=0.3nM Kd, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3631442/
% kon1 = 49e3; %1/M/s
% koff1 = 0.015e-3; %1/s
% %peptide that binds Mcl-1 with 6nM Kd
% kon2 = 49e3; %1/M/s
% koff2 = 1.5e-3; %1/s
% %peptide that binds Bcl-xL with 33nM Kd
% kon3 = 150e3; %1/M/s
% koff3 = 5e-3; %1/s
% %Fab, 11nM Kd
% kon4 = 2.3e5; %1/M/s
% koff4 = 2.5e-3; %1/s
% %Fab, 26nM Kd, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4697190/
% kon4 = 2.9e5; %1/M/s
% koff4 = 7.5e-3; %1/s
% %Fab, 133nM Kd
% kon4 = 6e4; %1/M/s
% koff4 = 8.1e-3; %1/s
% %Fab, 38nM Kd
% kon4 = 3.0e5; %1/M/s
% koff4 = 1.1e-2; %1/s
% %Fab, 89nM Kd
% kon4 = 2.7e5; %1/M/s
% koff4 = 2.4e-2; %1/s
%BIM-Mcl-1, 160pM Kd https://www.cell.com/fulltext/S0092-8674%2814%2900613-8#supplementaryMaterial
% kon1 = 317000; %1/M/s
% koff1 = 0.0000509; %1/s
kon1 = 31700; %1/M/s
koff1 = 0.0000509; %1/s
kon2 = 31700; %1/M/s
koff2 = 0.000509; %1/s
%BIM-Bfl-1, 544pM
% kon2 = 213000;
% koff2 = 0.000116;
% %BIM-Bcl-xL, 1.48nM
% kon3 = 221000;
% koff3 = 0.000326;
% %BIM-Bcl-w, 19.3nM
% kon4 = 122000;
% koff4 = 0.00235;
% %BIM-Bcl-2, 710pM
% kon5 = 183000;
% koff5 = 0.000130;
% P1 = 0e-9; %M
P2 = 100e-9; %M
%% 

%initial values
L_0 = 1e10/(6.022e23*25e-6); %M
LP1_0 = 0; %M
LP2_0 = 0; %M

%start time
t_0 = 0; %s
t_end = 3600*25; %s

%number of time steps
n=1000;

%time
time = linspace(t_0,t_end,5001);
%% plot + calculate LP1
figure
hold on
concs = [0,1,2.51,6.31,15.9,39.8,100];
Legend = cell(7,1);
for i = 1:7
    P1 = concs(i)*1e-9;
    [x,c] = ode15s(@(x,c)coupled(x,c,kon1,koff1,kon2,koff2,P1,P2),...
    time,[L_0,LP1_0,LP2_0]);

    if i == 1
        %LP1
        plot(x/60,c(:,2)/L_0)
        Legend{1} = sprintf('LP1, L1= %.2fnM',concs(i));

        %label axes
        xlabel('time (minutes)')
        ylabel('fraction target bound')

        %set graph limits
        axis([0 t_end/60 0 1])
        
    else
        plot(x/60,c(:,2)/L_0)
        Legend{i} = sprintf('LP1, L1= %.2fnM',concs(i));
    end
end
 
legend(Legend,'location','northwest')
set(gca,'Xscale','log')

%% plot + calculate LP2
figure
hold on
concs = [0,1,2.51,6.31,15.9,39.8,100];
Legend = cell(7,1);
for i = 1:7
    P1 = concs(i)*1e-9;
    [x,c] = ode15s(@(x,c)coupled(x,c,kon1,koff1,kon2,koff2,P1,P2),...
    time,[L_0,LP1_0,LP2_0]);

    if i == 1
        %LP2
        plot(x/60,c(:,3)/L_0)    
        Legend{1} = sprintf('LP2');

        %label axes
        xlabel('time (minutes)')
        ylabel('fraction off-target bound')

        %set graph limits
        axis([0 t_end/60 0 1])
        
    else
        plot(x/60,c(:,3)/L_0)
        Legend{i} = sprintf('LP2, L1= %.2fnM',concs(i));
    end
end
 
legend(Legend,'location','northwest')
set(gca,'Xscale','log')

%comment any of the following 'plot' commands to hide a graph (and
%change the legend accordingly)



% %% phase 1
% t=linspace(0,600,5000);
% LP1_LP2 = (1-exp(-(kon1*P1+kon2*P2)*t))*(kon1/kon2);
% figure
% xlim([t(1) t(end)])
% plot(t/60,LP1_LP2)
% set(gca,'XScale','log')
% 
% %% phase 2
% t1=5/2*(log(2)/(kon1*P1)+log(2)/(kon2*P2));
% t=linspace(t1,3600*25,5000);
% LP1_LP2 = kon1/kon2 + (koff2*kon1/(koff1*kon2)-kon1/kon2)*(1-exp(-(koff2-koff1)*(t-t1)));
% figure
% xlim([time(1) time(end)])
% plot(t/60,LP1_LP2)
% set(gca,'XScale','log')
%% 

function dcdt = coupled(x,c,kon1,koff1,kon2,koff2,P1,P2)
%differential equations governing peptide behavior
dcdt = [-kon1*c(1)*P1+koff1*c(2)-kon2*c(1)*P2+koff2*c(3);...
            kon1*c(1)*P1-koff1*c(2);...
            kon2*c(1)*P2-koff2*c(3)];
end