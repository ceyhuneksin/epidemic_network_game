%% Author: Ceyhun Eksin
% MIT License
% Copyright (c) 2016, Ceyhun Eksin
%% code for getting figures in pnas_disease_dynamics_network
% Uses output from epidemic_game_main.m, R_x_computation.m,
% R_0_computation.m

%% Figure 3: Run epidemic_runs_star_network_demo.m
%% and then run the following to get the figures on Figure 3 of the manuscript
infected_demo_one = sum(x_store_bounded_information(:,:,1,1,1));
infected_demo_two = sum(x_store_bounded_information(:,:,1,2,1));

infected_demo_one = sum(x(:,:,1,1,1));
infected_demo_two = sum(x(:,:,1,2,1));

figure
hold on
plot([infected_demo_one 0],'LineStyle', '-', 'Color','r','MarkerSize',4, 'LineWidth', 1.2)
plot([infected_demo_two 0],'LineStyle', '-', 'Color','k','MarkerSize',4, 'LineWidth', 1.2)
lll= ylabel('Total number of infected','FontSize',14);
set(lll,'Interpreter','Latex');
xlabel('Time','FontSize',14)
axis([1 length(x(1,:,1,2,1))+1 0 4])
% lll = legend('Weak empathy \& Weak averseness','Strong empathy \& Weak averseness');
% set(lll,'Interpreter','Latex');
% set(lll,'FontSize',12)
% legend boxoff
% set(lll,'FontSize',12);
set(gca,'XTick',[1:1:15]);
set(gca,'YTick',[1:1:4]);


%%
figure
hold on
% plot(aggregate_utility(:,1,1),'LineStyle', '-', 'Color','r','MarkerSize',4, 'LineWidth', 1.2)
% plot(aggregate_utility(:,1,2),'LineStyle', '-', 'Color','k','MarkerSize',4, 'LineWidth', 1.2)
plot([u(:,1,1); 4],'LineStyle', '-', 'Color','r','MarkerSize',4, 'LineWidth', 1.2)
plot([u(:,1,2); 4],'LineStyle', '-', 'Color','k','MarkerSize',4, 'LineWidth', 1.2)

lll = ylabel('Aggregate utility','FontSize',12);
set(lll,'Interpreter','Latex');
xlabel('Time','FontSize',12)
axis([1 length(x(1,:,1,2,1))+1 0 4])
lll = legend('Weak empathy \& Weak averseness','Strong empathy \& Weak averseness');
% lll = legend('$3c_1>c_0 >2c_1$, $c_0 > 3c_2$','$3c_1>c_0>2c_1$, $c_0<3c_2$');
legend boxoff
set(lll,'Interpreter','Latex');
set(lll,'FontSize',16)
set(gca,'YTick',[1:1:4]);
set(gca,'XTick',[1:1:15]);


%% Figure 4 plots %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('R_0data_new_computation_c_1_024.mat')  % or run R_0_computation_new.m

number_infected_by_patient_zero = zeros(c_1_values,c_2_values,trials);
for betaval = 1:c_1_values
for c_2_value = 1:c_2_values
    for trial = 1:trials
        tt = 1; count = 0; sick_person = initial_sick_individual(betaval,c_2_value,trial);
        while x_store_bounded_information(sick_person,tt,betaval,c_2_value,trial) == 1
            who_was_not_sick = x_store_bounded_information(find(contact_network(sick_person,:)==1),tt,betaval,c_2_value,trial) == 0;
            who_is_sick = x_store_bounded_information(find(contact_network(sick_person,:)==1),tt+1,betaval,c_2_value,trial) == 1;
            who_got_sick = who_is_sick.*who_was_not_sick;
            count = count + sum(who_got_sick);
            tt = tt + 1
        end
        number_infected_by_patient_zero(betaval,c_2_value,trial) = count;
    end
end
end
range = c_2_values-1;

for beta_val = 1:1:3
figure
hold on

plot((0:range)*0.02,z_values(beta_val,(1:range+1)),'Color','k','MarkerSize',8, 'LineWidth', 1.5)

%R_0 upper bound
plot((0:range)*0.02,0.5*(beta_val)*step_size_beta*log(1+min(c_0./((0:range)*0.02),N))/delta,'LineStyle', '-.','Color','k','MarkerSize',8, 'LineWidth', 1.5)
plot(1/(exp(2*delta/(step_size_beta*beta_val))-1)*ones(6,1),0:0.2:1,'LineStyle', ':','Color',[0 0 1],'MarkerSize',8, 'LineWidth', 1.5)
plot((0:range)*0.02,ones(range+1,1),'LineStyle', ':','Color',[1 0 0],'MarkerSize',8, 'LineWidth', 1.5)
scatter(0,0.5*(beta_val)*step_size_beta*log(N-1)/delta,'filled','MarkerEdgeColor','r','MarkerFaceColor',(beta_val/4)*[1 0 0])
scatter(1/(exp(2*delta/(step_size_beta*beta_val))-1),1,'filled','MarkerEdgeColor','b','MarkerFaceColor',(beta_val/4)*[0 0 1])
%%%%%

lll = xlabel('Empathy constant ($c_2$)','FontSize',16);
set(lll,'Interpreter','Latex');
lll = ylabel('$R_0$','FontSize',16);
set(lll,'Interpreter','Latex');
axis([0 0.5 0 3.5])
lll = legend(['Simulation - $\beta =$ ',num2str(beta_val*step_size_beta)],['$R_0$ upper bound - $\beta =$ ',num2str(beta_val*step_size_beta)],'$c_2$ critical value', '$R_0=1$', '$R_0=\beta \log(n)/\delta$');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',16);
legend boxoff
end

lll = xlabel('Empathy constant ($c_2$)','FontSize',12);
set(lll,'Interpreter','Latex');
% use these my prints if you would like to transfer them to pdf.
% my_print('R_0_empirical_theoretical_beta_01')
% my_print('R_0_empirical_theoretical_beta_02')
% my_print('R_0_empirical_theoretical_beta_03')


%% Figure 5 plots %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('R_xdata_new_computation_c_1_024.mat') % or run R_x_computation_new.m

number_infected_by_patient_zero = zeros(c_1_values,c_2_values,trials);
for betaval = 1:c_1_values
for c_2_value = 1:c_2_values
    for trial = 1:trials
        tt = 1; count = 0; sick_person = initial_sick_individual(betaval,c_2_value,trial);
        while x_store_bounded_information(sick_person,tt,betaval,c_2_value,trial) == 1
            who_was_not_sick = x_store_bounded_information(find(contact_network(sick_person,:)==1),tt,betaval,c_2_value,trial) == 0;
            who_is_sick = x_store_bounded_information(find(contact_network(sick_person,:)==1),tt+1,betaval,c_2_value,trial) == 1;
            who_got_sick = who_is_sick.*who_was_not_sick;
            count = count + sum(who_got_sick);
            tt = tt + 1
        end
        number_infected_by_patient_zero(betaval,c_2_value,trial) = count;
    end
end
end
range = c_2_values-1;

for beta_val = 1:1:3
figure
hold on

plot((0:range)*0.02,z_values(beta_val,(1:range+1)),'Color','k','MarkerSize',8, 'LineWidth', 1.5)
%R_x upper bound
plot((2:range)*0.02,(beta_val)*step_size_beta*(min(c_0./((2:range)*0.02),N))/(delta*log(N)),'LineStyle', '-.','Color','k','MarkerSize',8, 'LineWidth', 1.5)
plot((beta_val)*step_size_beta/(delta*log(N))*ones(6,1),0:0.2:1,'LineStyle', ':','Color',[0 0 1],'MarkerSize',8, 'LineWidth', 1.5)
plot((0:range)*0.02,ones(range+1,1),'LineStyle', ':','Color',[1 0 0],'MarkerSize',8, 'LineWidth', 1.5)
scatter((beta_val)*step_size_beta/(delta*log(N)),1,'filled','MarkerEdgeColor','b','MarkerFaceColor',(beta_val/4)*[0 0 1])
%%%%%
lll = xlabel('Empathy constant ($c_2$)','FontSize',16);
set(lll,'Interpreter','Latex');
lll = ylabel('$R_*$','FontSize',16);
set(lll,'Interpreter','Latex');
axis([0 0.5 0 3.5])
lll = legend(['Simulation - $\beta =$ ',num2str(beta_val*step_size_beta)],['$R_*$ upper bound - $\beta =$ ',num2str(beta_val*step_size_beta)],'$c_2$ critical value', '$R_*=1$');
set(lll,'Interpreter','Latex');
set(lll,'FontSize',16);
legend boxoff
end

lll = xlabel('Empathy constant ($c_2$)','FontSize',12);
set(lll,'Interpreter','Latex');

% use these my prints if you would like to transfer them to pdf.
% my_print('R_x_empirical_theoretical_beta_01')
% my_print('R_x_empirical_theoretical_beta_02')
% my_print('R_x_empirical_theoretical_beta_03')


%% Figure 6 plots %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% rate = 0.1, all infected initially
load('c_1vsc_2beta01_scalefree_manynetworks.mat')
%%
for c_2_value = 1:c_2_values
    for c_1_value = 1:c_1_values
        for trial = 1:trials
            if store_eradication_time(c_1_value,c_2_value,trial) ==0
                store_eradication_time(c_1_value,c_2_value,trial) =200;
            end
        end
    end
end
mean_eradication_frequency = zeros(c_1_values,c_2_values);
for c_2_value = 1:c_2_values
for c_1_value = 1:c_1_values
mean_eradication_frequency(c_1_value,c_2_value) = sum(store_eradication_time(c_1_value,c_2_value,:)<200)/trials;
end
end
mean_eradication_frequency;
%%
figure
mymap = flip(parula,1);
hold on
imagesc(0:0.05:0.5,0:0.1:1,1-mean_eradication_frequency)
colormap(mymap);
imagesc(0:0.05:0.5,0:0.1:1,mean_eradication_frequency)
lll = ylabel('Risk averseness constant ($c_1$)','FontSize',16);
set(lll,'Interpreter','Latex');
lll = xlabel('Empathy constant ($c_2$)','FontSize',16);
set(lll,'Interpreter','Latex');
axis([0 0.5 0 1]) 

plot(1/(exp(2*delta/(rate))-1)*ones(6,1),0:0.2:1,'LineStyle', '-.','Color','w','MarkerSize',8, 'LineWidth', 2.5)
plot((rate)/(delta*log(N))*ones(6,1),0:0.2:1,'LineStyle', '-','Color','r','MarkerSize',8, 'LineWidth', 2.5)

% my_print('c_1vsc_2_color_fiftytrial_fullinfected_scale_free_beta_01')

%% rate = 0.2, all infected initially
load('c_1vsc_2beta02_scalefree_manynetworks.mat')
%%
for c_2_value = 1:c_2_values
    for c_1_value = 1:c_1_values
        for trial = 1:trials
            if store_eradication_time(c_1_value,c_2_value,trial) ==0
                store_eradication_time(c_1_value,c_2_value,trial) =200;
            end
        end
    end
end
mean_eradication_frequency = zeros(c_1_values,c_2_values);
for c_2_value = 1:c_2_values
for c_1_value = 1:c_1_values
mean_eradication_frequency(c_1_value,c_2_value) = sum(store_eradication_time(c_1_value,c_2_value,:)<200)/trials;
end
end
mean_eradication_frequency;
%%
figure
mymap = flip(parula,1);
hold on
imagesc(0:0.05:0.5,0:0.1:1,1-mean_eradication_frequency)
 colormap(mymap);
imagesc(0:0.05:0.5,0:0.1:1,mean_eradication_frequency)
lll = ylabel('Risk averseness constant ($c_1$)','FontSize',16);
set(lll,'Interpreter','Latex');
lll = xlabel('Empathy constant ($c_2$)','FontSize',16);
set(lll,'Interpreter','Latex');
axis([0 0.5 0 1]) 

plot(1/(exp(2*delta/(rate))-1)*ones(6,1),0:0.2:1,'LineStyle', '-.','Color','w','MarkerSize',8, 'LineWidth', 2.5)
plot((rate)/(delta*log(N))*ones(6,1),0:0.2:1,'LineStyle', '-','Color','r','MarkerSize',8, 'LineWidth', 2.5)

% my_print('c_1vsc_2_color_fiftytrial_fullinfected_scale_free_beta_02')

%% rate = 0.3, all infected initially
load('c_1vsc_2beta03_scalefree_manynetworks.mat')
%%
for c_2_value = 1:c_2_values
    for c_1_value = 1:c_1_values
        for trial = 1:trials
            if store_eradication_time(c_1_value,c_2_value,trial) ==0
                store_eradication_time(c_1_value,c_2_value,trial) =200;
            end
        end
    end
end
mean_eradication_frequency = zeros(c_1_values,c_2_values);
for c_2_value = 1:c_2_values
for c_1_value = 1:c_1_values
mean_eradication_frequency(c_1_value,c_2_value) = sum(store_eradication_time(c_1_value,c_2_value,:)<200)/trials;
end
end
mean_eradication_frequency;
%%
figure
mymap = flip(parula,1);
hold on
imagesc(0:0.05:0.5,0:0.1:1,1-mean_eradication_frequency)
colormap(mymap);
imagesc(0:0.05:0.5,0:0.1:1,mean_eradication_frequency)
c = colorbar;
set(c, 'ylim', [0 1])
c.Label.String = 'Eradication frequency';
set(c,'FontSize',16)
lll = ylabel('Risk averseness constant ($c_1$)','FontSize',16);
set(lll,'Interpreter','Latex');
lll = xlabel('Empathy constant ($c_2$)','FontSize',16);
set(lll,'Interpreter','Latex');
axis([0 0.5 0 1]) 

plot(1/(exp(2*delta/(rate))-1)*ones(6,1),0:0.2:1,'LineStyle', '-.','Color','w','MarkerSize',8, 'LineWidth', 2.5)
plot((rate)/(delta*log(N))*ones(6,1),0:0.2:1,'LineStyle', '-','Color','r','MarkerSize',8, 'LineWidth', 2.5)

% my_print('c_1vsc_2_color_fiftytrial_fullinfected_scale_free_beta_03')


%% rate = 0.1, single infected initially
load('c_1vsc_2beta01_scalefree_manynetworks_one_individual_infected.mat')
%%
for c_2_value = 1:c_2_values
    for c_1_value = 1:c_1_values
        for trial = 1:trials
            if store_eradication_time(c_1_value,c_2_value,trial) ==0
                store_eradication_time(c_1_value,c_2_value,trial) =200;
            end
        end
    end
end
mean_eradication_frequency = zeros(c_1_values,c_2_values);
for c_2_value = 1:c_2_values
for c_1_value = 1:c_1_values
mean_eradication_frequency(c_1_value,c_2_value) = sum(store_eradication_time(c_1_value,c_2_value,:)<200)/trials;
end
end
mean_eradication_frequency;
%%
figure
mymap = flip(parula,1);
hold on
imagesc(0:0.05:0.5,0:0.1:1,1-mean_eradication_frequency)
 colormap(mymap);
imagesc(0:0.05:0.5,0:0.1:1,mean_eradication_frequency)
lll = ylabel('Risk averseness constant ($c_1$)','FontSize',16);
set(lll,'Interpreter','Latex');
lll = xlabel('Empathy constant ($c_2$)','FontSize',16);
set(lll,'Interpreter','Latex');
axis([0 0.5 0 1]) 

plot(1/(exp(2*delta/(rate))-1)*ones(6,1),0:0.2:1,'LineStyle', '-.','Color','w','MarkerSize',8, 'LineWidth', 2.5)
plot((rate)/(delta*log(N))*ones(6,1),0:0.2:1,'LineStyle', '-','Color','r','MarkerSize',8, 'LineWidth', 2.5)

% my_print('c_1vsc_2_color_fiftytrial_singleinfected_scale_free_beta_01')


%% rate = 0.2, single infected initially
load('c_1vsc_2beta02_scalefree_manynetworks_one_individual_infected.mat')
%%
for c_2_value = 1:c_2_values
    for c_1_value = 1:c_1_values
        for trial = 1:trials
            if store_eradication_time(c_1_value,c_2_value,trial) ==0
                store_eradication_time(c_1_value,c_2_value,trial) =200;
            end
        end
    end
end
mean_eradication_frequency = zeros(c_1_values,c_2_values);
for c_2_value = 1:c_2_values
for c_1_value = 1:c_1_values
mean_eradication_frequency(c_1_value,c_2_value) = sum(store_eradication_time(c_1_value,c_2_value,:)<200)/trials;
end
end
mean_eradication_frequency;
%%
figure
mymap = flip(parula,1);
hold on
imagesc(0:0.05:0.5,0:0.1:1,1-mean_eradication_frequency)
colormap(mymap);
imagesc(0:0.05:0.5,0:0.1:1,mean_eradication_frequency)
lll = ylabel('Risk averseness constant ($c_1$)','FontSize',16);
set(lll,'Interpreter','Latex');
lll = xlabel('Empathy constant ($c_2$)','FontSize',16);
set(lll,'Interpreter','Latex');
axis([0 0.5 0 1]) 

plot(1/(exp(2*delta/(rate))-1)*ones(6,1),0:0.2:1,'LineStyle', '-.','Color','w','MarkerSize',8, 'LineWidth', 2.5)
plot((rate)/(delta*log(N))*ones(6,1),0:0.2:1,'LineStyle', '-','Color','r','MarkerSize',8, 'LineWidth', 2.5)

% my_print('c_1vsc_2_color_fiftytrial_singleinfected_scale_free_beta_02')


%% rate = 0.3, single infected initially
load('c_1vsc_2beta03_scalefree_manynetworks_one_individual_infected.mat')
%%
for c_2_value = 1:c_2_values
    for c_1_value = 1:c_1_values
        for trial = 1:trials
            if store_eradication_time(c_1_value,c_2_value,trial) ==0
                store_eradication_time(c_1_value,c_2_value,trial) =200;
            end
        end
    end
end
mean_eradication_frequency = zeros(c_1_values,c_2_values);
for c_2_value = 1:c_2_values
for c_1_value = 1:c_1_values
mean_eradication_frequency(c_1_value,c_2_value) = sum(store_eradication_time(c_1_value,c_2_value,:)<200)/trials;
end
end
mean_eradication_frequency;
%%
figure
mymap = flip(parula,1);
hold on
imagesc(0:0.05:0.5,0:0.1:1,1-mean_eradication_frequency)
colormap(mymap);
imagesc(0:0.05:0.5,0:0.1:1,mean_eradication_frequency)
c = colorbar;
set(c, 'ylim', [0 1])
c.Label.String = 'Eradication frequency';
set(c,'FontSize',16)
lll = ylabel('Risk averseness constant ($c_1$)','FontSize',16);
set(lll,'Interpreter','Latex');
lll = xlabel('Empathy constant ($c_2$)','FontSize',16);
set(lll,'Interpreter','Latex');
axis([0 0.5 0 1]) 

plot(1/(exp(2*delta/(rate))-1)*ones(6,1),0:0.2:1,'LineStyle', '-.','Color','w','MarkerSize',8, 'LineWidth', 2.5)
plot((rate)/(delta*log(N))*ones(6,1),0:0.2:1,'LineStyle', '-','Color','r','MarkerSize',8, 'LineWidth', 2.5)

% my_print('c_1vsc_2_color_fiftytrial_singleinfected_scale_free_beta_03')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figures for the supplement %%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3 (Supplement)   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('c_1vsc_2beta01_scalefree_manynetworks.mat')

for c_2_value = 1:c_2_values
    for c_1_value = 1:c_1_values
        for trial = 1:trials
            if store_eradication_time(c_1_value,c_2_value,trial) ==0
                store_eradication_time(c_1_value,c_2_value,trial) =200;
            end
        end
    end
end

mean_eradication_time = zeros(c_1_values,c_2_values);

for c_2_value = 1:c_2_values
for c_1_value = 1:c_1_values
mean_eradication_time(c_1_value,c_2_value) = mean(store_eradication_time(c_1_value,c_2_value,:));
end
end
mean_eradication_time;

figure
hold on
mymap = parula;
imagesc(0.5,1,200)
imagesc(0,0,0)
colormap(mymap);
% c = colorbar;
imagesc([0:0.05:0.5],[0:0.1:1],mean_eradication_time)
lll = ylabel('Risk averseness constant ($c_1$)','FontSize',16);
% lll = ylabel('Infection rate ($\beta$)','FontSize',10);
set(lll,'Interpreter','Latex');
lll = xlabel('Empathy constant ($c_2$)','FontSize',16);
set(lll,'Interpreter','Latex');
% c.Label.String = 'Mean eradication time';
% set(c,'FontSize',16);
axis([0 0.5 0 1]) 
set(gca,'YDir','normal')

%%
load('c_1vsc_2beta02_scalefree_manynetworks.mat')

for c_2_value = 1:c_2_values
    for c_1_value = 1:c_1_values
        for trial = 1:trials
            if store_eradication_time(c_1_value,c_2_value,trial) ==0
                store_eradication_time(c_1_value,c_2_value,trial) =200;
            end
        end
    end
end

mean_eradication_time = zeros(c_1_values,c_2_values);

for c_2_value = 1:c_2_values
for c_1_value = 1:c_1_values
mean_eradication_time(c_1_value,c_2_value) = mean(store_eradication_time(c_1_value,c_2_value,:));
end
end
mean_eradication_time;


figure
hold on
mymap = parula;
imagesc(0.5,1,200)
imagesc(0,0,0)
colormap(mymap);
% c = colorbar;
imagesc([0:0.05:0.5],[0:0.1:1],mean_eradication_time)
lll = ylabel('Risk averseness constant ($c_1$)','FontSize',16);
% lll = ylabel('Infection rate ($\beta$)','FontSize',10);
set(lll,'Interpreter','Latex');
lll = xlabel('Empathy constant ($c_2$)','FontSize',16);
set(lll,'Interpreter','Latex');
% c.Label.String = 'Mean eradication time';
% set(c,'FontSize',16);
axis([0 0.5 0 1]) 
set(gca,'YDir','normal')

%% 
load('c_1vsc_2beta03_scalefree_manynetworks.mat')
%%
for c_2_value = 1:c_2_values
    for c_1_value = 1:c_1_values
        for trial = 1:trials
            if store_eradication_time(c_1_value,c_2_value,trial) ==0
                store_eradication_time(c_1_value,c_2_value,trial) =200;
            end
        end
    end
end

mean_eradication_time = zeros(c_1_values,c_2_values);

for c_2_value = 1:c_2_values
for c_1_value = 1:c_1_values
mean_eradication_time(c_1_value,c_2_value) = mean(store_eradication_time(c_1_value,c_2_value,:));
end
end
mean_eradication_time;


figure
hold on
mymap = parula;
imagesc(0.5,1,200)
imagesc(0,0,0)
colormap(mymap);
c = colorbar;
imagesc([0:0.05:0.5],[0:0.1:1],mean_eradication_time)
lll = ylabel('Risk averseness constant ($c_1$)','FontSize',16);
% lll = ylabel('Infection rate ($\beta$)','FontSize',10);
set(lll,'Interpreter','Latex');
lll = xlabel('Empathy constant ($c_2$)','FontSize',16);
set(lll,'Interpreter','Latex');
c.Label.String = 'Mean eradication time';
set(c,'FontSize',16);
axis([0 0.5 0 1]) 
set(gca,'YDir','normal')

%% Figure 4 (Supplement)   %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('c_1vsc_2beta01_scalefree_manynetworks.mat')

final_x = zeros(c_1_values,c_2_values,trials,5);
for c_2_value = 1:c_2_values
    for c_1_value = 1:c_1_values
        for trial = 1:trials
        final_x(c_1_value,c_2_value,trial,1) = sum(x_store_bounded_information(:,T,c_1_value,c_2_value,trial))/N;
        end
    end
end

final_x_averaged = zeros(c_1_values,c_2_values,2);
for c_2_value = 1:c_2_values
    for c_1_value = 1:c_1_values
        final_x_averaged(c_1_value,c_2_value,1) = mean(final_x(c_1_value,c_2_value,:,1));
    end
end

zz_values = zeros(c_1_values,c_2_values);

for c_2_value = 1:c_2_values
for c_1_value = 1:c_1_values
zz_values(c_1_value,c_2_value) = mean(final_x(c_1_value,c_2_value,:,1));
end
end
zz_values

%%
figure
hold on
mymap = parula;
imagesc(0:0.05:0.5,0:0.1:1,[0:0.1:1]'*[0:0.1:1])
colormap(mymap);
% c = colorbar;
imagesc(0:0.05:0.5,0:0.1:1,zz_values)

lll = ylabel('Risk averseness constant ($c_1$)','FontSize',16);
% lll = ylabel('Infection rate ($\beta$)','FontSize',10);
set(lll,'Interpreter','Latex');
lll = xlabel('Empathy constant ($c_2$)','FontSize',16);
set(lll,'Interpreter','Latex');
% c.Label.String = 'Ratio of infected';
axis([0 0.5 0 1]) 
set(gca,'YDir','normal')

%%
load('c_1vsc_2beta02_scalefree_manynetworks.mat')

final_x = zeros(c_1_values,c_2_values,trials,5);
for c_2_value = 1:c_2_values
    for c_1_value = 1:c_1_values
        for trial = 1:trials
        final_x(c_1_value,c_2_value,trial,1) = sum(x_store_bounded_information(:,T,c_1_value,c_2_value,trial))/N;
        end
    end
end

final_x_averaged = zeros(c_1_values,c_2_values,2);
for c_2_value = 1:c_2_values
    for c_1_value = 1:c_1_values
        final_x_averaged(c_1_value,c_2_value,1) = mean(final_x(c_1_value,c_2_value,:,1));
    end
end

zz_values = zeros(c_1_values,c_2_values);

for c_2_value = 1:c_2_values
for c_1_value = 1:c_1_values
zz_values(c_1_value,c_2_value) = mean(final_x(c_1_value,c_2_value,:,1));
end
end
zz_values

%%
figure
hold on
mymap = parula;
imagesc(0:0.05:0.5,0:0.1:1,[0:0.1:1]'*[0:0.1:1])
colormap(mymap);
% c = colorbar;
imagesc(0:0.05:0.5,0:0.1:1,zz_values)

lll = ylabel('Risk averseness constant ($c_1$)','FontSize',16);
% lll = ylabel('Infection rate ($\beta$)','FontSize',10);
set(lll,'Interpreter','Latex');
lll = xlabel('Empathy constant ($c_2$)','FontSize',16);
set(lll,'Interpreter','Latex');
% c.Label.String = 'Ratio of infected';
axis([0 0.5 0 1]) 
set(gca,'YDir','normal')

%%
load('c_1vsc_2beta03_scalefree_manynetworks.mat')
%%
final_x = zeros(c_1_values,c_2_values,trials,5);
for c_2_value = 1:c_2_values
    for c_1_value = 1:c_1_values
        for trial = 1:trials
        final_x(c_1_value,c_2_value,trial,1) = sum(x_store_bounded_information(:,T,c_1_value,c_2_value,trial))/N;
        end
    end
end

final_x_averaged = zeros(c_1_values,c_2_values,2);
for c_2_value = 1:c_2_values
    for c_1_value = 1:c_1_values
        final_x_averaged(c_1_value,c_2_value,1) = mean(final_x(c_1_value,c_2_value,:,1));
    end
end

zz_values = zeros(c_1_values,c_2_values);

for c_2_value = 1:c_2_values
for c_1_value = 1:c_1_values
zz_values(c_1_value,c_2_value) = mean(final_x(c_1_value,c_2_value,:,1));
end
end
zz_values

%%
figure
hold on
mymap = parula;
imagesc(0:0.05:0.5,0:0.1:1,[0:0.1:1]'*[0:0.1:1])
colormap(mymap);
c = colorbar;
imagesc(0:0.05:0.5,0:0.1:1,zz_values)

lll = ylabel('Risk averseness constant ($c_1$)','FontSize',16);
% lll = ylabel('Infection rate ($\beta$)','FontSize',10);
set(lll,'Interpreter','Latex');
lll = xlabel('Empathy constant ($c_2$)','FontSize',16);
set(lll,'Interpreter','Latex');
c.Label.String = 'Ratio of infected';
set(c,'FontSize',16);
axis([0 0.5 0 1]) 
set(gca,'YDir','normal')