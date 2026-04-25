clc;
clear;
close all;

data= readmatrix('est_1.csv');

P_lin = data(:,1);
P_H= data(:,2);
P_L = data(:,3);

Lin1 = data(:,4);
Lin2= data(:,5);
Lin3 = data(:,6);

H1 = data(:,7);
H2= data(:,8);
H3 = data(:,9);

L1 = data(:,10);
L2= data(:,11);
L3 = data(:,12);


H = 21;
h=linspace(0,H-1,H);
zero=zeros(1,H);



figure('Units','normalized','Position',[0.2,0.2,.6,.2])

subplot(1,3,1)

        hold on
        f1 = plot(h , P_lin', '-ok', 'Linewidth',1.5);
        hold on
        f2 =plot(h , P_H', '-^b', 'Linewidth',1.5);
        hold on
        f3 =plot(h , P_L', 's-r', 'Linewidth',1.5);
        hold on
        
        plot(h, zero,'-k', 'Linewidth',0.1);
        
                

        n=get(gca,'ytick');
        xlim([-0.1 20])
        ylim([-0.5 1.5])
        set(gca,'XTick', [0 5 10 15 20])
        
        box off
        hold off
        grid on


        legend([f1 f2 f3],{'Linear','High' ,'Low'}, 'FontSize', 10, 'Location', 'northwest');
        xlabel('Horizon (Quarter)','FontSize', 11)    
        title('Three Models','FontSize', 15)
        
    
        
subplot(1,3,2)        
        
        IRmed = Lin2;
        IRlo= Lin1;
        IRup= Lin3;
        
        patch([h fliplr(h)], [IRup',fliplr(IRlo')],[0.85 0.85 0.85], 'EdgeColor',[1 1 1]);     
        hold on
        f2 = plot(h , IRmed', '-ok', 'Linewidth',1.5);
        hold on
        plot(h , IRup', '--k', 'Linewidth',1.5);
        hold on
        plot(h , IRlo', '--k', 'Linewidth',1.5);
        hold on
        plot(h, zero,'-k', 'Linewidth',0.1);
        
         xlim([-0.1 20])
        %ylim([-2.5 3.5])
        set(gca,'XTick', [0 5 10 15 20])
        
        box off
        hold off
        grid on


        %legend([f1 f2 f3],{'Linear','High' ,'Low'}, 'FontSize', 10, 'Location', 'northwest');
        xlabel('Horizon (Quarter)','FontSize', 11)    
        title('Linear Model','FontSize', 15)
        
subplot(1,3,3)        
        
        IRmed = H2;
        IRlo= H1;
        IRup= H3;
        


        patch([h fliplr(h)], [IRup',fliplr(IRlo')],[0.85 0.85 1], 'EdgeColor',[1 1 1]);     
        hold on
        f2 = plot(h , IRmed', '-^b', 'Linewidth',1.5);
        hold on
        plot(h , IRup', '--b', 'Linewidth',1.5);
        hold on
        plot(h , IRlo', '--b', 'Linewidth',1.5);
        hold on
       
        
        IRmed = L2;
        IRlo= L1;
        IRup= L3;
        


        %patch([h fliplr(h)], [IRup',fliplr(IRlo')],[1.00 0.85 0.85], 'EdgeColor',[1 1 1]);     
        hold on
        f2 = plot(h , IRmed', 's-r', 'Linewidth',1.5);
        hold on
        plot(h , IRup', ':r', 'Linewidth',2.5);
        hold on
        plot(h , IRlo', ':r', 'Linewidth',2.5);
        hold on
        plot(h, zero,'-k', 'Linewidth',0.1);
        
       
        
        xlim([-0.1 20])
        %ylim([-2.5 3.5])
        set(gca,'XTick', [0 5 10 15 20])
        
        box off
        hold off
        grid on


        %legend([f1 f2 f3],{'Linear','High Connectedness' ,'Low Connectedness'}, 'FontSize', 10, 'Location', 'northwest');
         xlabel('Horizon (Quarter)','FontSize', 11)  
        title('State-Dependent Models','FontSize', 15)
        
        
        set(gcf,'PaperPositionMode','auto')
            print -depsc2 gdp_base.eps
            
%% PCE

data= readmatrix('est_2.csv');

P_lin = data(:,1);
P_H= data(:,2);
P_L = data(:,3);

Lin1 = data(:,4);
Lin2= data(:,5);
Lin3 = data(:,6);

H1 = data(:,7);
H2= data(:,8);
H3 = data(:,9);

L1 = data(:,10);
L2= data(:,11);
L3 = data(:,12);


H = 21;
h=linspace(0,H-1,H);
zero=zeros(1,H);


figure('Units','normalized','Position',[0.2,0.2,.6,.2])    

subplot(1,3,1)

        hold on
        f1 = plot(h , P_lin', '-ok', 'Linewidth',1.5);
        hold on
        f2 =plot(h , P_H', '-^b', 'Linewidth',1.5);
        hold on
        f3 =plot(h , P_L', 's-r', 'Linewidth',1.5);
        hold on
        
        plot(h, zero,'-k', 'Linewidth',0.1);
        
                

        n=get(gca,'ytick');
        xlim([-0.1 20])
        %ylim([-2.5 3.5])
        set(gca,'XTick', [0 5 10 15 20])
        
        box off
        hold off
        grid on


        %legend([f1 f2 f3],{'Linear','High Connectedness' ,'Low Connectedness'}, 'FontSize', 10, 'Location', 'northwest');
         xlabel('Horizon (Quarter)','FontSize', 11)   
        title('Three Models','FontSize', 15)
        
    
        
subplot(1,3,2)        
        
        IRmed = Lin2;
        IRlo= Lin1;
        IRup= Lin3;
        


        patch([h fliplr(h)], [IRup',fliplr(IRlo')],[0.85 0.85 0.85], 'EdgeColor',[1 1 1]);     
        hold on
        f2 = plot(h , IRmed', '-ok', 'Linewidth',1.5);
        hold on
        plot(h , IRup', '--k', 'Linewidth',1.5);
        hold on
        plot(h , IRlo', '--k', 'Linewidth',1.5);
        hold on
        plot(h, zero,'-k', 'Linewidth',0.1);
        
         xlim([-0.1 20])
         ylim([-0.5 1])
        set(gca,'XTick', [0 5 10 15 20])
        
        box off
        hold off
        grid on


        %legend([f1 f2 f3],{'Linear','High Connectedness' ,'Low Connectedness'}, 'FontSize', 10, 'Location', 'northwest');
        xlabel('Horizon (Quarter)','FontSize', 11)      
        title('Linear Model','FontSize', 15)
        
subplot(1,3,3)        
        
        IRmed = H2;
        IRlo= H1;
        IRup= H3;
        


        patch([h fliplr(h)], [IRup',fliplr(IRlo')],[0.85 0.85 1], 'EdgeColor',[1 1 1]);     
        hold on
        f2 = plot(h , IRmed', '-^b', 'Linewidth',1.5);
        hold on
        plot(h , IRup', '--b', 'Linewidth',1.5);
        hold on
        plot(h , IRlo', '--b', 'Linewidth',1.5);
        hold on
       
        
        IRmed = L2;
        IRlo= L1;
        IRup= L3;
        


        %patch([h fliplr(h)], [IRup',fliplr(IRlo')],[1.00 0.85 0.85], 'EdgeColor',[1 1 1]);     
        hold on
        f2 = plot(h , IRmed', 's-r', 'Linewidth',1.5);
        hold on
        plot(h , IRup', ':r', 'Linewidth',2.5);
        hold on
        plot(h , IRlo', ':r', 'Linewidth',2.5);
        hold on
        plot(h, zero,'-k', 'Linewidth',0.1);
               
        
        xlim([-0.1 20])
        ylim([-0.5 2.0])
        set(gca,'XTick', [0 5 10 15 20])
        
        box off
        hold off
        grid on


        %legend([f1 f2 f3],{'Linear','High Connectedness' ,'Low Connectedness'}, 'FontSize', 10, 'Location', 'northwest');
        xlabel('Horizon (Quarter)','FontSize', 11)     
        title('State-Dependent Models','FontSize', 15)
        
        
        set(gcf,'PaperPositionMode','auto')
            print -depsc2 cpi_base.eps
            
               
%% FFR

data= readmatrix('est_3.csv');

P_lin = data(:,1);
P_H= data(:,2);
P_L = data(:,3);

Lin1 = data(:,4);
Lin2= data(:,5);
Lin3 = data(:,6);

H1 = data(:,7);
H2= data(:,8);
H3 = data(:,9);

L1 = data(:,10);
L2= data(:,11);
L3 = data(:,12);



H = 21;
h=linspace(0,H-1,H);
zero=zeros(1,H);



figure('Units','normalized','Position',[0.2,0.2,.6,.2])    

subplot(1,3,1)

        hold on
        f1 = plot(h , P_lin', '-ok', 'Linewidth',1.5);
        hold on
        f2 =plot(h , P_H', '-^b', 'Linewidth',1.5);
        hold on
        f3 =plot(h , P_L', 's-r', 'Linewidth',1.5);
        hold on
        
        plot(h, zero,'-k', 'Linewidth',0.1);
        
                

        n=get(gca,'ytick');
        xlim([-0.1 20])
        %ylim([-2.5 3.5])
        set(gca,'XTick', [0 5 10 15 20])
        
        box off
        hold off
        grid on


        %legend([f1 f2 f3],{'Linear','High Connectedness' ,'Low Connectedness'}, 'FontSize', 10, 'Location', 'northwest');
         xlabel('Horizon (Quarter)','FontSize', 11)   
        title('Three Models','FontSize', 15)
        
    
        
subplot(1,3,2)        
        
        IRmed = Lin2;
        IRlo= Lin1;
        IRup= Lin3;
        


        patch([h fliplr(h)], [IRup',fliplr(IRlo')],[0.85 0.85 0.85], 'EdgeColor',[1 1 1]);     
        hold on
        f2 = plot(h , IRmed', '-ok', 'Linewidth',1.5);
        hold on
        plot(h , IRup', '--k', 'Linewidth',1.5);
        hold on
        plot(h , IRlo', '--k', 'Linewidth',1.5);
        hold on
        plot(h, zero,'-k', 'Linewidth',0.1);
        
         xlim([-0.1 20])
        %ylim([-2.5 3.5])
        set(gca,'XTick', [0 5 10 15 20])
        
        box off
        hold off
        grid on


        %legend([f1 f2 f3],{'Linear','High Connectedness' ,'Low Connectedness'}, 'FontSize', 10, 'Location', 'northwest');
        xlabel('Horizon (Quarter)','FontSize', 11)      
        title('Linear Model','FontSize', 15)
        
subplot(1,3,3)        
        
        IRmed = H2;
        IRlo= H1;
        IRup= H3;
        


        patch([h fliplr(h)], [IRup',fliplr(IRlo')],[0.85 0.85 1], 'EdgeColor',[1 1 1]);     
        hold on
        f2 = plot(h , IRmed', '-^b', 'Linewidth',1.5);
        hold on
        plot(h , IRup', '--b', 'Linewidth',1.5);
        hold on
        plot(h , IRlo', '--b', 'Linewidth',1.5);
        hold on
       
        
        IRmed = L2;
        IRlo= L1;
        IRup= L3;
        


        %patch([h fliplr(h)], [IRup',fliplr(IRlo')],[1.00 0.85 0.85], 'EdgeColor',[1 1 1]);     
        hold on
        f2 = plot(h , IRmed', 's-r', 'Linewidth',1.5);
        hold on
        plot(h , IRup', ':r', 'Linewidth',2.5);
        hold on
        plot(h , IRlo', ':r', 'Linewidth',2.5);
        hold on
        plot(h, zero,'-k', 'Linewidth',0.1);
        
       
        
        xlim([-0.1 20])
        ylim([-2 1])
        set(gca,'XTick', [0 5 10 15 20])
        
        box off
        hold off
        grid on


        %legend([f1 f2 f3],{'Linear','High Connectedness' ,'Low Connectedness'}, 'FontSize', 10, 'Location', 'northwest');
        xlabel('Horizon (Quarter)','FontSize', 11)     
        title('State-Dependent Models','FontSize', 15)
        
        
        set(gcf,'PaperPositionMode','auto')
            print -depsc2 ffr_base.eps