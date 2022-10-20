clc; clear all; close all;

N_pid   = 1;  % # of CPU for parallel computing
N       = 512;
shape   = 1;   % 0: square; 1: rectangular
Oxy_sat = 1;
T_char  = 2.6719;
t_mix   = 100;
H_coeff = 1./30;

dt_File = 0.1123*10; % 0.1123; 0.0244
Ini_t   = 100.0593; % 50.018; 40.0038;
Fin_t   = 593; % 499.8473; % 20.0-dt_File;

% home = ['C:\0_Code\1_bioreactor\basilisk_new\basilisk\src\bioreactor\oxygen_test\oxy23\Data_all'];
% home = ['E:\1_bioreactor\Data_save\bio_rec_new3_save\oxy26\Data_all'];
home  = ['/gpfs/scratch/mkim79/basilisk/bio_rec_new3/oxy28/oxy28_save'];
home1 = ['/gpfs/scratch/mkim79/basilisk/bio_rec_new3/oxy28/'];
home_data = [home,'/Data_all_'];
% home = ['E:\1_bioreactor\Data_save\bio_rot_sq_contact\N=',num2str(N),'_th=',num2str(th),'\Data_all'];


%% Plot options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

in_plot1 = 0; % ux field vs. x
in_plot2 = 0; % Plot: max/min tracer
in_plot3 = 0; % Plot: max/min volume fraction
in_plot4 = 0; % Plot: error vs n

in = 0; in1 = 0;
File_t  = [Ini_t:dt_File:floor(Fin_t/dt_File)*dt_File];
for i = File_t
    i
    in    = in+1;
    t(in) = i*T_char; % s
    
    %% Get fields: x, y, velocity, volume fraction, concentration %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % x,y,r_grid[],th_grid[],u.x[],u.y[],ur[],utheta[],f[],c[]
    
    x_field = [];    y_field = [];    ux_field = [];   uy_field = [];
    al_field = [];   tr_field = [];   so_field = [];   ox_field = [];
    
    for j = 0:N_pid-1
        
        % fID_field  = fopen([home_data,num2str(N),'_',num2str(i),'_',num2str(j),'.txt'],'r');
        fID_field  = fopen([home_data,num2str(N),'_',num2str(i),'_tot.txt'],'r');
        % Data_field = fscanf(fID_field,'%f'); fclose(fID_field);
        
        if      shape == 0 % square
            Data_tmp = textscan(fID_field,'%f %f %f %f %f %f'); fclose(fID_field);
            x_tmp  = Data_tmp{1};
            y_tmp  = Data_tmp{2};
            ux_tmp = Data_tmp{3};
            uy_tmp = Data_tmp{4};
            al_tmp = Data_tmp{5};
            tr_tmp = Data_tmp{6};
            
        elseif shape == 1 % rectangular
            Data_tmp = textscan(fID_field,'%f %f %f %f %f %f %f %f'); fclose(fID_field);
            x_tmp  = Data_tmp{1};
            y_tmp  = Data_tmp{2};
            ux_tmp = Data_tmp{3};
            uy_tmp = Data_tmp{4};
            al_tmp = Data_tmp{5};
            tr_tmp = Data_tmp{6};
            so_tmp = Data_tmp{7};
            ox_tmp = Data_tmp{8};
            
            so_field = [so_field; so_tmp];
        end
        
        x_field = [x_field; x_tmp];
        y_field = [y_field; y_tmp];
        ux_field = [ux_field; ux_tmp];
        uy_field = [uy_field; uy_tmp];
        al_field = [al_field; al_tmp];
        tr_field = [tr_field; tr_tmp];
        ox_field = [ox_field; ox_tmp];
        
        clear Data_field
    end
    
    if shape == 1
        in_solid = find(so_field>0);
        x_field  = x_field(in_solid);
        y_field  = y_field(in_solid);
        ux_field = ux_field(in_solid);
        uy_field = uy_field(in_solid);
        al_field = al_field(in_solid);
        tr_field = tr_field(in_solid);
        ox_field = ox_field(in_solid);
    end
    
    %%% extract liquid/gas field
    in_liq_field  = find(al_field>0);
    in_gas_field = find(al_field<1);
    al_field_liq = al_field(in_liq_field);
    al_field_gas = al_field(in_gas_field);
    tr_field_liq = tr_field(in_liq_field);
    tr_field_gas = tr_field(in_gas_field);
    ox_field_liq = ox_field(in_liq_field);
    ox_field_gas = ox_field(in_gas_field);
    leng_gas     = numel(ox_field_gas);
    leng_liq     = numel(ox_field_liq);
    
    %%% strict extraction
    in_liq_field2 = find(abs(al_field)>1-1e-10);
    in_gas_field2 = find(abs(al_field)<1e-10);
    in_inf_field2 = find(abs(al_field)<1-1e-10 & abs(al_field)>1e-10);
    al_field_liq2 = al_field(in_liq_field2);
    al_field_gas2 = al_field(in_gas_field2);
    al_field_inf2 = al_field(in_inf_field2);
    tr_field_liq2 = tr_field(in_liq_field2);
    tr_field_gas2 = tr_field(in_gas_field2);
    tr_field_inf2 = tr_field(in_inf_field2);
    ox_field_liq2 = ox_field(in_liq_field2);
    ox_field_gas2 = ox_field(in_gas_field2);
    ox_field_inf2 = ox_field(in_inf_field2);
    leng_gas2     = numel(ox_field_gas2);
    leng_liq2     = numel(ox_field_liq2);
    
    %%% tracer in the liquid/gas field
    max_tr_gas(in) = max(abs(tr_field_gas));
    max_tr_liq(in) = max(abs(tr_field_liq));
    sum_tr_gas(in) = sum(abs(tr_field_gas));
    sum_tr_liq(in) = sum(abs(tr_field_liq));
    sum_tr_tot(in) = sum(tr_field);
    
    %%% strict: tracer in liquid/gas
    max_tr_gas2(in) = max(abs(tr_field_gas2));
    max_tr_liq2(in) = max(abs(tr_field_liq2));
    sum_tr_gas2(in) = sum(abs(tr_field_gas2));
    sum_tr_liq2(in) = sum(abs(tr_field_liq2));
    sum_tr_inf2_liq(in) = sum(tr_field_inf2.*al_field_inf2);
    sum_tr_inf2_gas(in) = sum(tr_field_inf2.*(1-al_field_inf2));
    sum_tr_tot2(in) = sum(tr_field);
    
    max_ox_gas2(in) = max(abs(ox_field_gas2));
    max_ox_liq2(in) = max(abs(ox_field_liq2));
    sum_ox_gas2(in) = sum(abs(ox_field_gas2));
    sum_ox_liq2(in) = sum(abs(ox_field_liq2));
    sum_ox_inf2_liq(in) = sum(ox_field_inf2.*al_field_inf2);
    sum_ox_inf2_gas(in) = sum(ox_field_inf2.*(1-al_field_inf2));
    sum_ox_tot2(in) = sum(ox_field);
    
    %     %%% across the interface
    %     in_interf = find(al_field > 0 & al_field < 1);
    %     tr_center = tr_field(find(abs((x_field-x_plot)) <= 1e-4));
    %     y_center  = y_field (find(abs((x_field-x_plot)) <= 1e-4));
    %
    %     out_center = [y_center tr_center];
    %     out_center_sort = sortrows(out_center,1);
    %     out_center_save(in,:,:) = out_center_sort;
    %
    %     if rem(in,plot_frq) == 1
    %         plot(out_center_sort(:,1),out_center_sort(:,2)); hold on
    %     end
    
    al_field_sort = sort(al_field);
    x_field_sort  = sort(x_field);
    y_field_sort  = sort(y_field);
    
    %%% when t > t_mix
    if i > t_mix
        in1 = in1+1;
        
        %% The degree of mixing
        dA     = (1/N)^2; % uniform grid
        
        %%% consider interface
        c_avg(in1)  = sum(tr_field_liq*dA)/(dA*leng_liq);
        c2_avg(in1) = sum(tr_field_liq.^2*dA)/(dA*leng_liq);
        sig(in1)    = c2_avg(in1) - c_avg(in1)^2;
        if in1 == 1  c_mix_max = sig(in1);  end
        c_mix(in1) = 1-sig(in1)/c_mix_max;
        
        %%% strict restriction
        c_avg2(in1)  = sum(tr_field_liq2*dA)/(dA*leng_liq2);
        c2_avg2(in1) = sum(tr_field_liq2.^2*dA)/(dA*leng_liq2);
        sig2(in1)   = c2_avg2(in1) - c_avg2(in1)^2;
        if in1 == 1  c_mix_max2 = sig2(in1);  end
        c_mix2(in1) = 1-sig2(in1)/c_mix_max2;
        
        %% Oxygen transfer
        dA       = (1/N)^2;
        A_tot    = dA*leng_liq;
        A_tot2   = dA*leng_liq2;
        ox_con_liq (in1) = sum(ox_field_liq *dA)/A_tot;
        ox_con_liq2(in1) = sum(ox_field_liq2*dA)/A_tot2;
        
        %% Oxygen transfer: volume simulation concentration
        % saturation concentration in liquid
        if in1 == 1
            A_tot_gas   = dA*leng_gas;
            ox_sat_sim  = sum(ox_field_gas *dA)/A_tot_gas *H_coeff;
            A_tot_gas2  = dA*leng_gas2;
            ox_sat_sim2 = sum(ox_field_gas2*dA)/A_tot_gas2*H_coeff;
        end
        if in1==1
            kLa_sim(in1)  = NaN;
            kLa_sim2(in1) = NaN;
        else
            kLa_sim(in1)  = log((ox_sat_sim  - ox_con_liq(1)) /(ox_sat_sim  - ox_con_liq (in1)))/((t(in)-t(1))/3600);
            kLa_sim2(in1) = log((ox_sat_sim2 - ox_con_liq2(1))/(ox_sat_sim2 - ox_con_liq2(in1)))/((t(in)-t(1))/3600);
        end
    end
end

out_cmix = [t' c_avg'  c2_avg' sig' c_mix' c_avg2' c2_avg2' sig2' c_mix2'];
out_oxy  = [t' ox_con_liq' kLa_sim' ox_con_liq2' kLa_sim2'];

out2 = [t' max_tr_gas2' max_tr_liq2' sum_tr_gas2' sum_tr_liq2' sum_tr_inf2_gas' sum_tr_inf2_liq' sum_tr_tot2'];
out3 = [t' max_ox_gas2' max_ox_liq2' sum_ox_gas2' sum_ox_liq2' sum_ox_inf2_gas' sum_ox_inf2_liq' sum_ox_tot2'];
out  = [out2 out3];

save([home1,'out_cmix.mat'],'out_cmix');
save([home1,'out_oxy.mat'] ,'out_oxy');
save([home1,'out_tr_oxy.mat'],'out');
    
%     %%% calculate the height of interface
%     x_left  = min(x_field);
%     x_right = max(x_field);
%     
%     in_l = find(x_field==x_left);
%     in_r = find(x_field==x_right);
%     
%     al_l = al_field(in_l);  % left alpha
%     al_r = al_field(in_r);  % right alpha
%     y_l  = y_field(in_l);
%     y_r  = y_field(in_r);
%     leng = numel(x_field);
% 
%     for j = 1:numel(al_l)-1
%         % t = 0
%         if     i==0 & al_l(j)==1   & al_l(j+1)==0 
%             y_fl(in) = 0.5*(y_l(j)+y_l(j+1)); break;
%         % t > 0
%         elseif i>0  & al_l(j)>1e-6 & al_l(j)<1-1e-6
%             y_fl(in) = 0.5*(y_l(j)+y_l(j+1)); break;
%         end        
%     end
%     for j = 1:numel(al_r)-1
%         if     i==0 & al_r(j)==1   & al_r(j+1)==0 
%             y_fr(in) = 0.5*(y_r(j)+y_r(j+1)); break;
%         % t > 0
%         elseif i>0  & al_r(j)>1e-6 & al_r(j)<1-1e-6
%             y_fr(in) = 0.5*(y_r(j)+y_r(j+1)); break;
%         end 
%     end
%     
%     if i == File_t(end)
%         figure,plot(File_t,y_fl,'-o',File_t,y_fr,'-d','linewidth',1);
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('Interface at the end'); xlabel('time'); ylabel('y');
%         legend('left end','right end');
%         f = gcf; exportgraphics(f,'time_vs_y_interf.png','Resolution',300);
%         
%         out1 = [File_t' y_fl' y_fr'];
%     end
%     
%     %%% extract liquid field
%     in_al_field  = find(al_field>0);
%     in_gas_field = find(al_field<1);
%     al_field_liq = al_field(in_al_field);
%     ux_field_liq = ux_field(in_al_field);
%     uy_field_liq = uy_field(in_al_field);
%     al_field_gas = al_field(in_al_field);
%     ux_field_gas = ux_field(in_al_field);
%     uy_field_gas = uy_field(in_al_field);
%     ox_field_liq = ox_field(in_al_field);
%     leng_liq     = numel(ox_field_liq);
%     
%     if in == 1 ox_field_gas = ox_field(in_gas_field); end
%     
%     %%% calculate average velocity of liquid
%     ux_ave_liq(in) = mean(ux_field_liq.*al_field_liq);
%     uy_ave_liq(in) = mean(uy_field_liq.*al_field_liq);
%     
%     %%% calculate average velocity of gas
%     ux_ave_gas(in) = mean(ux_field_gas.*al_field_gas);
%     uy_ave_gas(in) = mean(uy_field_gas.*al_field_gas);
%     
%     %%% calculate standard deviation of liquid
%     ux_dev_liq(in) = sqrt(mean(((ux_field_liq.*al_field_liq)-ux_ave_liq(in)).^2));
%     uy_dev_liq(in) = sqrt(mean(((uy_field_liq.*al_field_liq)-uy_ave_liq(in)).^2));
%     
%     %%% calculate standard deviation of gas
%     ux_dev_gas(in) = sqrt(mean(((ux_field_gas.*al_field_gas)-ux_ave_gas(in)).^2));
%     uy_dev_gas(in) = sqrt(mean(((uy_field_gas.*al_field_gas)-uy_ave_gas(in)).^2));
%     
%     %%% calculate maximum local velocity of liquid
%     ux_max_liq(in) = max(abs(ux_field_liq).*al_field_liq);
%     uy_max_liq(in) = max(abs(uy_field_liq).*al_field_liq);
% 
%     %%% calculate maximum local velocity of gas
%     ux_max_gas(in) = max(abs(ux_field_gas).*al_field_gas);
%     uy_max_gas(in) = max(abs(uy_field_gas).*al_field_gas);
%     
%     
%     %% The degree of mixing
%     dA     = (1/N)^2; % uniform grid
%     c_avg  = sum(tr_field*dA)/(dA*leng);
%     c2_avg = sum(tr_field.^2*dA)/(dA*leng);
%     
%     sig2   = c2_avg-c_avg^2;
%     if i == Ini_t  c_mix_max = sig2;  end
%     c_mix(in)  = 1-sig2/c_mix_max;
%     
%     %% Oxygen transfer
%     dA       = (1/N)^2;
%     A_tot    = dA*leng_liq;
%     ox_con_liq(in)     = sum(ox_field_liq*dA)/A_tot;
%     ox_con_vol_liq(in) = sum(ox_field_liq.*dA)/sum(ox_field_gas.*dA);
%     
%     %% Oxygen transfer: volume simulation concentration
%     oxy_sat_sim = 1;    % saturation concentration in liquid
%     if in==1  kLa_sim(in) = 0;
%     else      kLa_sim(in) = log((oxy_sat_sim-ox_con_liq(1))/(oxy_sat_sim-ox_con_liq(in)))/((t(in)-t(1))/3600); end
%     
%     %% Oxygen transfer: volume actual concentration: only oxygen
%     oxy_sat_vol = 0.034*5;    % volume for saturated oxygen in liquid
%     if in==1  kLa_vol(in) = 0;
%     else
%         if ox_con_vol_liq(in)/oxy_sat_vol < 0.99
%             kLa_vol(in) = log((oxy_sat_vol-ox_con_vol_liq(1))/(oxy_sat_vol-ox_con_vol_liq(in)))/((t(in)-t(1))/3600);
%         else
%             kLa_vol(in) = 0;
%         end
%     end
%     
%     %% Oxygen transfer: mass actual concentration: 2L water/2L oxygen
%     rho_g     = 1.2;        % 1.2g/L
%     Liq_vol   = 2;          % L
%     Gas_vol   = 2;          % 2/5; divided by 5 for oxygen
%     oxy_sat_mas = 8.2e-3*5; % mass for saturated oxygen in liquid
%     ox_con_mas_liq(in) = (rho_g*Gas_vol)*ox_con_vol_liq(in)/Liq_vol;  % mass concentration
%     if in==1  kLa_mas(in) = 0;
%     else
%         if ox_con_mas_liq(in)/oxy_sat_mas < 0.99
%             kLa_mas(in) = log((oxy_sat_mas-ox_con_mas_liq(1))/(oxy_sat_mas-ox_con_mas_liq(in)))/((t(in)-t(1))/3600);
%         else
%             kLa_mas(in) = 0;
%         end
%     end
%     
%     if i == File_t(end)
%         figure,plot(File_t,ux_ave_liq,'-o',File_t,uy_ave_liq,'-d','linewidth',1);
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('Mean velocity'); xlabel('time'); ylabel('mean vel');
%         legend('x','y');
%         f = gcf; exportgraphics(f,'time_vs_mean_vel.png','Resolution',300);
%         
%         figure,plot(File_t,ux_dev_liq,'-o',File_t,uy_dev_liq,'-d','linewidth',1);
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('standard deviation'); xlabel('time'); ylabel('mean vel');
%         legend('x','y');
%         
%         figure,plot(File_t,ux_max_liq,'-o',File_t,uy_max_liq,'-d','linewidth',1);
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('Max velocity'); xlabel('time'); ylabel('mean vel');
%         legend('x','y');
%         
%         figure,plot(File_t,ux_ave_gas,'-o',File_t,uy_ave_gas,'-d','linewidth',1);
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('Mean velocity gas'); xlabel('time'); ylabel('mean vel');
%         legend('x','y');
%         
%         figure,plot(File_t,c_mix,'-o','linewidth',1);
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('The degree of mixing');
%         xlabel('time'); ylabel('mixing');
%         f = gcf; exportgraphics(f,'time_vs_mixing.png','Resolution',300);
% 
%         figure,plot(File_t,ox_con_liq,'-o','linewidth',1);
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('Oxygen concentration: simul'); xlabel('time'); ylabel('oxy concentraiton');
%         
%         figure,plot(File_t,ox_con_vol_liq,'-o','linewidth',1);
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('Oxygen concentration: volume'); xlabel('time'); ylabel('oxy concentraiton');
%         
%         figure,plot(File_t,ox_con_mas_liq,'-o','linewidth',1);
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('Oxygen concentration: mass'); xlabel('time'); ylabel('oxy concentraiton');
%         
%         figure,plot(File_t,kLa_sim,'-o','linewidth',1);
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('Oxygen transfer: simul'); xlabel('time'); ylabel('kLa');
%         
%         figure,plot(File_t,kLa_vol,'-o','linewidth',1);
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('Oxygen transfer: volume'); xlabel('time'); ylabel('kLa');
%         
%         figure,plot(File_t,kLa_mas,'-o','linewidth',1);
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('Oxygen transfer: mass'); xlabel('time'); ylabel('kLa');
%         
%         out2 = [File_t' t' ux_ave_liq' uy_ave_liq' ux_ave_gas' uy_ave_gas' ...
%             ux_dev_liq' uy_dev_liq' ux_dev_gas' uy_dev_gas' ux_max_liq' uy_max_liq' ...
%             ux_max_gas' uy_max_gas' ox_con_liq' ox_con_vol_liq' ox_con_mas_liq' kLa_sim' kLa_vol' kLa_mas' c_mix'];
%     end

    %%% extract interface
%     Interf_x  = Dal_field(1:2:end);
%     Interf_y  = Dal_field(2:2:end);
%     Interf    = [Interf_x Interf_y];
    
    
%     %% Plot x- and y-velocity  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %%% theta=fixed
%     [tmp in_plotth] = min(abs(th_tar_plot-th_field));
%     th_plot = th_field(in_plotth);
%     
%     in_find  = find(th_field==th_plot);
%     uth_plot(:,1) = r_field(in_find);
%     uth_plot(:,2) = uth_field(in_find);
%     ur_plot_th0(:,1) = r_field(in_find);
%     ur_plot_th0(:,2) = ur_field(in_find);
%     uth_plot      = sortrows(uth_plot);
%     ur_plot_th0   = sortrows(ur_plot_th0);
%     uth_plot_s(in+1,:,:) = uth_plot;
%     ur_plot_th0_s(in+1,:,:) = ur_plot_th0;
%     
%     
%     %%% y=0
%     [tmp in_plotx] = min(abs(y_tar_plot-y_field));
%     y_plot = y_field(in_plotx);
%     
%     in_find  = find(y_field==y_plot);
%     uy_plot(:,1) = x_field(in_find);
%     uy_plot(:,2) = uy_field(in_find);
%     ux_plot_y0(:,1) = x_field(in_find);
%     ux_plot_y0(:,2) = ux_field(in_find);
%     uy_plot      = sortrows(uy_plot);
%     ux_plot_y0   = sortrows(ux_plot_y0);
%     uy_plot_s(in+1,:,:) = uy_plot;
%     ux_plot_y0_s(in+1,:,:) = ux_plot_y0;
%         
%     if acc == 1
%         [tmp in_plotx2] = min(abs(y_tar_plot-ya_field));
%         ya_plot = ya_field(in_plotx2);
%         
%         in_find  = find(y_field==y_plot);        
%         ax_plot(:,1) = xa_field(in_find);
%         ay_plot(:,1) = xa_field(in_find);
%         ax_plot(:,2) = ax_field(in_find);
%         ay_plot(:,2) = ay_field(in_find);
%         ax_plot(:,3) = cox_field(in_find);
%         ay_plot(:,3) = coy_field(in_find);
%         ax_plot(:,4) = cex_field(in_find);
%         ay_plot(:,4) = cey_field(in_find);
%         ax_plot(:,5) = uxa_field(in_find);
%         ay_plot(:,5) = uya_field(in_find);
%         ax_plot(:,6) = uxa2_field(in_find);
%         ay_plot(:,6) = uya2_field(in_find);        
%         ax_plot      = sortrows(ax_plot);
%         ay_plot      = sortrows(ay_plot);
%         ax_plot_s(in+1,:,:) = ax_plot;
%         ay_plot_s(in+1,:,:) = ay_plot;
%     end
%     
%     %%% x=0
%     [tmp in_ploty] = min(abs(x_tar_plot-x_field));
%     x_plot = x_field(in_ploty);
%         
%     in_find  = find(x_field==x_plot);
%     ux_plot(:,1) = y_field(in_find);
%     ux_plot(:,2) = ux_field(in_find);
%     ux_plot   = sortrows(ux_plot);
%     ux_plot_s(in+1,:,:) = ux_plot;
    
    
%     %% Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if in_plot1 == 1
%         
%         % Interface
%         figure(1); plot(Interf_x,Interf_y,'o'); hold on; legend('1','2','3','4','5','6','7');
%         
%         % u_theta - 45 degree
%         figure(2); plot(uth_plot_s(in,:,1),uth_plot_s(in,:,2),'-o','linewidth',1); hold on
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('uth field vs r; th=45'); xlabel('x'); ylabel('uth');
%         legend('1','2','3','4','5','6','7');
%         f = gcf; exportgraphics(f,'uth_vs_r.png','Resolution',300); 
%         disp('Plot: uth_vs_r');
%         
%         % uy @ y=0
%         figure(3); plot(uy_plot_s(in,:,1),uy_plot_s(in,:,2),'-o','linewidth',1); hold on
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('uy field vs x; y=0'); xlabel('x'); ylabel('uy');
%         legend('1','2','3','4','5','6','7');
%         f = gcf; exportgraphics(f,'uy_vs_x_y=0.png','Resolution',300); 
%         disp('Plot: uy_vs_x');
%         
%         % ux @ x=0
%         figure(4); plot(ux_plot_s(in,:,1),ux_plot_s(in,:,2),'-o','linewidth',1); hold on
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('ux field vs y; x=0'); xlabel('y'); ylabel('ux');
%         legend('1','2','3','4','5','6','7');
%         f = gcf; exportgraphics(f,'ux_vs_y_x=0.png','Resolution',300); 
%         disp('Plot: ux_vs_y');
%         
%         % ux @ y=0
%         figure(5); plot(ux_plot_y0_s(in,:,1),ux_plot_y0_s(in,:,2),'-o','linewidth',1); hold on
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],'fontsize',15,'FontName','Times');
%         title('ux field vs x; y=0'); xlabel('x'); ylabel('ux');
%         legend('1','2','3','4','5','6','7');
%         f = gcf; exportgraphics(f,'ux_vs_x_y=0.png','Resolution',300); 
%         disp('Plot: ux_vs_x');
%     end
%     
%     clear uth_plot ur_plot_th0 uy_plot ux_plot_y0 ux_plot
%     clear ux_plot_y0_s ux_plot_s uy_plot_s uth_plot_s Interf_x Interf_y ur_plot_th0_s
% end

% % Max/Min tracer
% if in_plot2 == 1
%     figure(11),plot(tr_err(:,1),tr_err(:,2),'-o','linewidth',2);
%     set(gca,'linewidth',1,'ticklength',[0.015 0.008],...
%         'fontsize',15,'FontName','Times');
%     title('Max tr vs i');
%     xlabel('i');
%     ylabel('Max tr');
%     f = gcf; exportgraphics(f,'Max_tr_i.png','Resolution',300); 
%     disp('Plot: Max tr');
%     
%     figure(12),plot(tr_err(:,1),tr_err(:,3),'-o','linewidth',2);
%     set(gca,'linewidth',1,'ticklength',[0.015 0.008],...
%         'fontsize',15,'FontName','Times');
%     title('Min tr vs i');
%     xlabel('i');
%     ylabel('Min tr');
%     f = gcf; exportgraphics(f,'Min_tr_i.png','Resolution',300); 
%     disp('Plot: Min tr');
% end
% 
% % Max/Min volume fraction
% if in_plot3 == 1
%     figure(13),plot(al_err(:,1),al_err(:,2),'-o','linewidth',2);
%     set(gca,'linewidth',1,'ticklength',[0.015 0.008],...
%         'fontsize',15,'FontName','Times');
%     title('Max al vs i');
%     xlabel('i');
%     ylabel('Max al');
%     f = gcf; exportgraphics(f,'Max_al_i.png','Resolution',300); 
%     disp('Plot: Max al');
%     
%     figure(14),plot(al_err(:,1),al_err(:,3),'-o','linewidth',2);
%     set(gca,'linewidth',1,'ticklength',[0.015 0.008],...
%         'fontsize',15,'FontName','Times');
%     title('Min al vs i');
%     xlabel('i');
%     ylabel('Min al');
%     f = gcf; exportgraphics(f,'Min_al_i.png','Resolution',300); 
%     disp('Plot: Min al');
% end

    
%     if Rotating
%         for i = 1:in
%             figure(4),plot(ax_plot_s(i,:,1),ax_plot_s(i,:,2),'-o', ...
%                 ax_plot_s(i,:,1),ax_plot_s(i,:,3),'-x', ...
%                 ax_plot_s(i,:,1),ax_plot_s(i,:,4),'-d', ...
%                 'linewidth',1); hold on
%         end
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],...
%             'fontsize',15,'FontName','Times');
%         title('accX field at y=0');
%         xlabel('x');
%         ylabel('accX');
%         legend('ax','cox','cex');
%         f = gcf; exportgraphics(f,'accX_vs_x.png','Resolution',300);
%         disp('Plot: accXY_vs_x');
        
%         for i = 1:in
%             figure(5), plot(ay_plot_s(i,:,1),ay_plot_s(i,:,2),'-o', ...
%                 ay_plot_s(i,:,1),ay_plot_s(i,:,3),'-x', ...
%                 ay_plot_s(i,:,1),ay_plot_s(i,:,4),'-d', ...
%                 'linewidth',1); hold on
%         end
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],...
%             'fontsize',15,'FontName','Times');
%         title('accy field at y=0');
%         xlabel('x');
%         ylabel('accY');
%         legend('ay','coy','cey');
%         f = gcf; exportgraphics(f,'accy_vs_x.png','Resolution',300);
%         disp('Plot: accY_vs_x');
        
%         for i = 1:in
%             figure(6),plot(ax_plot_s(i,:,1),ax_plot_s(i,:,5),'-o', ...
%                 ax_plot_s(i,:,1),ax_plot_s(i,:,6),'-x', ...
%                 ux_plot_y0_s(i,:,1),ux_plot_y0_s(i,:,2),'-d', ...
%                 'linewidth',1); hold on
%         end
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],...
%             'fontsize',15,'FontName','Times');
%         title('ux field at y=0');
%         xlabel('x');
%         ylabel('ux');
%         legend('ux[]','ux[-1]','ux[]out');
%         
%         for i = 1:in
%             figure(7),plot(ay_plot_s(i,:,1),ay_plot_s(i,:,5),'-o', ...
%                 ay_plot_s(i,:,1),ay_plot_s(i,:,6),'-x', ...
%                 uy_plot_s(i,:,1),uy_plot_s(i,:,2),'-d', ...
%                 'linewidth',1); hold on
%         end
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],...
%             'fontsize',15,'FontName','Times');
%         title('uy field at y=0');
%         xlabel('x');
%         ylabel('uy');
%         legend('uy[]','uy[-1]','uy[]out');
        
%         for i = 1:in
%             figure(7),plot(ax_plot_s(i,:,1),ax_plot_s(i,:,5),'-o', ...
%                 ax_plot_s(i,:,1),ax_plot_s(i,:,6),'-x', ...
%                 ux_plot_y0_s(i,:,1),ux_plot_y0_s(i,:,2),'-d', ...
%                 'linewidth',1); hold on
%         end
%         set(gca,'linewidth',1,'ticklength',[0.015 0.008],...
%             'fontsize',15,'FontName','Times');
%         title('ux field at y=0');
%         xlabel('x');
%         ylabel('ux');
%         legend('ux[]','ux[-1]','ux[]out');
% 
%     end






