clc; clear all; close all;

N_pid   = 40;  % # of CPU for parallel computing
N       = 512;
shape   = 1;   % 0: square; 1: rectangular

dt_File = 0.1123; % 0.0244/4; 0.1123; 0.0281
Ini_t   = 429.8844; % 200.0158
Fin_t   = 593-dt_File; % 500-dt_File;

home =  ['/gpfs/scratch/mkim79/basilisk/bio_rec_new3/oxy28/Data_all/'];
home2 = ['/gpfs/scratch/mkim79/basilisk/bio_rec_new3/oxy28/oxy28_save']; mkdir(home2);
home_data  = [home,'/Data_all_'];
home_data2 = [home2,'/Data_all_'];

in = 0;
File_t  = [Ini_t:dt_File:floor(Fin_t/dt_File)*dt_File];
for i = File_t
    i
    in = in+1;
        
    %% Get fields: x, y, velocity, volume fraction, concentration %%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x_field = [];    y_field = [];    ux_field = [];   uy_field = [];
    al_field = [];   tr_field = [];   so_field = [];   ox_field = [];
    
    for j = 0:N_pid-1
        
        fID_field  = fopen([home_data,num2str(N),'_',num2str(i),'_',num2str(j),'.txt'],'r');
        Data_field = fscanf(fID_field,'%f'); fclose(fID_field);
        
        if      shape == 0 % square
            x_tmp  = Data_field(1:6:end);
            y_tmp  = Data_field(2:6:end);
            ux_tmp = Data_field(3:6:end);
            uy_tmp = Data_field(4:6:end);
            al_tmp = Data_field(5:6:end);
            tr_tmp = Data_field(6:6:end);
            
        elseif shape == 1 % rectangular
            x_tmp  = Data_field(1:8:end);
            y_tmp  = Data_field(2:8:end);
            ux_tmp = Data_field(3:8:end);
            uy_tmp = Data_field(4:8:end);
            al_tmp = Data_field(5:8:end);
            tr_tmp = Data_field(6:8:end);
            so_tmp = Data_field(7:8:end);
            ox_tmp = Data_field(8:8:end);
            
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
    
    if     shape == 0
        out = [x_field y_field ux_field uy_field al_field tr_field];
    elseif shape == 1
        out = [x_field y_field ux_field uy_field al_field tr_field so_field ox_field];
    end
    
    %% save files
    writematrix(out,[home_data2,num2str(N),'_',num2str(i),'_tot.txt'],'Delimiter','tab');
end



