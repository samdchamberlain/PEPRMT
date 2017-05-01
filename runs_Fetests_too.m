%%Looper for model runs
clear all

%files and parameter settings
siteFiles = {'WP_xdata.csv'};
OC_init = [0 3e9];% (3e9)/2 3e9]; %intital C stocks to use
Fe_pools = [0 1]; %simulate both presence

%Parameter changes; for now use only one set
                %alpha %Ea          %Km_2
% params(:,:,1) = [9e10 ((65.2)*1000) 1.8e-5; %original values
%                 9e11 ((68)*1000)   1.3e-6];
params(:,:,1) = [5e10 ((66)*1000)   1.8e-5; %balanced intermediate
                 6e10 ((66)*1000)   1.8e-5];             

%init for looping program             
nDataFiles = size(siteFiles,1);
OC_len = size(OC_init,2);
nParam = size(params,3);
nFe = length(Fe_pools);
count = 0;

%loop for model runs
for fi = 1:nDataFiles
    count = count+1;
    
    filename = siteFiles{fi};
    importdata(filename, ',',1);
    xdata = ans.data;
    
    %To use input GPP, first set positive values to zero
    ind = xdata(:,9) > 0;
    xdata(ind,9) = 0; % +GEP set to zero
    
    if strmatch('EE',filename) == 1
        xdata(:,4) =  xdata(:,4)+12; %make sure East End doesn't go dry
        site = 'EE';
    elseif strmatch('MB',filename) == 1
        xdata = xdata;
        site = 'MB';
    else
        xdata = xdata;
        site = 'WP';
    end
    
    %run GPP model and compare to known fluxes
    %GPP = PEPRMT_sys_CO2_GPP(xdata);
    GPP = xdata(:,9)'*60*60*24; % or use eddy GPP (umol m-2 d-1)
    
    %set output files
    M.CH4 = NaN(OC_len, length(GPP), nFe);
    M.S1 = NaN(OC_len, length(GPP), nFe);
    M.S2 = NaN(OC_len, length(GPP), nFe);
    M.M1 = NaN(OC_len, length(GPP), nFe);
    M.M2 = NaN(OC_len, length(GPP), nFe);
    M.GPP = NaN(OC_len, length(GPP), nFe);
    
    for i = 1:OC_len
        count = count+1;
        
        for j = 1:length(Fe_pools)
            
        %param_set = params(:,:,j); %load parameter set
        Fe = Fe_pools(j);
        
        %run Reco and soil pool model
        SOM = OC_init(:,i);
        [NEE_mod, S1, S1sol, S2, Reco_1, priming_re] = PEPRMT_DAMM_sys_Reco_prime(xdata, SOM, GPP);
        
        %C model info
        run{i} = strcat('C_',num2str(SOM)); 
      
            
            %run FCH4 model
            [FCH4, M1, M2, M_Vmax1, M_Vmax2, priming_coef, WT_2_adj, CH4_water] = PEPRMT_TP_sys_CH4_SDC(xdata, S1, S2, GPP, params, Fe);           
        
            %save output
            M.CH4(i,:,j) = FCH4;
            M.S1(i,:,j) = S1sol;
            M.S2(i,:,j) = S2;
            M.M1(i,:,j) = M1;
            M.M2(i,:,j) = M2;
            M.NEE(i,:,j) = NEE_mod;
            M.ER(i,:,j) = Reco_1;
            M.GPP(i,:,j) = GPP;
        end

    end

end

%% plotting section
% figure(1);
% plot(xdata(:,11),'.'); hold on;
% plot(M.NEE(1,:,1)); hold on; plot(M.NEE(2,:,1)); hold on;
% plot(M.NEE(1,:,2)); hold on; plot(M.NEE(2,:,2)); 
% legend('true flux', '0 kgC: no Fe', '36 kgC: no Fe', '0 kgC: Fe', '36kgC: Fe')
% 
% figure(2);
% plot(xdata(:,10),'.'); hold on;
% plot(M.ER(1,:,1)); hold on; plot(M.ER(2,:,1)); 
% legend('true flux', '0 kg C', '36 kg C')

figure(3);
plot(xdata(:,9),'.'); hold on;
plot(M.GPP(1,:,1)*(12.01/10^6)); hold on; plot(M.GPP(2,:,1)*(12.01/10^6)); hold on;
plot(M.GPP(1,:,2)*(12.01/10^6)); hold on; plot(M.GPP(2,:,2)*(12.01/10^6)); 
legend('true flux', '0 kgC: no Fe', '36 kgC: no Fe', '0 kgC: Fe', '36kgC: Fe')
% 
% figure(4);
% plot(xdata(:,12)*(60*60*24/(1e6)), '.'); hold on; 
% plot(M.CH4(1,:,1)./1000); hold on; plot(M.CH4(2,:,1)./1000); hold on;
% plot(M.CH4(1,:,2)./1000); hold on; plot(M.CH4(2,:,2)./1000); 
% %legend('36 kg C m-3', '18 kg C m-3', '0 kg C m-3');
% legend('true flux', '0 kgC: no Fe', '36 kgC: no Fe', '0 kgC: Fe', '36kgC: Fe')


%% Finalized plots for Mayberry
time = xdata(:,1)+4015;

c_conv = 1.037664; %conversion from nmol or umol m-2 s-1 to mgC-CH4 m-2 d-1 or g C-CO2 m-2 d-1, respectively

figure
subplot(2,1,1)
plot(time, xdata(:,12).*c_conv, '.'); hold on; %conversion to mgC m-2 d-1
plot(time, M.CH4(1,:,1).*(12.01/1000)); hold on; plot(time, M.CH4(2,:,1).*(12.01/1000)); hold on;
plot(time, M.CH4(2,:,2).*(12.01/1000)); 
legend('observed flux', '0 kg C: no Fe', '36 kg C: no Fe', '36 kg C: Fe');
ylabel('F_{CH4} (g CH_{4}-C m^{-2} d^{-1})')
datetick('x',2,'keeplimits')

subplot(2,1,2)
plot(time, xdata(:,10)*c_conv,'.'); hold on;
plot(time, M.ER(1,:,1)*(12.01/10^6)); hold on; plot(time, M.ER(2,:,1)*(12.01/10^6)); hold on;
plot(time, M.ER(2,:,2)*(12.01/10^6));
%legend('observed flux', '0 kg C', '36 kg C')
ylabel('ER (g CO_{2}-C m^{-2} d^{-1})')
datetick('x',2,'keeplimits')
%% Regressions CH4
mdl = fitlm(M.CH4(2,:,1).*(1000/(24*60*60)), xdata(:,12))

figure
scatter(M.CH4(2,:,1).*(1000/(24*60*60)), xdata(:,12), '.'); ylabel('Obs F_{CH4}'); xlabel('Model F_{CH4}'); refline([1 0]);
