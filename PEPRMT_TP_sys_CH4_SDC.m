%DAMM CH4 model 
%Based on Davidson 2012 DAMM model
%parameterized for restored wetlands in the Sacramento-San Joaquin River
%Delta
%Created by Patty Oikawa, modified by Sam Chamberlain (2017)
%patty.oikawa@gmail.com, schamberlain@berkeley.edu

%all PEPRMT models use the same input structure (xdata) for CH4 models
%however not all models use all variables in the structure
%all variables are at the daily time step

function  [pulse_emission_total, M1, M2, M_Vmax1, M_Vmax2, priming_coef, WT_2_adj, CH4water] = PEPRMT_final_sys_CH4_alt(xdata, S1, S2, GPP, param_set, Fe_pools)
%This is run at daily time step--all variables in daily time step

%Constants
R = 8.314;                  %J K-1 mol-1
Time_2 =xdata(:,1);         %day of year
DOY_disc_2=xdata(:,2);      %day of year that starts over every year

%Exogenous Variables
TA_2 = xdata(:,3);          %Air temperature- measured (C)
WT_2 = xdata(:,4);          %Water table height (m) equals 0 when water table at soil surface
GPP_2 = abs(GPP);           %Modeled or real GPP - use output from PEPRMT-GPP (umol m-2 d-1)
S1_2=S1;                    %Modeled SOC pool - use output from PEPRMT Reco model -cumulative umol m-2
S2_2=S2;                    %Modeled labile C pool - use output from PEPRMT Reco model -cumulative umol m-2
wetland_age_2=xdata(:,7);   %age of wetland in years (whole numbers only)

WT_2_adj=(WT_2/100)+1;      %makes a new variable where wt=1 at soil surface

%CH4 PARAMETERS
%SOC pool
M_alpha1 = param_set(1,1);  %umol m-3 s-1
M_ea1 = param_set(1,2);     %parameter in kJ mol-1; multiplied by 1000 = J mol-1
M_km1 = param_set(1,3);     %umol m-3 
%Labile C pool
M_alpha2 = param_set(2,1);
M_ea2 = param_set(2,2);
M_km2 = param_set(2,3);

%CH4 oxidation parameters
M_alpha3 = 9e10;
M_ea3 = (59.6)*1000;
M_km3 =253e-5;

%parameter for hydrodynamic flux
k=0.04;                     %gas transfer velocity(m d-1)

%parameters for plant-mediated transport
Vtrans=0.24;                %gas transfer velocity through plants(m d-1)
Oxi_factor=0.35;            %percent oxidized during transport

%empirical factors for inhibition of ch4 production when WT drops
beta1=0.48;
beta2=-0.18;
beta3=0.0042;

%empirical factor for decaying inhibition of CH4 production following first
%flooding of wetland
zeta1=5.1e-6;
zeta2=0.00058;
zeta3=0.11;
  
GPPmax=max(GPP_2);          %plant-mediated transport/priming parameter

%Time Invariant
RT = R .* (TA_2 + 274.15);  %T in Kelvin - all units cancel out
M_Vmax1 = M_alpha1 .* exp(-M_ea1./RT); %umol m-2 s-1 
M_Vmax2 = M_alpha2 .* exp(-M_ea2./RT);
M_Vmax3 = M_alpha3 .* exp(-M_ea3./RT);

%priming coefficient (SDC addition)
priming_coef = (2*(GPP_2/GPPmax)); %priming scales with GPP
priming_coef = priming_coef - min(priming_coef) + 1; %set min to 1 (no priming)

%preallocating space
S1sol = zeros(1,length(Time_2));
S2sol = zeros(1,length(Time_2));

M1 = zeros(1,length(Time_2));
M2 = zeros(1,length(Time_2));
M_full=zeros(1,length(Time_2));
M_percent_reduction=zeros(1,length(Time_2));
M_percent_reduction_2=zeros(1,length(Time_2));
M_percent_reduction_3=zeros(1,length(Time_2));

CH4water=zeros(1,length(Time_2));
Hydro_flux=zeros(1,length(Time_2));
Plant_flux=zeros(1,length(Time_2));
Plant_flux_net=zeros(1,length(Time_2));
CH4water_store=zeros(1,length(Time_2));
CH4water_0=zeros(1,length(Time_2));
Oxi_full=zeros(1,length(Time_2));
R_Oxi=zeros(1,length(Time_2));
CH4water_0_2=zeros(1,length(Time_2));
trans2=zeros(1,length(Time_2));

for t = 1:length(Time_2)
    
    %parameter for plant-mediated transport--function of GPP
    trans2(t)= GPP_2(t)/GPPmax;

    if trans2(t)<0 %set between 0 and 1
        trans2(t)=0;
    end
    if trans2(t)>1
        trans2(t)=1;
    end


    %following Davidson and using multiple eq for different substrate pools
    M1(t) = M_Vmax1(t).*S1_2(t)./(M_km1+S1_2(t)); %umol m2 sec rxn velocity
   
    if S2_2(t)==0   %in winter, no CH4 from plant C
        M2(t)=0;
    else
        M2(t) = M_Vmax2(t).*S2_2(t)./(M_km2+S2_2(t)) ; %umol m2 sec
    end
    
    if M1(t)<0 %make sure CH4 pools cant go negative
        M1(t)=0;
    end
    
    if M2(t)<0
        M2(t)=0;
    end
    
%Empirical eq Oikawa for CH4 inhibition when WT falls below soil
%surface--if WT below soil surface any time in 10 days previous--CH4
%production reduced

   if t<=20 && WT_2_adj(t)<1
        M_percent_reduction(t)=(beta1*WT_2_adj(t).^2)+(beta2*WT_2_adj(t))+beta3;
   else
        M_percent_reduction(t)=1;
   end
    
   Sel=nan(20,1);
   
   if t>20
        Sel=WT_2_adj(21-19:21);
   end
   
   if nanmin(Sel)<1
        M_percent_reduction(t)=(beta1*WT_2_adj(t).^2)+(beta2*WT_2_adj(t))+beta3;
   else
        M_percent_reduction(t)=1;
   end
   
    
   if WT_2_adj(t)<0
        M_percent_reduction(t)=0;
   end
    
   if M_percent_reduction(t)<0
        M_percent_reduction(t)=0;
   end
   
   if M_percent_reduction(t)>1
        M_percent_reduction(t)=1;
   end

%Empirical eq Oikawa for CH4 inhibition following restoration
   if wetland_age_2(t)<2
        M_percent_reduction_2(t)=(zeta1*DOY_disc_2(t).^2)+(zeta2*DOY_disc_2(t))+zeta3;
   else
        M_percent_reduction_2(t)=1;
   end
   
   if  M_percent_reduction_2(t)>1
        M_percent_reduction_2(t)=1;
   end

%Simulate Fe present in the soil; only impacts M2 as oxygen is input to
%root zone only
    if Fe_pools == 1 
        M_percent_reduction_3(t)=0; %reduced methane flux to 0% (likely extreme)
    else
        M_percent_reduction_3(t)=1; %no affect w/o Fe presence
    end
     
    M1(t) = M1(t)*M_percent_reduction(t) ;  %umol m2 sec
    M2(t) = M2(t)*M_percent_reduction(t) ;  %umol m2 sec
    M1(t) = M1(t)*M_percent_reduction_2(t); %umol m2 sec
    M2(t) = M2(t)*M_percent_reduction_2(t); %umol m2 sec
    M2(t) = M2(t)*M_percent_reduction_3(t); %umol m2 sec
   
    %incorporate a priming effect to CH4 production
    M1(t) = M1(t)*priming_coef(t);
    M2(t) = M2(t)*priming_coef(t);

    %S1sol and S2sol are the new SOC and labile pools adjusted for C lost
    %thru CH4
        S1sol(t) = S1_2(t) - (M1(t)*60*60*24); %accounts for depletion of soil C due to Reco and FCH4
        S2sol(t) = S2_2(t) - (M2(t)*60*60*24);
    
    if S1sol(t)<0   %make sure values don't go below zero
        S1sol(t)=0;
    end
    
    if S2sol(t)<0
        S2sol(t)=0;
    end
   
    M_full(t)=(M1(t)*60*60*24)+(M2(t)*60*60*24); %total CH4 produced at this time step in umol m-3 soil day-1
    
end

%NOW COMPUTE CH4 TRANSPORT %%%%%%%%%%%%%%
for t = 1:length(Time_2)

    %make sure WT_2_adj is never negative
    if WT_2_adj(t)<=0
        WT_2_adj(t)=0;
    end
    
    %now start ch4 transport loop
    if t==1
        if WT_2_adj(t)>1 %WT_2_adj = 1 at soil surface or above
            
            %Methane oxidation is zero when WT above surface
            R_Oxi(t) = 0;
            Oxi_full(t)=0; %umol m-3 d-1
      
            %This assumes you start out the year with no CH4 in water
            %where CH4water_0 is the initial concentration of CH4 in water
            CH4water_0(t)=0;    %umol m-3
            CH4water_0_2(t)=0;  %CH4 produced in previous time step
      
            %Only model CH4 in 1m3 of water 
            CH4water(t)= M_full(t)+CH4water_0(t); %umol per m^3 - concentration in water doesn't change

            %based on the concentrations in the soil and water, you get hydro and
            %plant-mediated fluxes
            Hydro_flux(t)=k*CH4water(t);                %umol m-2 d-1; Hydrodynamic flux Poindexter
            Plant_flux(t)=Vtrans*CH4water(t)*trans2(t); %umol m-2 d-1; Plant mediated transport Tian 2001--which just uses CH4 in soil
            Plant_flux_net(t)=Plant_flux(t)*Oxi_factor; %umol m-2 d-1; Plant mediated transport after oxidation

            %subtract moles of methane lost from the pools (soil and
            %water) to the atm
            CH4water_store(t)=CH4water(t)-Hydro_flux(t)-Plant_flux(t); %umol m-3 stored in the system

        else
            %If you start out the year with no water above the surface
            %gives CH4 concentration in water which is now less than 1 m3
            CH4water_0(t) = M_full(t); %umol per m^3 - concentration in soil and water are the same

            %Methane oxidation turns on when WT falls below the surface
            R_Oxi(t) = M_Vmax3(t) .* CH4water_0(t) ./ (M_km3 + CH4water_0(t)) ; %umol m2 sec Reaction velocity
            Oxi_full(t) = R_Oxi(t)*60*60*24;%umol m-2 d-1    
            CH4water(t) = CH4water_0(t) - Oxi_full(t); %now you have less ch4
      
            if CH4water(t)<0    %never let CH4 conc go negative   
                CH4water(t)=0;
            end
            
            CH4water_0_2(t)=0;
            
            Hydro_flux(t)=k*CH4water(t);  %umol m-2 d-1; Hydrodynamic flux Poindexter
            Plant_flux(t)=(Vtrans*CH4water(t))*trans2(t);  %umol m-2 d-1; Plant mediated transport Tian 2001--which just uses CH4 in soil
            Plant_flux_net(t)=Plant_flux(t)*Oxi_factor;  %umol m-2 d-1; Plant mediated transport after oxidation

            CH4water_store(t)=CH4water(t)-Plant_flux(t)-Hydro_flux(t); %umol m-3 stored in the system       
        end
        
    else
        if WT_2_adj(t)>1  %if you have water, then the CH4 should mix between the 2 layers and concentrations should be the same in water and soil
      
            %current CH4 concentration
            CH4water_0(t)= CH4water_store(t-1);
      
            %Methane oxidation is zero when WT above surface
            R_Oxi(t) = 0;
            Oxi_full(t)=0;%umol m-3 d-1
            CH4water_0_2(t)=0;
      
            %Now add the new CH4 produced today to the soil
            CH4water(t)= M_full(t)+CH4water_0(t); %umol m-3
      
            %again compute fluxes
            Hydro_flux(t)=k*CH4water(t);  %umol m-2 d-1; Hydrodynamic flux Poindexter      
            Plant_flux(t)=(Vtrans*CH4water(t))*trans2(t);  %umol m-2 d-1; Plant mediated transport Tian 2001--which just uses CH4 in soil    
            Plant_flux_net(t)=Plant_flux(t)*Oxi_factor;  %umol m-2 d-1; Plant mediated transport after oxidation

            %subtract the moles of methane lost from the pools
            CH4water_store(t)=CH4water(t)-Plant_flux(t)-Hydro_flux(t); %umol m-3 stored in the system
            
        else  %if you don't have WT_2_adj above the soil, then all the CH4 in the water goes to the soil
            
            CH4water_0(t)=CH4water_store(t-1);

            %now add new CH4 to this new concentrated pool of CH4 in the soil
            CH4water_0_2(t)=M_full(t)+CH4water_0(t);   %umol per m-3

            %Methane oxidation turns on when WT falls below the surface
            R_Oxi(t) = M_Vmax3(t) .* CH4water_0_2(t) ./(M_km3 + CH4water_0_2(t)) ; %umol m2 sec Reaction velocity     
            Oxi_full(t)=R_Oxi(t)*60*60*24;  %umol m-2 d-1
            CH4water(t)=CH4water_0_2(t)-Oxi_full(t);
      
            if CH4water(t)<0
                CH4water(t)=0;
            end
            
            Hydro_flux(t)=k*CH4water(t);  %umol m-2 d-1; Hydrodynamic flux Poindexter
            Plant_flux(t)=(Vtrans*CH4water(t))*trans2(t);  %umol m-2 d-1; Plant mediated transport Tian 2001--which just uses CH4 in soil
            Plant_flux_net(t)=Plant_flux(t)*Oxi_factor;  %umol m-2 d-1; Plant mediated transport after oxidation
            CH4water_store(t)=CH4water(t)-Plant_flux(t)-Hydro_flux(t); %umol m-3 stored in the system
      
        end    
    end  
end

pulse_emission_total=Plant_flux_net+Hydro_flux;%umol CH4 m-2 day-1; total CH4 flux to atm
