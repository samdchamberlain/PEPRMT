%Light Use Efficiency model  to predict GPP
%following Ruimy 1999
%parameterized for restored wetlands in the Sacramento-San Joaquin River
%Delta

%Written by Patty Oikawa
%patty.oikawa@gmail.com

%all PEPRMT models use the same input structure (xdata) for CO2 models
%however not all models use all variables in the structure
%all variables are at the 30min time step


function  [GPP] = PEPRMT_sys_CO2_GPP(xdata)

%Exogenous Variables
Time_2 =xdata(:,1);% day of year
TA_2 = xdata(:,3);%Air temperature- measured (C)
PAR_2 = xdata(:,5);%incoming PAR - measured
LAI_2 = xdata(:,6);%LAI - calculated from GCC following Oikawa et al. 2016

%COMPUTE GPP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%PARAMETERS
k=0.6;
Ha =139.7;
Hd=143;

%CONSTANTS
LUE_mean=0.93;%computed a mean across each growing season (g C m-2 MJ-1)
vcopt = 1.0;
R_t=0.00831;
T_opt = 25 + 274.15;%our Temp opt for Ps is 25C

%EQUATIONS
PAR_2_MJ=(PAR_2*0.0002186)*0.001;   %convert PAR umol m-2 s-1 to MJ m-2 s-1
fPAR_2=0.95*(1-exp(-k*LAI_2));
APAR_2=fPAR_2.*PAR_2_MJ;            %MJ m-2

AirT_K =TA_2+ 274.15;               %C to Kelvin
AirT_K=AirT_K';

max_time = length(TA_2);
vct=zeros(1,length(Time_2));
NPP_FPAR_T=zeros(1,length(Time_2));

for t = 1:max_time
    exponent1=(Ha*(AirT_K(t)-T_opt))./(AirT_K(t)*R_t*T_opt);
    exponent2=(Hd*(AirT_K(t)-T_opt))./(AirT_K(t)*R_t*T_opt);
    top = Hd*exp(exponent1);   
    bottom=Hd-(Ha*(1-exp(exponent2)));
    vct(t) = vcopt*(top./bottom);
    NPP_FPAR_T(t)=((vct(t)*(APAR_2(t)*LUE_mean)));%g C m-2 s-1
end

GPP=(NPP_FPAR_T/12)*10^6*-1*60*60*24;% umol m-2 d-1
    
end














