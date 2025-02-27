
function [NPP, Rd, T_lf, InsPP, InsRd, InsEP] = RP2_forFlux(DOY, Lat, lai, netRAD, Tmax, Tmin, Ta, WND, CO2a,P, C3C4, Leaf_width)

% Descrption:
%
% The radiation is used in SWAT, EPIC, and DCaPS
% This version remove the turbulance r_t, which is influenced by height and
% wind
% Date: 06-Jan-2020 15:57:30
% Author: Tongxi Hu
% Username: hu.1555

% ------------------------------ constants --------------------------------
sc = 1367; %4.896; %  4.921 MJ/m^2/hr  = 1367 W/m^2
pi = 3.1415926;
omega = pi/12;  % rad/hr

% ----------------------------- environmental -----------------------------
%DOY = 298;      % day of year
%Lat = -27.5;
Latr = Lat/180*pi;

% C3C4 = 0;
C3 = C3C4;  C4 = 1- C3C4;

%Tmax = 30;
%Tmin = 20;
%WND = 1.5; % wind speed
Atmo_Tran = 0.75;

% replace missing values with default values
CO2a(CO2a < -9000) = 400; % umol / mol   ambient air CO2 partial pressure     mol/m3
WND(WND < -9000) = 1.5;
P(P<-9000) = 101.325;
% ------------------------------  crop related ----------------------------
% lai = 6;
% Accroding to CLM5.0  Table 3.1
% for C3 Crop chil = -0.3
% for Temp corn chil = -0.3
% for Temp soybean chil = -0.5
chil = -0.5;
ref(1) = 0.11;  % leaf reflectance for PAR
ref(2) = 0.35 ;   % leaf reflectance for NIR
trans(1)= 0.05 ; % leaf transmittance for PAR
trans(2) = 0.34; % leaf tansmittance for NIR
albg = 0.3;   % soil albedo
%Height = 1.5 ;% m
%Leaf_width = 0.1 ; % m

% -------------------------% parameters for cal inter CO2 precessure ----
% parameters for CSD method
O2 = 210   ;   %(mmol/mol)
CiCa_ratio = 0.12 * C3 + 0.19 * C4;

% ----------------------- photosynthesis rate -----------------------------
% po2m = 20900; % Pa   0.2*Air Pressure was used in CLM  DaCPS: 209000 mubar

% ==========================================================
% rstfac = 1.0;  % water stress
Kd = 0.719;   % diffuse and scattered diffuse PAR extinction coefficient
Vcmax25 = 100.7; % mumol/m2/s  Ref CLM4.5/5.0
% effcon = 0.08; %  according to USGS  initial value
Kn = 0.11;      % extinction coeff in N


%% ******************* I. Prepare hourly environmental factors *************

% ====================== (1). SWAT radaition calculation ====================

%E0 = (r0/r)^2;  % r0 : the mean earth-sun distance
% r: the earh-sun distance for any given day of a year
E0 = 1+0.033*cos(2*pi*DOY/365);

SolarDec = 23.45*sin(pi*2*(DOY+284)/365)*pi/180;   % in radians

T_SR = acos(-tan(Latr)*tan(SolarDec))/omega;  % duriation from noon to the hour of sunrise
T_SS = -acos(-tan(Latr)*tan(SolarDec))/omega; % duriation from noon to the hour of sunset


h_SR = 12 - T_SR;  % sunrise hour
h_SS = 12 - T_SS;  % sunset hour

DayLength = 2*acos(-tan(Latr)*tan(SolarDec))/omega;

% The maximum Radiation reaches the earth on a clear day
totR = 24*E0*sc*3600*(omega*T_SR*sin(SolarDec)*sin(Latr)+ cos(SolarDec)*cos(Latr)*sin(omega*T_SR))/pi;
% J/m2/day
PreDawn =floor(12 - T_SR);
PostDusk = ceil(12 - T_SS);

Tot_Sh = DayLength*sin(Latr)*sin(SolarDec) + 24*cos(Latr)*cos(SolarDec)*sin(pi*DayLength/24)/pi;
AccSolar  = zeros(PostDusk - PreDawn +1, 0);
AccFrac = zeros(PostDusk - PreDawn +1, 0);
for t = 12 : PostDusk
    
    if t < h_SS
        
        AccSolar(t) = (t-12)*sin(Latr)*sin(SolarDec) + 12*cos(Latr)*cos(SolarDec)*sin(pi*(t - 12)/12)/pi;
        AccFrac(t)= AccSolar(t)/Tot_Sh;
    elseif t>= h_SS
        AccSolar(t) = (h_SS-12)*sin(Latr)*sin(SolarDec) + 12*cos(Latr)*cos(SolarDec)*sin(pi*(h_SS - 12)/12)/pi;
        AccFrac(t)= AccSolar(t)/Tot_Sh;
    end
    
end
% the fraction of radiation in each hr after predown and before postdusk
Frac = zeros(24,1);
Frac(13:PostDusk) = AccFrac(13:PostDusk)- AccFrac(12: PostDusk-1);
Frac(1:12) = flip(Frac(13:24));

% Diffuse light fraction (FRDF) from atmospheric transmission (Atmo_Tranmiss)
if Atmo_Tran < 0.22
    
    Diff_frac = 1;
elseif Atmo_Tran > 0.22 && Atmo_Tran < 0.35
    Diff_frac = 1.-6.4*(Atmo_Tran-0.22)^2 ;
else
    Diff_frac = 1.47-1.66*Atmo_Tran;
end

totRAD = totR*Frac*Atmo_Tran/3600;
idR = find(netRAD < -9900);
if sum(idR) > 0
    netRAD = totRAD;
end
% total incoming PAR and NIR reaching the earch during

PAR = 0.5*netRAD;    % J/m2/s
NIR = 0.5*netRAD;

% ============= (2). de Wit hourly temperature simulation ================

simuTa = zeros(24,1);   % in degreee Cesiu

for t = 1 : 24
    Tss = Tmin + (Tmax - Tmin)*sin(pi*(h_SS - 12 + DayLength/2)/(DayLength + 3));
    if t >= h_SR && t<= h_SS
        simuTa(t) = Tmin + (Tmax - Tmin)*sin(pi*(t - 12 + DayLength/2)/(DayLength + 3));
        
    elseif t > h_SS
        
        simuTa(t) = (Tmin - Tss*exp(-(24-DayLength)/4) + (Tss - Tmin) *exp(-(t-h_SS)/4))/...
            (1- exp(-(24- DayLength)/4));
    elseif t < h_SR
        simuTa(t) = (Tmin - Tss*exp(-(24-DayLength)/4)+ (Tss-Tmin)*exp(-(t+h_SR)/4))/...
            (1- exp(-(24- DayLength)/4));
    end
    % SWAT method
    % T_air(t) = (Tmax + Tmin)/2 + (Tmax - Tmin)/2*cos(0.2618*(t-15));
end

% replace missing values with simulated Ta
idT = find(Ta < -9900);
Ta(idT) = simuTa(idT);
%save('T_airSWAT.mat','T_air')
% =============== (3) vapor pressure ===========================
% SVP = 0.611 * exp(17.4 * T_air ./ (239 + T_air));     % kPa
obsTmin = min(Ta);
AVP = 0.611 * exp(17.4 * obsTmin / (239 + obsTmin));

%% **************** II. Biomass Accumulation by Hour **********************
InsPP = zeros(24, 1);
InsRd = zeros(24, 1);
InsEP = zeros(24, 1);
T_lf = zeros(24, 1)-9999.0; 

for t = 1:24
    
    %% ================ (1). Radiation absorbed by Canopy =====================
    if t >= h_SR && t <= h_SS
        
        CosTheta = sin(Latr)*sin(SolarDec)+ cos(Latr)*cos(SolarDec)*cos(pi*(t - 12)/12);
        
        fDiff = max(Diff_frac,0.15+0.85*(1.- exp(-0.1/CosTheta)));
        
        % ---incoming diffuse PAR (PARDF) and direct PAR (PARDR)
        PAR_diff = PAR(t) .* fDiff ;               % J/m2/s
        PAR_dir = PAR(t) - PAR_diff;
        %*---incoming diffuse NIR (NIRDF) and direct NIR (NIRDR)
        NIR_diff = NIR(t) .* fDiff ;
        NIR_dir = NIR(t) - NIR_diff;
        
        [Idir_sun, ~, Idif_sun, Idif_sha, Kb] = Rad_TwoStream_Dai(lai, chil, albg, ref, trans, CosTheta);
        
        absPAR_Sun = PAR_dir * Idir_sun(1) + PAR_diff* Idif_sun(1);
        absPAR_Sha = PAR_diff* Idif_sha(1);
        
        absRad_Sun = absPAR_Sun + NIR_dir * Idir_sun(2) + NIR_diff * Idif_sun(2);
        absRad_Sha = absPAR_Sha + NIR_diff* Idif_sha(2);
        
        %absRad_tot = absRad_Sun + absRad_Sha;
        
        %% =================== (2). CO2, H2O and Heat ==============================
        Kw = 0.5491  ;%  according to CSD , leaf angle 60,  LAI = 6, scatter = 0.2
        %*---fraction of sunlit and shaded components in canopy
        Sunlit = 1./Kb/lai*(1.-exp(-Kb*lai));
        Shaded = 1.-Sunlit;
        
        % scaling-up coefficients from leaf to canopy--2014
        % cintsha(1) = (1.-exp(-extkn*lai))/extkn - cintsun(1)
        % Eq. 37a, 37b, 38a, 38b
        cintsun(1) = (1.-exp(-(Kn + Kb).*lai))/(Kn+ Kb); % 0.11 was used by Bonan
        cintsun(2) = (1.-exp(-(Kb+Kd).*lai))/(Kb+Kd);
        cintsun(3) = (1.-exp(-Kb.*lai))/Kb;
        
        cintsha(1) = (1.-exp(-Kn.*lai))/Kn - cintsun(1);
        cintsha(2) = (1.-exp(-Kd.*lai))/Kd - cintsun(2);
        cintsha(3) = lai - cintsun(3);
        
        %*---turbulence resistance for canopy (r_t) and for soil (RTS)
        % r_t = 0.74*(log((2.-0.7*Height)/(0.1*Height)))^2/(0.4^2*WND);
        
        %*---boundary layer resistance for canopy, sunlit and shaded leaves
        rbh = 100*sqrt(Leaf_width/WND(t)); %Leaf boundary layer conductance for heat transfer  m/s
        gb_canopy = (1.-exp(- 0.5*Kw *lai))/(0.5*Kw )/rbh; % Canopy boundary layer conductance for heat transfer
        % Kw: wind extinction coeff in the canopy  ?????????????????
        
        gb_Sun = (1.-exp(-(0.5*Kw+Kb)*lai))/(0.5*Kw+Kb)/rbh;
        gb_Sha = gb_canopy - gb_Sun;
        
        rbh_Sun = 1./gb_Sun;    % !boundary layer resistance to heat,sunlit part
        rbw_Sun = 0.93*rbh_Sun; % boundary layer resistance to H2O, sunlit part
        rbh_Sha = 1./gb_Sha;    % !boundary layer resistance to heat,shaded part
        rbw_Sha = 0.93*rbh_Sha; % !boundary layer resistance to H2O, shaded part
        
        % r_t_Sun = r_t*Sunlit;
        % r_t_Sha = r_t*Shaded;
        
        %% ==================== leaf N contents====================================
        % N is not considered becasue it was included in Vmax25
        
        %% ================(3) Instaneious Photosynthesis =========================
        % cal photosynthesis and respiration for sunlit and shaded leaves
        % !!! Estimation of leaf temperature is coupled within A-gs photosynthesis model
        % However, the CSD method require vapor pressure input, the CLM method
        % requires humidity input.
        % -------------------- 3.1 cal T_lf, CO2i, Ass_n (CSD)-----------------------
        % ref: P13 Eq 2, 6, 7, 8  in CSD book
        
        % --------------------- first round --------------------------------
        % using Ta to approximate T_lf
        T_lf(t) = Ta(t);
        
        %         % test Tf steability
        %         Tf_stable = zeros(4,10);
        %
        %         for it = 1:10
        %
        %          if it == 1
        [SVPL1, CO2i1] = cal_CO2(C3C4, T_lf(t), AVP, CiCa_ratio, CO2a(t), O2);
        
        % % CoLM method
        % qt = (T_lf - 298.16)/10; % qt is Q10
        % kc = 30*2.1^qt;       % Michaelis-Menton constant for carboxylation by Rubisco (Pa)
        % ko = 30000 *1.2^qt;   % Michaelis-Menton constant for oxygenation by Rubisco (Pa)
        % gammas = 0.5*po2m/(2600*0.57^qt) *c3;
        
        % ----------------------------3.2 Photosynthesis Rate ------------------------
        %  Ref Dai. 2004 and CLM
        %  !!!!!! Warning !!!!!: Photo is net photosynthesis here
        [Photo1_Sun, ~] = photo_rate(C3C4,Vcmax25,T_lf(t), cintsun,absPAR_Sun,O2,CO2i1, P(t));
        
        [Photo1_Sha, ~] = photo_rate(C3C4,Vcmax25,T_lf(t), cintsha,absPAR_Sha,O2, CO2i1, P(t));
        
        VPD = max(0., SVPL1 - AVP);
        % -------------------------------- 3.3 Conductance -----------------------
        % cal s (Eq.5) in Eq. (2) P11 CSD book
        % the first round: s is the deviration of saturated VP of T_air
        SlopeL1 = 4158.6 * SVPL1/(T_lf(t) + 239.)^2; % s: slope of saturated vapour pressure curve kPa/degree
        
        % --- cal leaf conductance for CO2 (gc) and the stomatal resistance to water (r_swp)-----
        
        r_s1_Sun = cal_Cdtan(Photo1_Sun, T_lf(t), CO2a(t), CO2i1, rbw_Sun);
        r_s1_Sha = cal_Cdtan(Photo1_Sha, T_lf(t), CO2a(t), CO2i1, rbw_Sha);
        
        %----- cal leaf transpiration using Penman-Monteith equation--------------
        % !*---net absorbed radiation and Ep
        [Ep_Sun1, NetRd_Sun1] = cal_EP(r_s1_Sun, rbw_Sun, rbh_Sun, absRad_Sun,Atmo_Tran,Sunlit, T_lf(t), AVP, SlopeL1, VPD);
        [Ep_Sha1, NetRd_Sha1] = cal_EP(r_s1_Sha, rbw_Sha, rbh_Sha, absRad_Sha,Atmo_Tran,Shaded, T_lf(t), AVP, SlopeL1, VPD);
        
        % ---- cal the difference between leaf temperature and air temperature
        del_t_Sun = cal_deltaT(NetRd_Sun1, Ep_Sun1, rbh_Sun);
        del_t_Sha = cal_deltaT(NetRd_Sha1, Ep_Sha1, rbh_Sha);
        
        T_lf_Sun = Ta(t) + del_t_Sun;
        T_lf_Sha = Ta(t) + del_t_Sha;
        %          Tf_stable(1,it)= del_t_Sun;
        %          Tf_stable(2,it)= del_t_Sha;
        %          Tf_stable(3,it)= T_lf_Sun;
        %          Tf_stable(4,it)= T_lf_Sha;
        T_lf(t)=  (T_lf_Sun + T_lf_Sha)/2;
        % ============ Second round to determine the photosynthesis and T_leaf ====
        % cal internal CO2 pressure first
        %           else
        [SVPL2_Sun, CO2i2_Sun] = cal_CO2(C3C4,T_lf_Sun, AVP, CiCa_ratio, CO2a(t), O2);
        [SVPL2_Sha, CO2i2_Sha] = cal_CO2(C3C4,T_lf_Sha, AVP, CiCa_ratio, CO2a(t), O2);
        
        [Photo2_Sun, Rd2_Sun]  = photo_rate(C3C4,Vcmax25,T_lf_Sun, cintsun,absPAR_Sun,O2,CO2i2_Sun, P(t));
        [Photo2_Sha, Rd2_Sha]  = photo_rate(C3C4,Vcmax25,T_lf_Sha, cintsha,absPAR_Sha,O2,CO2i2_Sha, P(t));
        
        SVP = 0.611 * exp(17.4* Ta(t)/(Ta(t) + 239.0));    % kPa    saturated VPD of air
        
        if T_lf_Sun ~= Ta(t) && T_lf_Sha ~= Ta(t)
            Slope_lf2_Sun = (SVPL2_Sun - SVP)/(T_lf_Sun - Ta(t));   % slope of saturated VPD
            Slope_lf2_Sha = (SVPL2_Sha - SVP)/(T_lf_Sha - Ta(t));
        end
        
        % --- cal leaf conductance for CO2 (gc) and the stomatal resistance to water (r_s)-----
        r_s2_Sun = cal_Cdtan(Photo2_Sun, T_lf_Sun, CO2a(t), CO2i2_Sun, rbw_Sun);
        r_s2_Sha = cal_Cdtan(Photo2_Sha, T_lf_Sha, CO2a(t), CO2i2_Sha, rbw_Sha);
        
        %----- cal leaf transpiration using Penman-Monteith equation--------------
        % !*---net absorbed radiation
        [Ep2_Sun, ~] = cal_EP(r_s2_Sun, rbw_Sun, rbh_Sun, absRad_Sun,Atmo_Tran , Sunlit, T_lf_Sun, AVP, Slope_lf2_Sun, VPD);
        [Ep2_Sha, ~] = cal_EP(r_s2_Sha, rbw_Sha, rbh_Sha, absRad_Sha,Atmo_Tran , Shaded, T_lf_Sha, AVP, Slope_lf2_Sha, VPD);
        
        %           % ------------- for testing the stability of Tf ---------------
        %           % ---- cal the difference between leaf temperature and air temperature
        %            del_t_Sun = cal_deltaT(NetRd2_Sun, Ep2_Sun, rb_Sun, r_t_Sun);
        %            del_t_Sha = cal_deltaT(NetRd2_Sha, Ep2_Sha, rb_Sha, r_t_Sha);
        %
        %            T_lf_Sun = T_air(t) + del_t_Sun;
        %            T_lf_Sha = T_air(t) + del_t_Sha;
        %
        %            Tf_stable(1,it)= del_t_Sun;
        %            Tf_stable(2,it)= del_t_Sha;
        %            Tf_stable(3,it)= T_lf_Sun;
        %            Tf_stable(4,it)= T_lf_Sha;
        %          end
        %         end
        %         save([num2str(t),'_clock.mat'], 'Tf_stable');
        % ***************** End of T_lf, An, CO2i *********************************
        
        InsPP(t) = (Photo2_Sun + Photo2_Sha)*3600;   % umol CO2 /m2/hr
        InsRd(t) = (Rd2_Sun + Rd2_Sha)*3600;
        InsEP(t) = (Ep2_Sun + Ep2_Sha)*3600;
      else
        InsPP(t) = 0.00; 
        InsRd(t) = 0.00; 
        T_lf(t) = -9999.0;
    end % during the day
    
end
%Bio = sum(InsPP)*0.409;
NPP = sum(InsPP);%*1E6/44.0/3600;       % g cO2/m2/hr ---> umol CO2/m2/s
Rd = sum(InsRd); %*1E6/44.0/3600;  
end


function [SVP_lf, CO2i] = cal_CO2(C3C4, T_lf, AVP, CiCa_ratio, CO2a, O2)
E_Vcmax = 65330.; % !energy of activation for Vcmx(J/mol)
Rd_Vcmax25 = 0.0089;   % !ratio of dark respiration to Vcmax at 25oC
E_Rd = 46390. ;   % !energy of activation for dark respiration(J/mol)

E_KMC = 79430. ;  % !energy of activation for KMC (J/mol)
E_KMO = 36380. ;  % !energy of activation for KMO (J/mol)
KmC25 =  404.9;   % umol / mol    This only applied to C3 for PP in CLM
KmO25 =  278.4;     % mmol / mol   but different values are used for C4 in CSD

SVP_lf = 0.611 * exp(17.4* T_lf/(T_lf + 239.0));  % kPa

VPD_lf = max(0, SVP_lf - AVP);

KmC = KmC25*exp((1./298.-1./(T_lf+273.))*E_KMC/8.314);
KmO = KmO25*exp((1./298.-1./(T_lf+273.))*E_KMO/8.314);

% compensation point in the absence of dark respiration

Ga_Max = 0.5* exp(-3.3801+5220./(T_lf+273.)/8.314)*O2*KmC/KmO;

% !*---CO2 compensation point (GAMMA)
Rd_Vcmax = Rd_Vcmax25* exp((1./298.-1./(T_lf+273.))*(E_Rd-E_Vcmax)/8.314);

Ga0 = (Ga_Max+Rd_Vcmax*KmC*(1.+O2/KmO))/(1.-Rd_Vcmax); % for C3 crops

if C3C4 == 1
    
    Ga = Ga0 ;
elseif C3C4 ==0
    
    Ga= Ga0/10.0 ;
end
% !*---internal/ambient CO2 ratio, based on data of Morison & Gifford (1983)
CiCa = 1.-(1.-Ga/CO2a)*(0.14+CiCa_ratio*VPD_lf);
% !*---intercellular CO2 concentration
CO2i = CiCa * CO2a;
end


function [Photo, Rd] = photo_rate(C3C4, Vcmax25, T_leaf, cint, PAR, O2, CO2i, P)

% Descrption:
%
%  Modified from Dai and CLM 5.0
%
% Date: 14-Feb-2020 15:03:12
% Author: Tongxi Hu
% Username: hu.1555
s1= 0.3;     %  slope of high temperature inhibition function     (0.3)
s2 = 313.16; % ! 1/2 point of high temperature inhibition function (313.16)
s3 = 0.2;    %! slope of low temperature inhibition function      (0.2)
s4 = 288.15; %! 1/2 point of low temperature inhibition function  (288.16)
s5 = 1.3;    % ! temperature coefficient in gs-a model             (1.3)
s6 = 328.16; % ! temperature coefficient in gs-a model             (328.16)
%trop = 298.16; % ! temperature coefficient in gs-a model             (298.16)
% binter = 0.01; % for C3 plants and 0.04 for C4
P = P*1000;     %101325; % Pa
% effcon = 0.8;
Ea_Vcmax = 65330.; % !energy of activation for Vcmx              (J/mol)
Ed_Vcmax = 149250;
S_Vcmax = 485;

Ea_Jmax = 43540;
Ed_Jmax = 152040;
S_Jmax = 495;

Ea_Rdmax = 46390; % !energy of activation for dark respiration(J/mol)
Ed_Rdmax = 150650;
S_Rdmax = 490;

Ea_Gamma = 37830;

E_KmC = 79430. ;  % !energy of activation for KMC (J/mol)
E_KmO = 36380. ;  % !energy of activation for KMO (J/mol)

KmC25 = 404.9;   % umol /mol
KmO25 = 278.4;   % mmol /mol
Gamma25 = 42.75; % umol /mol

Rgas = 8.314;  % universal gas constant(J/mol/K)

C3 = C3C4;
C4 = 1-C3C4;

T_leaf = T_leaf + 273;

% *************************** CLM4.5/5.0 ************************
% Ref. Chapter 8 in CLM 4.5
Rt = 298.15 * Rgas;
Jmax25 = 1.97*Vcmax25;
Rd25 = (0.015 * C3 + 0.025 * C4) * Vcmax25;

% only for C4 plants
kp25 = 20000*Vcmax25;
qt = (T_leaf - 298.15)/10; % for Q10

% ---------------- 1. Vcmax ---------------------------
fT_V = exp(Ea_Vcmax*(1-298.15/T_leaf)/Rt);
fH_V = (1 + exp((298.15*S_Vcmax - Ed_Vcmax)/Rt))/(1 + exp((T_leaf*S_Vcmax - Ed_Vcmax)/Rt));

fh_V = 1+exp(s1*(T_leaf-s2));
fl_V = 1+exp(s3*(s4-T_leaf));

Vcmax = Vcmax25*fT_V*fH_V * C3 + Vcmax25 * 2.0^qt/fh_V/fl_V*C4;

% scale to the canopy
Vcmax = Vcmax * cint(1);

% Vcmax = Vcmax * WaterStress; % add water stress here

% ----------------2. Jmax ----------------------------
fT_J = exp(Ea_Jmax*(1-298.15/T_leaf)/Rt);
fH_J = (1 + exp((298.15*S_Jmax - Ed_Jmax)/Rt))/(1 + exp((T_leaf*S_Jmax - Ed_Jmax)/Rt));
Jmax = Jmax25*fT_J*fH_J;

Jmax = Jmax * cint(2);

% Jmax = Jmax * WaterStress; % add water stress here
% ----------------3. Rd ------------------------------
fT_R = exp(Ea_Rdmax*(1-298.15/T_leaf)/Rt);
fH_R = (1 + exp((298.15*S_Rdmax - Ed_Rdmax)/Rt))/(1 + exp((T_leaf*S_Rdmax - Ed_Rdmax)/Rt));

Rd = Rd25 * fT_R * fH_R * C3 + Rd25*2.0^qt/(1 + exp(s5*(T_leaf - s6)))*C4;

Rd = Rd * cint(1);
% Rd = Rd * WaterStress;
Rd = Rd * 44*1e-6;
% =============== cal net assimilation rate ====================

% -------------------------- Ac (umol CO2/m2/s)--------------------
% RuBP carboxylase limited rate of carboxylation Ac
Kc = KmC25 * exp(E_KmC*(1-298.15/T_leaf)/Rt);               % umol /mol
Ko = KmO25 * exp(E_KmO*(1-298.15/T_leaf)/Rt);               % mmol /mol
Gamma = Gamma25*exp(Ea_Gamma*(1-298.15/T_leaf)/Rt);         % umol /mol

rrkk = Kc *(1+O2/Ko)*C3;

Ac = Vcmax * ( CO2i-Gamma ) / ( CO2i + rrkk ) * C3 + Vcmax * C4 ;

% ----------------------- Aj, umol CO2/m2/s ----------------------------
% The maximum rate of carboxylation allowed by the capacity to regenerate
% RuBP (the light limited rate )
ePAR = 4.6*PAR;   % PAR  J/m2/s = W/m2  -----> 4.6 umol photon/J

Iphi2 = 0.5*0.85*ePAR;
Phi2 = 0.7;
% the small root of a quadratic function related to Jmax
J = (Iphi2 + Jmax -sqrt((Iphi2+ Jmax).^2 - 4*Iphi2.*Jmax.*Phi2))./2./Phi2;

Aj = J * ( CO2i- Gamma ) / ( 4*CO2i+8.*Gamma) * C3 + ePAR*0.05 * C4;

% ----------------------- Ap, product-limited rate ---------------------
kp = kp25*2.0^qt;  % for C4

% CO2i is the internal leaf CO2 partial pressure
% P is atmosphere pressure (Pa)
Ap = 0.5*Vcmax * C3 + kp*CO2i/P/10 * C4;

% --------------------- co-limitation A --------------------------------
% Ref. Eq . 8.8 in CLM 4.5

Thetaj = 0.98* C3 + 0.8 * C4;
Thetap = 0.95;

sqrtin = max( 0., ( (Ac+Aj).^2 - 4.*Thetaj.*Ac.*Aj ) );

Ai   = (( Ac + Aj ) - sqrt(sqrtin)) ./ ( 2.*Thetaj);

sqrtin= max( 0., ( (Ap+Ai).^2 - 4.*Thetap.*Ap.*Ai));

A = ( ( Ai + Ap ) - sqrt( sqrtin ) ) / ( 2.* Thetap);

Photo = max(0, 44*1e-6*A-Rd)  ;   % net photosynthesis rate
                                      % umol m-2 s-1 --> g CO2 m-2 s-1
% ***************************************************************
end


% ----------------9. conductance of CO2 ---------------------------------
function r_s = cal_Cdtan(Photo, T_, CO2a, CO2i, rb_w)

    gc = Photo*(273.+T_)/0.53717/(CO2a-CO2i);
    
    r_s = max(1E-10, 1./gc - rb_w*1.4)/1.6;    % r_bw has sunlit and shaded leaves

end

% ------------------10. cal leaf transpiration ----------------------------
function [Ep_S, NetRd_Abs] = cal_EP(r_s_S, rbw_S, rb_S, Ab_Rad, Trans,Frac, T_, DVP, Slope, VPD)
    
    Boltz = 5.668E-8;  % Stefan-Boltzmann constant                  J/m2/s/K4
    Latent = 2.4E+6;   % latent heat of water vaporization          J/kg
    V_ca = 1200.0;     % Volumetric heat capacity                   J/m3/degree
    Psy = 0.067;        %psychrometric constant                     kPa/degree
    
    Clear = max(0., min(1., (Trans-0.25)/0.45));    %  !sky clearness
    Blackbody = Boltz*(T_ +273.)^4;             %  !black body radiation
    
    Longwave_net = Blackbody*(0.56-0.079*sqrt(DVP*10.))*(0.1+0.9*Clear)*Frac;
    
    NetRd_Abs = Ab_Rad - Longwave_net;
    
    % !*---intermediate variable related to resistances
    Psr = Psy*(rbw_S +r_s_S)/rb_S;
    
    % !*---radiation-determined term
    Ptr = NetRd_Abs*Slope /(Slope+Psr)/Latent;
    
    % !*---vapour pressure-determined term
    Ptd = (V_ca*VPD/rb_S)/(Slope+Psr)/Latent;  % r_bh : leaf boundary layer resistane to heat (s/m)
    
    %!*---potential evaporation or transpiration
    Ep_S = max(1.E-10,Ptr+Ptd);

end

% ------------------ 11. cal difference T ------------------------------------
function delta_t = cal_deltaT(Net_Rd, EP, rb_s)

    V_ca = 1200.0;     % Volumetric heat capacity             (J/m3/degree)
    Latent = 2.4E+6;   % latent heat of water vaporization          J/kg
    
    delta_t = (Net_Rd -Latent* EP)*rb_s/V_ca;
    
    if delta_t < -25
        delta_t = -25;
    elseif delta_t > 25
        delta_t = 25;
    end

end