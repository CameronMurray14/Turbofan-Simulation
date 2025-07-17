clear
clc
global M_0R height m0R alphaR gamc cpc gamt cpt pi_d pi_b pi_n pi_fn pi_cR pi_cLR pi_fR eta_b eta_mL eta_mH e_cH e_cL e_tL e_tH e_f P_Pstd T_Tstd gc T0R hpr M R_c R_t a_0R V_0R trR eta_r ...
    pi_rR Tt0R Tt1R Tt2R Tt3R Tt4R Tt5R Tt6R Tt7R Tt8R Tt9R tcLR eta_cLR tlambR Tt2_5R pi_cHR tcHR eta_cHR tfR fR Tt13R eta_fR ttHR pi_tHR eta_tHR Tt4_5R ttLR pi_tLR eta_tLR pt9_p9R p0_p9R M9R T9_T0R ...
    V9_a0R V9_V0R pt19_p19R p0_p19R M19R Tt19R T19_T0R V19_a0R V19_V0R pi_d_max 
M_0R = .8; Tt4R = 3000; height = [0,20000,40000]; m0R = 600; alphaR = 8; gamc = 1.4; cpc = .24; gamt = 1.33; cpt = .276; 
pi_d_max = .99; pi_b = .99; pi_n = .99; pi_fn = .99; pi_cR = 36; pi_cLR = 1.7; pi_fR = 1.7; eta_b = .99; eta_mL = .997; ONE = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
eta_mH = .9915; e_cH = .9; e_cL = .89; e_tL = .91; e_tH = .89; e_f = .89; P_Pstd = [1,.4599,.1858]; T_Tstd = [1,.8626,.7519]; Tt4_set = [2600,3200,3537,3714];gc=32.174; theta2 = [.8626*ONE]; theta3 = [.7519*ONE];
T0R = 390;hpr = 18400; M = [0:.9/29:.9]; err =1;
format long
i = 1; j = 1; k = 1;


%¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯

%Preliminary
R_c = 778.16*cpc*(gamc-1)/gamc;
R_t = 778.16*cpt*(gamt-1)/gamt;
a_0R = sqrt(gamc*R_c*gc*T0R);
V_0R = a_0R*M_0R;

%Ram
trR = 1 + (gamc-1) * M_0R^2 / 2;
eta_r = 1;
pi_rR = trR^(gamc/(gamc-1));
Tt0R = T0R*trR;
Tt1R=Tt0R;

%Diffuser
pi_d = pi_d_max * eta_r;

tlambR = (cpt*Tt4R)/(cpc*T0R);

%Compressor
Tt2R = Tt1R;
tcLR = pi_cLR ^ ((gamc-1)/(gamc*e_cL));
eta_cLR = (pi_cLR ^ ((gamc-1)/(gamc)) - 1) / (tcLR-1);
Tt2_5R=T0R*trR*tcLR;
pi_cHR = pi_cR/pi_cLR;
tcHR = pi_cHR ^ ((gamc-1)/(gamc*e_cH));
eta_cHR = (pi_cHR ^ ((gamc-1)/(gamc)) - 1) / (tcHR-1);

%Fan
tfR = pi_fR ^((gamc-1)/(gamc*e_f));
Tt13R = Tt2R*tfR;
eta_fR = (pi_fR ^ ((gamc-1)/(gamc)) - 1) / (tfR-1);

Tt3R = Tt0R*tcLR*tcHR;

%Burner
fR = (cpt*Tt4R-cpc*Tt3R)/(eta_b*hpr-cpt*Tt4R);

%Turbine
ttHR = 1 - (1/(eta_mH*(1+fR)))*(trR/tlambR)*tcLR*(tcHR-1);
pi_tHR = ttHR^(gamt/((gamt-1)*e_tH));
eta_tHR = (1-ttHR)/(1-ttHR^(1/e_tH));
Tt4_5R = Tt4R*ttHR;
ttLR1 = 1/(eta_mL*(1+fR));
ttLR2 = trR/(tlambR*ttHR);
ttLR3 = ((tcLR-1)+alphaR*(tfR-1));
ttLR = 1 - ( 1/(eta_mL*(1+fR)) ) * (trR/(tlambR*ttHR))*(tcLR-1+alphaR*(tfR-1));
    
pi_tLR = ttLR^(gamt/((gamt-1)*e_tL));
eta_tLR = (1-ttLR)/(1-(ttLR^(1/e_tL)));
Tt5R = Tt4_5R*ttLR;
Tt6R = Tt5R;
Tt7R = Tt6R;
Tt8R = Tt7R;

%Exhaust Nozzle
p0_p9_critR = (((gamt+1)/2)^(gamt/(gamt-1)))/(pi_rR*pi_d*pi_cLR*pi_cHR*pi_b*pi_tHR*pi_tLR*pi_n);

if (p0_p9_critR <=1)
 pt9_p9R = (((gamt+1)/2)^(gamt/(gamt-1)));
 p0_p9R = (((gamt+1)/2)^(gamt/(gamt-1)))/(pi_rR*pi_d*pi_cLR*pi_cHR*pi_b*pi_tHR*pi_tLR*pi_n);
 M9R = 1;
else
 pt9_p9R = (pi_rR*pi_d*pi_cLR*pi_cHR*pi_b*pi_tHR*pi_tLR*pi_n);
 p0_p9R = 1;
 M9R = sqrt((2/(gamt-1))*(pt9_p9R^((gamt-1)/gamt) -1));
end
Tt9R = Tt8R;
T9_T0R = (Tt9R/T0R)/(pt9_p9R^((gamt-1)/gamt));
V9_a0R = M9R*sqrt((gamt*R_t)/(gamc*R_c)*T9_T0R);
V9_V0R = V9_a0R*(a_0R/V_0R);

%Fan Nozzle

p0_p19_critR = (((gamc+1)/2)^(gamc/(gamc-1)))/(pi_rR*pi_d*pi_fR*pi_fn);

if(p0_p19_critR <= 1)
   pt19_p19R = ((gamc+1)./2)^(gamc/(gamc-1));
   p0_p19R = (((gamc+1)/2)^(gamc/(gamc-1)))/(pi_rR*pi_d*pi_fR*pi_fn);
   M19R = 1;
else
   pt19_p19R = pi_rR*pi_d*pi_fR*pi_fn;
   p0_p19R = 1;
   M19R = sqrt((2/(gamc-1))*(pt19_p19R^((gamc-1)/gamc) -1));
end
Tt19R = Tt13R;
T19_T0R = (trR*tfR)/(pt19_p19R)^((gamc-1)/gamc);
V19_a0R = M19R*sqrt(T19_T0R);
V19_V0R = V19_a0R*(a_0R/V_0R);

%Performance
ST_CR = (1 ./ (1+alphaR)) * (a_0R/gc) * ( (1+fR)*V9_a0R - M_0R + (1+fR)*(R_t/R_c)*(T9_T0R/V9_a0R)*(1-p0_p9R)/gamc);
ST_BR = (alphaR/(1+alphaR))*(a_0R/gc) * (V19_a0R - M_0R + (T19_T0R/V19_a0R)*((1-p0_p19R)/gamc));
STR = ST_CR + ST_BR;
SR = fR ./ ((1+alphaR)*STR);
eta_thR = (a_0R^2) * ( (1+fR)*(V9_a0R^2) + alphaR * (V19_a0R^2) - (1+alphaR)*M_0R^2    ) / (2*gc*fR * hpr*778.16);
eta_pR = 2 * M_0R * ( (1+fR)*V9_a0R - (1+alphaR)*M_0R + alphaR*V19_a0R ) / ( (1+fR)*V9_a0R^2 + alphaR*V19_a0R^2 - (1+alphaR)*M_0R^2);
eta_OR = eta_thR * eta_pR;
Thrust_R = STR * m0R;

%Alt Loop
while k < 4
T0 = T_Tstd(k) * 518.69;
p0 = P_Pstd(k) * 14.7;

%Throttle Loop
 
while j < 5 
Tt4 = Tt4_set(j);

%Mach Loop

while i < 31
F = OffDesign(M(i),T0,Tt4,p0);
T_CH(i) = F(1);
PICH(i) = F(2);
PI_F(i) = F(3);
M_19(i) = F(4);
M_9(i) = F(5);
Alpha(i) = F(6);
TF(i) = F(7);
TCL(i) = F(8);
PICL(i) = F(9);
TTL(i) = F(10);
PITL(i) = F(11);
f_(i) = F(12);
Thrust(i) = F(13);
S_(i) = F(14);
m0_(i) = F(15);
T4(i) = F(16);
Pic(i) = F(17);
i = i+1;
end
i = 1;
if j == 1 && k == 1
    ThrustT1A1 = Thrust;
    ST1A1 = S_;
    pi_cHT1A1 = PICH;
    pi_cLT1A1 = PICL;
    pi_fT1A1 = PI_F;
    alphaT1A1 = Alpha;
    m0T1A1 = m0_;
    T4T1A1 = T4;
    PicT1A1 = Pic;
elseif j ==2 && k ==1
    ThrustT2A1 = Thrust;
    ST2A1 = S_;
    pi_cHT2A1 = PICH;
    pi_cLT2A1 = PICL;
    pi_fT2A1 = PI_F;
    alphaT2A1 = Alpha;
    m0T2A1 = m0_;
    T4T2A1 = T4;
    PicT2A1 = Pic;
elseif j == 3 && k ==1
    ThrustT3A1 = Thrust;
    ST3A1 = S_;
    pi_cHT3A1 = PICH;
    pi_cLT3A1 = PICL;
    pi_fT3A1 = PI_F;
    alphaT3A1 = Alpha;
    m0T3A1 = m0_;
    T4T3A1 = T4;
     PicT3A1 = Pic;
elseif j ==4 && k ==1
    ThrustT4A1 = Thrust;
    ST4A1 = S_;
    pi_cHT4A1 = PICH;
    pi_cLT4A1 = PICL;
    pi_fT4A1 = PI_F;
    alphaT4A1 = Alpha;
    m0T4A1 = m0_;
    T4T4A1 = T4;
    PicT4A1 = Pic;
elseif j == 1 && k == 2
    ThrustT1A2 = Thrust;
    ST1A2 = S_;
    pi_cHT1A2 = PICH;
    pi_cLT1A2 = PICL;
    pi_fT1A2 = PI_F;
    alphaT1A2 = Alpha;
    m0T1A2 = m0_;
    T4T1A2 = T4;
     PicT1A2 = Pic;
elseif j ==2 && k ==2
    ThrustT2A2 = Thrust;
    ST2A2 = S_;
    pi_cHT2A2 = PICH;
    pi_cLT2A2 = PICL;
    pi_fT2A2 = PI_F;
    alphaT2A2 = Alpha;
    m0T2A2 = m0_;
    T4T2A2 = T4;
     PicT2A2 = Pic;
elseif j == 3 && k ==2
    ThrustT3A2 = Thrust;
    ST3A2 = S_;
    pi_cHT3A2 = PICH;
    pi_cLT3A2 = PICL;
    pi_fT3A2 = PI_F;
    alphaT3A2 = Alpha;
    m0T3A2 = m0_;
    T4T3A2 = T4;
    PicT3A2 = Pic;
elseif j ==4 && k ==2
    ThrustT4A2 = Thrust;
    ST4A2 = S_;
    pi_cHT4A2 = PICH;
    pi_cLT4A2 = PICL;
    pi_fT4A2 = PI_F;
    alphaT4A2 = Alpha;
    m0T4A2 = m0_;
    T4T4A2 = T4;
    PicT4A2 = Pic;
elseif j == 1 && k == 3
    ThrustT1A3 = Thrust;
    ST1A3 = S_;
    pi_cHT1A3 = PICH;
    pi_cLT1A3 = PICL;
    pi_fT1A3 = PI_F;
    alphaT1A3 = Alpha;
    m0T1A3 = m0_;
    T4T1A3 = T4;
    PicT1A3 = Pic;
elseif j ==2 && k ==3
    ThrustT2A3 = Thrust;
    ST2A3 = S_;
    pi_cHT2A3 = PICH;
    pi_cLT2A3 = PICL;
    pi_fT2A3 = PI_F;
    alphaT2A3 = Alpha;
    m0T2A3 = m0_;
    T4T2A3 = T4;
    PicT2A3 = Pic;
elseif j == 3 && k ==3
    ThrustT3A3 = Thrust;
    ST3A3 = S_;
    pi_cHT3A3 = PICH;
    pi_cLT3A3 = PICL;
    pi_fT3A3 = PI_F;
    alphaT3A3 = Alpha;
    m0T3A3 = m0_;
    T4T3A3 = T4;
    PicT3A3 = Pic;
elseif j ==4 && k ==3
    ThrustT4A3 = Thrust;
    ST4A3 = S_;
    pi_cHT4A3 = PICH;
    pi_cLT4A3 = PICL;
    pi_fT4A3 = PI_F;
    alphaT4A3 = Alpha;
    m0T4A3 = m0_;
    T4T4A3 = T4;
    PicT4A3 = Pic;
else

end


delete = size(Thrust);
 Thrust = zeros(delete);S_ = zeros(delete);PICH = zeros(delete);PICL = zeros(delete);PI_F = zeros(delete);Alpha = zeros(delete);m0_ = zeros(delete); T4 = zeros(delete); Pic = zeros(delete);
   j = j+1; 

end
j = 1;
k = k+1;
end

figure(1) 
plot(M,ThrustT1A1,'r')
hold on
plot(M,ThrustT1A2,'b')
hold on
plot(M,ThrustT1A3,'g')
hold on
plot(M,ThrustT2A1, 'r--')
hold on
plot(M,ThrustT2A2,'b--')
hold on
plot(M,ThrustT2A3,'g--')
hold on
plot(M,ThrustT3A1,'r:')
hold on
plot(M,ThrustT3A2,'b:')
hold on
plot(M,ThrustT3A3,'g:')
hold on
plot(M,ThrustT4A1,'r-.')
hold on
plot(M,ThrustT4A2,'b-.')
hold on
plot(M,ThrustT4A3,'g-.')
hold on
xlabel('Mach Number')
ylabel('Thrust lbf')
title('Thrust(lbf) vs. Mach')
figure(2)
plot(M,ST1A1,'r')
hold on
plot(M,ST1A2,'b')
hold on
plot(M,ST1A3,'g')
hold on
plot(M,ST2A1,'r-')
hold on
plot(M,ST2A2,'b-')
hold on
plot(M,ST2A3,'g-')
hold on
plot(M,ST3A1,'r--')
hold on
plot(M,ST3A2,'b--')
hold on
plot(M,ST3A3,'g--')
hold on
plot(M,ST4A1,'r-.')
hold on
plot(M,ST4A2,'b-.')
hold on
plot(M,ST4A3,'g-.')
hold on
xlabel('Mach Number')
ylabel('TSFC lbm/lbfhr')
title('TSFC vs. Mach')
figure(3)
plot(ONE,PicT1A1,'r')
hold on
plot(theta2,PicT1A2,'b')
hold on
plot(theta3,PicT1A3,'g')
hold on
plot(ONE,PicT2A1,'r--')
hold on
plot(theta2,PicT2A2,'b--')
hold on
plot(theta3,PicT2A3,'g--')
hold on
plot(ONE,PicT3A1,'r:')
hold on
plot(theta2,PicT3A2,'b:')
hold on
plot(theta3,PicT3A3,'g:')
hold on
plot(ONE,PicT4A1,'r-.')
hold on
plot(theta2,PicT4A2,'b-.')
hold on
plot(theta3,PicT4A3,'g-.')
hold on
xlabel('θ')
ylabel('πc')
title('πc vs θ')
figure(4)
plot(M,pi_fT1A1,'r')
hold on
plot(M,pi_fT1A2,'b')
hold on
plot(M,pi_fT1A3,'g')
hold on
plot(M,pi_fT2A1,'r--')
hold on
plot(M,pi_fT2A2,'b--')
hold on
plot(M,pi_fT2A3,'g--')
hold on
plot(M,pi_fT3A1,'r:')
hold on
plot(M,pi_fT3A2,'b:')
hold on
plot(M,pi_fT3A3,'g:')
hold on
plot(M,pi_fT4A1,'r-.')
hold on
plot(M,pi_fT4A2,'b-.')
hold on
plot(M,pi_fT4A3,'g-.')
hold on
xlabel('Mach Number')
ylabel('πf')
title('πf vs. Mach')
figure(5)
plot(M,pi_cHT1A1,'r')
hold on
plot(M,pi_cHT1A2,'b')
hold on
plot(M,pi_cHT1A3,'g')
hold on
plot(M,pi_cHT2A1,'r--')
hold on
plot(M,pi_cHT2A2,'b--')
hold on
plot(M,pi_cHT2A3,'g--')
hold on
plot(M,pi_cHT3A1,'r:')
hold on
plot(M,pi_cHT3A2,'b:')
hold on
plot(M,pi_cHT3A3,'g:')
hold on
plot(M,pi_cHT4A1,'r-.')
hold on
plot(M,pi_cHT4A2,'b-.')
hold on
plot(M,pi_cHT4A3,'g-.')
xlabel('Mach Number')
ylabel('πcH')
title('πcH vs. Mach')
figure(6) 
plot(M,T4T1A1,'r')
hold on
plot(M,T4T1A2,'b')
hold on
plot(M,T4T1A3,'g')
hold on
plot(M,T4T2A1,'r--')
hold on
plot(M,T4T2A2,'b--')
hold on
plot(M,T4T2A3,'g--')
hold on
plot(M,T4T3A1,'r:')
hold on
plot(M,T4T3A2,'b:')
hold on
plot(M,T4T3A3,'g:')
hold on
plot(M,T4T4A1,'r-.')
hold on
plot(M,T4T4A2,'b-.')
hold on
plot(M,T4T4A3,'g-.')
hold on
xlabel('Mach Number')
ylabel('Tt4(R)')
title('Tt4(R) vs. Mach')
figure(7)
plot(M,alphaT1A1,'r')
hold on
plot(M,alphaT1A2,'b')
hold on
plot(M,alphaT1A3,'g')
hold on
plot(M,alphaT2A1,'r--')
hold on
plot(M,alphaT2A2,'b--')
hold on
plot(M,alphaT2A3,'g--')
hold on
plot(M,alphaT3A1,'r:')
hold on
plot(M,alphaT3A2,'b:')
hold on
plot(M,alphaT3A3,'g:')
hold on
plot(M,alphaT4A1,'r-.')
hold on
plot(M,alphaT4A2,'b-.')
hold on
plot(M,alphaT4A3,'g-.')
hold on
xlabel('Mach Number')
ylabel('α')
title('α vs Mach')
figure(8)
plot(M,m0T1A1,'r')
hold on
plot(M,m0T1A2,'b')
hold on
plot(M,m0T1A3,'g')
hold on
plot(M,m0T2A1,'r--')
hold on
plot(M,m0T2A2,'b--')
hold on
plot(M,m0T2A3,'g--')
hold on
plot(M,m0T3A1,'r:')
hold on
plot(M,m0T3A2,'b:')
hold on
plot(M,m0T3A3,'g:')
hold on
plot(M,m0T4A1,'r-.')
hold on
plot(M,m0T4A2,'b-.')
hold on
plot(M,m0T4A3,'g-.')
hold on
xlabel('Mach Number')
ylabel('m0 (lbm/s')
title('mass flow rate(lbm/s) vs. Mach')
%F [tcH,pi_cH,pi_f,M19,M9,alpha,tf,tcL,pi_cL,ttL,pi_tL,f,Thrust,S]


function [Data] = OffDesign(M0,T0,Tt4,p0)

global M_0R height m0R alphaR gamc cpc gamt cpt pi_d pi_b pi_n pi_fn pi_cR pi_cLR pi_fR eta_b eta_mL eta_mH e_cH e_cL e_tL e_tH e_f P_Pstd T_Tstd gc T0R hpr M R_c R_t a_0R V_0R trR eta_r ...
    pi_rR Tt0R Tt1R Tt2R Tt3R Tt4R Tt5R Tt6R Tt7R Tt8R Tt9R tcLR eta_cLR tlambR Tt2_5R pi_cHR tcHR eta_cHR tfR fR Tt13R eta_fR ttHR pi_tHR eta_tHR Tt4_5R ttLR pi_tLR eta_tLR pt9_p9R p0_p9R M9R T9_T0R ...
    V9_a0R V9_V0R pt19_p19R p0_p19R M19R Tt19R T19_T0R V19_a0R V19_V0R pi_d_max tlamb a0 V0 tr pi_r Tt0 Tt1 eta_r

T0 = T0;
p0 = p0;
err = 0;

%Test Params
p0R = 14.7;

while err < 1
    pi_tL = pi_tLR;
    ttL = ttLR;
    tf = tfR;
    a0 = sqrt(gamc*R_c*gc*T0);
    V0 = M0*a0;
    tr = 1 + (gamc-1)*(M0^2)/2;
    pi_r = tr^(gamc/(gamc-1));
    Tt0 = T0 * tr;
    Tt1 = Tt0;
    eta_r = 1;
    pi_d = pi_d_max;
    tlamb = (cpt*Tt4)/(cpc*T0);
    tf0 = tfR; pi_tL0 = pi_tLR;
     x0 = [tf0,pi_tL0];
    Iter = @(x) iterate(x,T0);
   
    g = fsolve(Iter,x0,optimset('Diagnostics','off', 'Display','off'));
    tf = g(1); pi_tL = g(2);
    tcH  = 1 + (( (tlamb/tr) / (tlambR/trR) )*(tfR/tf)* (tcHR -1));
    %ttL x(2)
    ttL= 1 - eta_tLR * (1-(pi_tL)^((gamt-1)/gamt) );
    %pi_cH x(3)
    pi_cH = (1 + eta_cHR*(tcH-1))^(gamc/(gamc-1));
    %pi_f x(4)
    pi_f = (1+eta_fR*(tf-1))^(gamc/(gamc-1));
    pi_tH = pi_tHR;
    ttH = ttHR;
    pi_cL = pi_f;
    pi_cR = pi_cHR*pi_cLR;
    Tt2 = T0*tr;
    K_1 = .7519/Tt4R * (pi_cR ^ (gamc-1)/gamc -1);
    pi_cMax = (1 + Tt4/(Tt2/520) *K_1)^ (gamc / (gamc-1));
    if pi_cMax > 36
        Tt4 = Tt4 - 5;
    else 
      err = 1;
    end

end
err =0;
    %M19 x(5)
    p0_p9_crit = (((gamt+1)/2)^(gamt/(gamt-1)))/(pi_r*pi_d*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_n);
    if p0_p9_crit <= 1
        pt9_p9 = (((gamt+1)/2)^(gamt/(gamt-1)));
        p0_p9 = (((gamt+1)/2)^(gamt/(gamt-1)))/(pi_r*pi_d*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_n);
        M9 = 1;
    else
        pt9_p9 = pi_r*pi_d*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_n;
        p0_p9 = 1;
        M9 = sqrt( (2/(gamt-1)) * ((pt9_p9)^((gamt-1)/gamt) - 1));
    end

    T9_T0 = tlamb*ttH*ttL*cpc/(cpt*(pt9_p9)^((gamt-1)/gamt));
    a9 = sqrt(gamt*R_t*gc*T9_T0*T0);
    V9 = M9*a9;
    V9_a0 = M9*sqrt(T9_T0*gamt*R_t/(gamc*R_c));
    V9__V0R = V9_a0*(a_0R/V_0R);

    p0_p19_crit = ((gamc+1)/2)^(gamc/(gamc-1))/(pi_r*pi_d*pi_f*pi_fn);
    if p0_p19_crit <= 1
        pt19_p19 = ((gamc+1)/2)^(gamc/(gamc-1));
        p0_p19 = ((gamc+1)/2)^(gamc/(gamc-1))/(pi_r*pi_d*pi_f*pi_fn);
        M19 = 1;
    else
        pt19_p19 = pi_r*pi_d*pi_f*pi_fn;
        p0_p19 = 1;
        M19 = sqrt( (2/(gamc-1)) * ((pt19_p19)^((gamc-1)/gamc) - 1));
    end
    T19_T0 = tr*tf/(pt19_p19^((gamt-1)/gamt));
    V19_a0 = M19*sqrt(gamt*R_t*T19_T0/(gamc*R_c));
    alp1 = alphaR*pi_cHR*pi_cLR*pi_f/(pi_cH*pi_cL*pi_fR);
    alp2 = sqrt(((tlamb/(tf*tr))/(tlambR/(tfR*trR)))*(tfR/tf) );
    alp3 = MFP(gamc,R_c,M19) / MFP(gamc,R_c,M19R);
    alpha = alp1*alp2*alp3;
    tf1 = (1-ttL)/(1-ttLR); tf2 = ((tlamb/tr)/(tlambR/trR)); tf3 = (1+alphaR)/(1+alpha); tf4 = (tfR-1);
    tcL = 1+(tf-1)*(tcLR-1)/(tfR-1);

    %Performance Conditions
    f = (tlamb-tr*tf*tcH)/((eta_b*hpr)/(cpc*T0)-tlamb);
    m0 = m0R*((1+alpha)/(1+alphaR))*((p0*pi_r*pi_d*pi_f*pi_cH)/(p0R*pi_rR*pi_d*pi_fR*pi_cHR))*sqrt(Tt4R/Tt4);
    STC = (a0/((1+alpha)*gc))*((1+f)*V9_a0-M0+(1+f)*R_t*T9_T0*(1-p0_p9)/(R_c*V9_a0*gamc));
    STB = (alpha*a0)/((1+alpha)*gc)*(V19_a0-M0+((T19_T0*(1-p0_p19))/(V19_a0*gamc)));
    ST = STC+STB;
    Thrust = m0*ST;
    S = f*3600/((1+alpha)*ST);
    Tt4Real = Tt4;
    pi_c = pi_cL*pi_cH;

    Data = [tcH,pi_cH,pi_f,M19,M9,alpha,tf,tcL,pi_cL,ttL,pi_tL,f,Thrust,S,m0,Tt4Real,pi_c];
end

 %¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯¯\_(ツ)_/¯

function I = iterate(x0,T0)
global M_0R height m0R alphaR gamc cpc gamt cpt pi_d pi_b pi_n pi_fn pi_cR pi_cLR pi_fR eta_b eta_mL eta_mH e_cH e_cL e_tL e_tH e_f P_Pstd T_Tstd gc T0R hpr M R_c R_t a_0R V_0R trR eta_r ...
    pi_rR Tt0R Tt1R Tt2R Tt3R Tt4R Tt5R Tt6R Tt7R Tt8R Tt9R tcLR eta_cLR tlambR Tt2_5R pi_cHR tcHR eta_cHR tfR fR Tt13R eta_fR ttHR pi_tHR eta_tHR Tt4_5R ttLR pi_tLR eta_tLR pt9_p9R p0_p9R M9R T9_T0R ...
    V9_a0R V9_V0R pt19_p19R p0_p19R M19R Tt19R T19_T0R V19_a0R V19_V0R pi_d_max tlamb a0 V0 tr pi_r Tt0 Tt1 eta_r 
     tf = x0(1);
     pi_tL = x0(2);
     
    tcH  = 1 + (( (tlamb/tr) / (tlambR/trR) )*(tfR/tf)* (tcHR -1));
    %ttL x(2)
    ttL= 1 - eta_tLR * (1-(pi_tL)^((gamt-1)/gamt) );
    %pi_cH x(3)
    pi_cH = (1 + eta_cHR*(tcH-1))^(gamc/(gamc-1));
    %pi_f x(4)
    pi_f = (1+eta_fR*(tf-1))^(gamc/(gamc-1));
    pi_tH = pi_tHR;
    ttH = ttHR;
    pi_cL = pi_f;

    %M19 x(5)
    p0_p9_crit = (((gamt+1)/2)^(gamt/(gamt-1)))/(pi_r*pi_d*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_n);
    if p0_p9_crit <= 1
        pt9_p9 = (((gamt+1)/2)^(gamt/(gamt-1)));
        p0_p9 = (((gamt+1)/2)^(gamt/(gamt-1)))/(pi_r*pi_d*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_n);
        M9 = 1;
    else
        pt9_p9 = pi_r*pi_d*pi_cL*pi_cH*pi_b*pi_tH*pi_tL*pi_n;
        p0_p9 = 1;
        M9 = sqrt( (2/(gamt-1)) * ((pt9_p9)^((gamt-1)/gamt) - 1));
    end

    T9_T0 = tlamb*ttH*ttL*cpc/(cpt*(pt9_p9)^((gamt-1)/gamt));
    a9 = sqrt(gamt*R_t*gc*T9_T0*T0);
    V9 = M9*a9;
    V9_a0 = M9*sqrt(T9_T0*gamt*R_t/(gamc*R_c));
    V9__V0R = V9_a0*(a_0R/V_0R);

    p0_p19_crit = ((gamc+1)/2)^(gamc/(gamc-1))/(pi_r*pi_d*pi_f*pi_fn);
    if p0_p19_crit <= 1
        pt19_p19 = ((gamc+1)/2)^(gamc/(gamc-1));
        p0_p19 = ((gamc+1)/2)^(gamc/(gamc-1))/(pi_r*pi_d*pi_f*pi_fn);
        M19 = 1;
    else
        pt19_p19 = pi_r*pi_d*pi_f*pi_fn;
        p0_p19 = 1;
        M19 = sqrt( (2/(gamc-1)) * ((pt19_p19)^((gamc-1)/gamc) - 1));
    end
    T19_T0 = tr*tf/(pt19_p19^((gamt-1)/gamt));
    V19_a0 = M19*sqrt(gamt*R_t*T19_T0/(gamc*R_c));
    alp1 = alphaR*pi_cHR*pi_cLR*pi_f/(pi_cH*pi_cL*pi_fR);
    alp2 = sqrt(((tlamb/(tf*tr))/(tlambR/(tfR*trR)))*(tfR/tf) );
    alp3 = MFP(gamc,R_c,M19) / MFP(gamc,R_c,M19R);
    alpha = alp1*alp2*alp3;
    tf1 = (1-ttL)/(1-ttLR); tf2 = ((tlamb/tr)/(tlambR/trR)); tf3 = (1+alphaR)/(1+alpha); tf4 = (tfR-1);
    tf_new = 1+(tf1*tf2*tf3*tf4);
    pi_tL_new = pi_tLR * sqrt(ttL/ttLR)*MFP(M9R,gamt,R_t)/MFP(M9,gamt,R_t);
    I(1) = tf-tf_new;
    I(2) = pi_tL - pi_tL_new;
    end
function [mfp] = MFP(M,gam,R)
global gc
mfp = sqrt(gam*gc/R)*M*(1+((gam-1)/2)*M^2)^(-(gam+1)/(2*(gam-1)));
end




