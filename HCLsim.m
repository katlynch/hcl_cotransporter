%%Code for an HCL TRANSPORTER

%%Things that are still needed:
%%% Concentration variables and ODEs:
%H Internal
%H External
%CL Internal
%Cl external

%%%Terms that will need to be included into the concentration ODES.
%any term that has the following rate accompanying it:
%Cl-p--,Cl---p,r203p, r148p: these will be negative terms in conc. ODES
%Cl-m--,Cl---m,r203d, r148d: these will be positive terms in conc. ODES

%%%Some Transition Rates must be modified by Corresponding Concentrations
%The following Rates will be conc. dep. and should be multiplied by concentration variable:
%All Cl-p-- rates: external CL conc.
%All Cl---p rates: internal CL conc.
%All r203p rates: Internal H conc.
%All r148p rates: External H conc. 

function [time,S]=voltageclamp

%%%Initial Conditions

%%%Sets all states to 0 except the full empty/down state.
S0=zeros(1,48);
S0(1)=1;
%%%



%%%%Set up time interval and solve the ODE
t_int=[0:0.1:10];
[time,S]=ode45(@deRHS,t_int,S0);
S(end,:);



end


function s_prime=deRHS(t,s) 


%%%%%% Parameters for Model
%%%Descriptions:
%%%initial tag on rate:
%cl-chloride ion rate, four indexes are (148,external,central, internal)
%r203-H ion rate at site 203, indexes are(203 site, cent. cl site) 
%r148-H ion rate at site 148, indexes are (148 site, ext.l cl, cent. cl)
%h- H ion transfer between sites, indexes (148, 203, ext.l cl, cent. cl)
%u/d: movement of 148 up/down, indexes are (148, cent. cl)

%%%Rate index possibilities
%m/d: minus ion/deprotonation of site
%p: gain ion/protonation of site
%0/1: Site is empty/occupied


cl0m00=520000;
cl1m00=260000;
cl000p=88;
cl100p=123;
cl0p00=279;
cl1p00=596;
cl000m=370000;
cl100m=720000;
cl0p01=174;
cl1p01=415;
cl010p=62;
cl110p=88;
cl0m01=540000;
cl1m01=270000;
cl010m=390000;
cl110m=780000;
cl001p=83;
cl101p=91;
cl0p10=111;
cl1p10=211;
cl0m10=560000;
cl1m10=350000;
cl001m=980000;
cl101m=1900000;
cl0p11=179;
cl1p11=210;
cl0m11=580000;
cl1m11=290000;
cl011p=77;
cl111p=80;
cl011m=1100000;
cl111m=2000000;
cl00pm=680000;
cl10pm=1000000;
cl00mp=10000;
cl10mp=22000;
cl1pm0=730;
cl1pm1=870;
cl1mp0=410;
cl1mp1=9500;
cl01pm=360000;
cl11pm=310000;
cl01mp=2400;
cl11mp=1200;
r203p0=10000;
r203p1=10000;
r203d0=1000;
r203d1=1000;
r148d00=1100;
r148d01=1.8;
r148p00=430;
r148p01=3000;
r148d10=1.4;
r148d11=.0023;
r148p10=430;
r148p11=3000;
hpm00=240;
hmp00=.00029;
hpm01=640000;
hmp01=.77;
hpm10=790;
hmp10=.000088;
hpm11=2200000;
hmp11=.23;
u10=18;
d10=2500000;
u00=49;
d00=490000000;
u11=16777216;



%%%%%Variable descriptions
%6 indexes-all either 0's or 1's
%Index 1: E148 up 1 or down 0
%The remainder indexes are occupied 1 or not 0
%Index 2: E148- H external
%Index 3: E203- H Internal
%Index 4: Sout- CL external
%Index 5: SCen- Cl central site
%Index 6: Sin-  Cl internal

%%%%Initialization of State variables
%E148 down
%  the state p0xxx1x is not allowed
% the state pxx11xxx is also not allowed
p000000=s(1);
p000001=s(2);
p000100=s(3);
p000101=s(4);
p001000=s(5);
p001001=s(6);
p001100=s(7);
p001101=s(8);
p010000=s(9);
p010001=s(10);
p010100=s(11);
p010101=s(12);
p011000=s(13);
p011001=s(14);
p011100=s(15);
p011101=s(16);

%E148 UP
p100000=s(17);
p100001=s(18);
p100010=s(19);
p100011=s(20);
p100100=s(21);
p100101=s(22);
p100110=s(23);
p100111=s(24);
p101000=s(25);
p101001=s(26);
p101010=s(27);
p101011=s(28);
p101100=s(29);
p101101=s(30);
p101110=s(31);
p101111=s(32);
p110000=s(33);
p110001=s(34);
p110010=s(35);
p110011=s(36);
p110100=s(37);
p110101=s(38);
p110110=s(39);
p110111=s(40);
p111000=s(41);
p111001=s(42);
p111010=s(43);
p111011=s(44);
p111100=s(45);
p111101=s(46);
p111110=s(47);
p111111=s(48);




%%%RHS of the ODE equations
s_prime=[p000001*cl000m+p000100*cl0m00+p001000*r203d0+p100000*d00-p000000*(cl000p+cl0p00+r203p0+u00),...
         p000000*cl000p+p000101*cl0m01+p001001*r203d0+p100001*d00-p000001*(cl000m+cl0p01+r203p0+u00),...
         p000000*cl0p00+p000101*cl010m+p001100*r203d0+p100100*d00-p000100*(cl0m00+cl010p+r203p0+u00),...
         p000001*cl0p01+p000100*cl010p+p001101*r203d0+p100101*d00-p000101*(cl0m01+cl010m+r203p0+u00),...
         p000000*r203p0+p001001*cl000m+p001100*cl0m00+p010000*hmp00+p101000*d00-p001000*(r203d0+cl000p+cl0p00+hpm00+u00),...
         p000001*r203p0+p001000*cl000p+p001101*cl0m01+p010001*hmp00+p101001*d00-p001001*(r203d0+cl000m+cl0p01+hpm00+u00),...
         p000100*r203p0+p001000*cl0p00+p001101*cl010m+p010100*hmp10+p101100*d00-p001100*(r203d0+cl0m00+cl010p+hpm10+u00),...
         p000101*r203p0+p001001*cl0p01+p001100*cl010p+p010101*hmp10+p101101*d00-p001101*(r203d0+cl0m01+cl010m+hpm10+u00),...
         p001000*hpm00+p010001*cl100m+p010100*cl1m00+p011000*r203d0+p110000*d10-p010000*(hmp00+cl100p+cl1p00+r203p0+u10),...
         p001001*hpm00+p010000*cl100p+p010101*cl1m01+p011001*r203d0+p110001*d10-p010001*(hmp00+cl100m+cl1p01+r203p0+u10),...
         p001100*hpm10+p010000*cl1p00+p010101*cl110m+p011100*r203d0+p110100*d10-p010100*(hmp10+cl1m00+cl110p+r203p0+u10),...
         p001101*hpm10+p010001*cl1p01+p010100*cl110p+p011101*r203d0+p110101*d10-p010101*(hmp10+cl1m01+cl110m+r203p0+u10),...
         p010000*r203p0+p011001*cl100m+p011100*cl1m00+p111000*d10-p011000*(r203d0+cl100p+cl1p00+u10),...
         p010001*r203p0+p011000*cl100p+p011101*cl1m01+p111001*d10-p011001*(r203d0+cl100m+cl1p01+u10),...
         p010100*r203p0+p011000*cl1p00+p011101*cl110m+p111100*d10-p011100*(r203d0+cl1m00+cl110p+u10),...
         p010101*r203p0+p011001*cl1p01+p011100*cl110p+p111101*d10-p011101*(r203d0+cl1m01+cl110m+u10),...
         p100001*cl000m+p100100*cl0m00+p101000*r203d0+p110000*r148d00+p000000*u00-p100000*(d00+r148p00+cl000p+cl0p00+r203p0),...
         p100000*cl000p+p100010*cl00mp+p100101*cl0m01+p101001*r203d0+p110001*r148d00+p000001*u00-p100001*(cl000m+cl0p01+r203p0+r148p00+d00+cl00pm),...
         p100001*cl00pm+p100011*cl001m+p100110*cl0m10+p101010*r203d1+p110010*r148d01-p100010*(cl00mp+cl001p+cl0p10+r203p1+r148p01),...
         p100010*cl001p+p100111*cl0m11+p101011*r203d1+p110011*r148d01-p100011*(cl001m+cl0p11+r203p1+r148p01),...
         p000100*u00+p100000*cl0p00+p100101*cl010m+p101100*r203d0+p110100*r148d10-p100100*(d00+cl0m00+cl010p+r203p0+r148p10),...
         p100001*cl0p01+p000101*u00+p100100*cl010p+p100110*cl01mp+p101101*r203d0+p110101*r148d10-p100101*(d00+cl010m+cl01pm+r203p0+r148p10+cl0m01),...
         p100010*cl0p10+p100101*cl01pm+p100111*cl011m+p101110*r203d1+p110110*r148d11-p100110*(cl0m10+cl01mp+cl011p+r203p1+r148p11),...
         p100011*cl0p11+p100110*cl011p+p101111*r203d1+p110111*r148d11-p100111*(cl0m11+cl011m+r203p1+r148p11),...
         p001000*u00+p100000*r203p0+p101001*cl000m+p101100*cl0m00+p111000*r148d00-p101000*(d00+r203d0+cl000p+cl0p00+r148p00),...
         p001001*u00+p100001*r203p0+p101000*cl000p+p101010*cl00mp+p101101*cl0m01+p111001*r148d00-p101001*(d00+r203d0+cl000m+cl00pm+cl0p01+r148p00),...
         p100010*r203p1+p101001*cl00pm+p101011*cl001m+p101110*cl0m10+p110010*hmp01+p111010*r148d01-p101010*(r203d1+cl00mp+cl001p+cl0p10+hpm01+r148p01),...
         p100011*r203p1+p101010*cl001p+p101111*cl0m11+p110011*hmp01+p111011*r148d01-p101011*(r203d1+cl001m+cl0p11+hpm01+r148p01),...
         p001100*u00+p100100*r203p0+p101000*cl0p00+p101101*cl010m+p110100*hmp10+p111100*r148d10-p101100*(d00+r203d0+cl0m00+cl010p+r148p10),...
         p001101*u00+p100101*r203p0+p101001*cl0p01+p101100*cl010p+p101110*cl01mp+p110101*hmp10+p111101*r148d10-p101101*(d00+r203d0+cl0m01+cl010m+cl01pm+r148p10),...
         p100110*r203p1+p101010*cl0p10+p101101*cl01pm+p101111*cl011m+p110110*hmp11+p111110*r148d11-p101110*(r203d1+cl0m10+cl01mp+cl011p+hpm11+r148p11),...
         p100111*r203p1+p101011*cl0p11+p101110*cl011p+p110111*hmp11+p111111*r148d11-p101111*(r203d1+cl0m11+cl011m+hpm11+r148p11),...
         p010000*u10+p100000*r148p00+p110001*cl100m+p110100*cl1m00+p111000*r203d0-p110000*(d10+r148d00+cl100p+cl1p00+r203p0),...
         p010001*u10+p100001*r148p00+p110000*cl100p+p110010*cl10mp+p110101*cl1m01+p111001*r203d0-p110001*(d10+r148d00+cl100m+cl10pm+cl1p01+r203p0),...
         p100010*r148p01+p101010*hpm01+p110001*cl10pm+p110011*cl101m+p110100*cl1mp0+p110110*cl1m10+p111010*r203d1-p110010*(r148d01+hmp01+cl10mp+cl101p+cl1pm0+cl1p10+r203p1),...
         p100011*r148p01+p101011*hpm01+p110010*cl101p+p110101*cl1mp1+p110111*cl1m11+p111011*r203d1-p110011*(r148d01+hmp01+cl101m+cl1pm1+cl1p11+r203p1),...
         p010100*u10+p100100*r148p10+p110000*cl1p00+p110010*cl1pm0+p110101*cl110m+p111100*r203d0-p110100*(d10+r148d10+cl1m00+cl1mp0+cl110p+r203p0+hmp10),...
         p010101*u10+p100101*r148p10+p110001*cl1p01+p110011*cl1pm1+p110100*cl110p+p111101*r203d0+p110110*cl11mp-p110101*(hmp10+d10+r148d10+cl1m01+cl1mp1+cl110m+r203p0+cl11pm),...
         p100110*r148p11+p101110*hpm11+p110010*cl1p10+p110101*cl11pm+p111110*r203d1+p110111*cl111m-p110110*(r148d11+hmp11+cl1m10+cl11mp+r203p1+cl111p),...
         p100111*r148p11+p101111*hpm11+p110011*cl1p11+p110110*cl111p+p111111*r203d1-p110111*(r148d11+hmp11+cl1m11+cl111m+r203p1),...
         p011000*u10+p101000*r148p00+p110000*r203p0+p111001*cl100m+p111100*cl1m00-p111000*(d10+r148d00+r203d0+cl100p+cl1p00),...
         p011001*u10+p101001*r148p00+p110001*r203p0+p111000*cl100p+p111010*cl10mp+p111101*cl1m01-p111001*(d10+r148d00+r203d0+cl100m+cl10pm+cl1p01),...
         p101010*r148p01+p110010*r203p1+p111001*cl10pm+p111011*cl101m+p111100*cl1mp0+p111110*cl1m10-p111010*(r148d01+r203d1+cl10mp+cl101p+cl1pm0+cl1p10),...
         p101011*r148p01+p110011*r203p1+p111101*cl1mp1+p111111*cl1m11+p111010*cl101p-p111011*(r148d01+r203d1+cl1pm1+cl1p11+cl101m),...
         p011100*u10+p101100*r148p10+p110100*r203p0+p111000*cl1p00+p111010*cl1pm0+p111101*cl110m-p111100*(d10+r148d10+r203d0+cl1m00+cl1mp0+cl110p),...
         p011101*u10+p101101*r148p10+p110101*r203p0+p111001*cl1p01+p111011*cl1pm1+p111100*cl110p+p111110*cl11mp-p111101*(d10+r148d10+r203d0+cl1m01+cl1mp1+cl110m+cl11pm),...
         p101110*r148p11+p110110*r203p1+p111010*cl1p10+p111101*cl11pm+p111111*cl111m-p111110*(r148d11+r203d1+cl1m10+cl11mp+cl111p),...
         p101111*r148p11+p110111*r203p1+p111011*cl1p11+p111110*cl111p-p111111*(r148d11+r203d1+cl1m11+cl111m)]';
  
  
 




end