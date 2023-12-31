% reaction rates

% loaded by into matlab by hand (JPK) and paired with states
% separate file makes easy to edit / easy to adapt code if we get better
% values

% binary states resorted (KL) maintaining pairing

% now enter the data  
% these are the reactions in the current code
X=[110000;
    110010;
    110100;
    110110;
    100000;
    100010;
    100100;
    100110;
    110000;
    100000;
    010000;
    000000;
  %  010010;  %this reaction is not relevant
    001000;
    101010;
    000000;
    100010;
    110010;
    110011;
    100010;
    110010;
    100110;
    110110;
    000000;
    010000;
    000001;
    010001;
    100010;
    110010;
    100011;
    110011;
    000100;
    010100;
    000101;
    010101;
    100110;
    110110;
    100111;
    110111;
    010100;
    010101;
    000000;
    010000;
    000100;
    010100;
    100010;
    110010;
    100110;
    110110;
    000001;
    010001;
   000101;
   010101;
   100011;
   110011;
   100111;
   110111;
   000001;
   010001;
   000101;
   010101;
   001000;
   010000;
   101010;
   110010;
   001100;
   010100;
   101110;
   110110];
   
Y=[100000;
    100010;
    100100;
    100110;
    110000;
    110010;
    110100;
    110110;
    010000;
    000000;
    110000;
    100000;
  %  110010;  %this reaction is not relevant
    000000;
    100010;
    001000;
    101010;
    010100;
    010101;
    000001;
    010001;
    000101;
    010101;
    000100;
    010100;
    000101;
    010101;
    100110;
    110110;
    100111;
    110111;
    000000;
    010000;
    000001;
    010001;
    100010;
    110010;
    100011;
    110011;
    110010;
    110011;
    000001;
    010001;
    000101;
    010101;
    100011;
    110011;
    100111;
    110111;
    000000;
    010000;
    000100;
    010100;
    100010;
    110010;
    100110;
    110110;
    100010;
    110010;
    100110;
    110110;
    010000;
    001000;
    110010;
    101010;
    010100;
    001100;
    110110;
    101110];
     
XY=[X,Y];

% now to get the reaction rates
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
%cl101m=1900000;
cl101m=190000;
cl0p11=179;
cl1p11=210;
cl0m11=580000;
cl1m11=290000;
cl011p=77;
cl111p=80;
%cl011m=1100000;
cl011m=110000;
%cl111m=2000000;
cl111m=200000;
cl00pm=680000;
%cl10pm=1000000;
cl10pm=100000;
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
%hpm11=2200000;
hpm11=220000;
hmp11=.23;
u10=18;
d10=2500000;
u00=49;
%d00=490000000;
d00=490000;
u11=16777216;

% put these into the right position:
rxn =[r148d00;
    r148d01;
    r148d10;
    r148d11;
    r148p00;
    r148p01;
    r148p10;
    r148p11;
    d10;
    d00;
    u10;
    u00;
    %u11;%this reaction is not relevant
    r203d0;
    r203d1;
    r203p0;
    r203p1;
    cl1pm0;
    cl1pm1;
    cl00mp;
    cl10mp;
    cl01mp;
    cl11mp;
    cl0p00;
    cl1p00;
    cl0p01;
cl1p01;
cl0p10;
cl1p10;
cl0p11;
cl1p11;
cl0m00;
cl1m00;
cl0m01;
cl1m01;
cl0m10;
cl1m10;
cl0m11;
cl1m11;
cl1mp0;
cl1mp1;
cl000p;
cl100p;
cl010p;
cl110p;
cl001p;
cl101p;
cl011p;
cl111p;
cl000m;
cl100m;
cl010m;
cl110m;
cl001m;
cl101m;
cl011m;
cl111m;
cl00pm;
cl10pm;
cl01pm;
cl11pm;
hpm00;
hmp00;
hpm01;
hmp01;
hpm10;
hmp10;
hpm11;
hmp11];


XYr = [X Y rxn];

% sortrows([X Y rxn])

% writematrix(XYr,'hcl_data.xls') 
 