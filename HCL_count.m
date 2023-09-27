% calculate HCl states
global   pw10 discC
pw = [5,4,3,2,1,0];
pw10=[10^5,10^4,10^3,10^2,10,1];
 

clear Xl
clear nperj

nper=0; % permitted reactions
np=0;  % not permitted counter
  

clear Xnp
clear npj

np=0;

%    'transitions not permitted'
for j1=0:1
    for j2=0:1
        for j3 =0:1

                x(1) = j1;
                y(1) = j1;
                x(2) = 0;
                y(2) = 0;
                x(3) = j2;
                y(3) = j2;
                x(4) = 1;
                y(4) = 0;
                x(5) = 0;
                y(5) =1;
                x(6) = j3;
                y(6) = j3;
 nx=getindex(x);

ny = getindex(y);
 if ((nx>0)&(ny>0))
                   np = np+1;
               
              npj(np,1) = nx;
             Xnp(np,1) = getxn(x);
              npj(np,2)=ny;
              Xnp(np,2) = getxn(y);
 end
        end    
    end

end

np
%    'reactions not permitted'
for j1=0:1
                x(1) = 1;
                y(1) = 1;
                x(2) = 1;
                y(2) = 0;
                x(3) = 0;
                y(3) = 1;
                x(4) = 0;
                y(4) = 0;
                x(5) = 0;
                y(5) = 0;
                x(6) = j1;
                y(6) = j1;
               nx=getindex(x);

ny = getindex(y);
 if ((nx>0)&(ny>0))
                   np = np+1;
               
              npj(np,1) = nx;
             Xnp(np,1) = getxn(x);
              npj(np,2)=ny;
              Xnp(np,2) = getxn(y);
 end
end
np
% 
%   'reactions not permitted'
for j1=0:1
    for j2=0:1


                x(1) = 1;
                y(1) = 1;
                x(2) = 0;
                y(2) = 1;
                x(3) = 1;
                y(3) = 0;
                x(4) = j1;
                y(4) = j1;
                x(5) = 0;
                y(5) =0;
                x(6) = j2;
                y(6) = j2;
               nx=getindex(x);

ny = getindex(y);
 if ((nx>0)&(ny>0))
                   np = np+1;
               
              npj(np,1) = nx;
             Xnp(np,1) = getxn(x);
              npj(np,2)=ny;
              Xnp(np,2) = getxn(y);
 end
    end
end
np

%  'reactions not permitted'
for j1=0:1
    for j2=0:1


                x(1) = 1;
                y(1) = 1;
                x(2) = 0;
                y(2) = 1;
                x(3) = 1;
                y(3) = 0;
                x(4) = 1;
                y(4) =  1;
                x(5) = j1;
                y(5) = j1;
                x(6) = j2;
                y(6) = j2;
               nx=getindex(x);

ny = getindex(y);
 if ((nx>0)&(ny>0))
                   np = np+1;
               
              npj(np,1) = nx;
             Xnp(np,1) = getxn(x);
              npj(np,2)=ny;
              Xnp(np,2) = getxn(y);
 end
    end
end
  'transitions not permitted'
  Xnp
  npj
  np
% 


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
 
XY=[X,Y]

for j = 1:length(X)
 
vdx=getv(X(j));
rX(j,1)=getindex(vdx);
vdx=getv(Y(j));
rX(j,2)=getindex(vdx);
end
rX


% check to see if there are any not permitted reactions
notperm=0;
for j=1:np  % number of not permitted reactions
    xt = Xnp(j,1);
    yt = Xnp(j,2);

    for k = 1:length(X)
        if ((xt==X(k)) & (yt==Y(k)))
            'not permitted'
            k
            xt
            yt
            notperm = notperm+1;
        end
    end
end
notperm

 
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

rX

%rxn=rxn.*(rxn>500)

% get rid of small numbers;


% are there missing reactions?'missing reactants'
nr=0;

for j=1:36
    dx = find(rX(:,1)==j);
    if(isempty(dx)==1)
        nr=nr+1;
       missrxn(nr)= getxn(codebreak(j));

       njmissr(nr) = j;
    end
end

nr

missrxn
njmissr
'missing products'
np=0;

for j=1:36
    dx = find(rX(:,2)==j);
    if(isempty(dx))
       np=np+1;

        missprod(np)=getxn(codebreak(j));
        njmissp(np) = j;
    end
end
np
missprod
njmissp

 

% % now populate the transition matrix
clear Amat B C
Amat=zeros(36,36);
length(X)
for j = 1:length(X)

    Amat(rX(j,2),rX(j,1)) = rxn(j); 
end
 

incldcols=find(sum(Amat)>0)
incldrows = find(sum(Amat'>0))
B=Amat(incldrows,:);
  
C=B(:,incldrows)
 
Cd = sum(C)
nmr = length(incldrows)
for j = 1:nmr
    C(j,j) = -Cd(j);
end
 C;

 % check for reducibility

 sc = max(Cd)+10000;
 discC= diag(ones(nmr,1)) +C/sc;
 size(discC);
 mc = dtmc(discC);
 isreducible(mc)
% 
%  C(nmr,:)=ones(1,nmr);
% RHS = zeros(nmr,1);
% RHS(nmr) = 1;
% 
% 
%  [L,U] = lu(C)


figure(1)
graphplot(mc)

codebreak(incldrows(20))
codebreak(incldrows(6))
incldrows(20)
incldrows(6)
p0=zeros(24,1);
p0(1) = 1;

p1=iterateC(p0,50000000)

% check for convergence
C*p1

% now find the flux
% next...

figure(2)
semilogy(p1,'*')
 function gv=getv(V)
xt=V;
for j = 1:6
    gv(j) = fix(xt/(10^(6-j)));
xt=xt-gv(j)*10^(6-j);
end
 end
     
function ndx = getindex(x)
 
ndx=12*(2*x(2)+x(3))+4*(x(1)+x(5))+2*x(4)+x(6)+1; 

if(x(1)==0&x(5)==1)
    ndx = 0;
end

if(x(2)==1&x(3)==1)
    ndx = 0;
end
 
end  

function xn = getxn(x)
global pw10

xn=sum(x.*pw10);

end

function xout = codebreak(ndx)
clear x

for j1=0:1
    for j2=0:1
        for j3=0:1
            for j4=0:1
                for j5=0:1
                    for j6=0:1
     x(1) = j1;
     x(2) = j2;
     x(3) = j3;
     x(4) = j4;
     x(5) = j5;
     x(6) = j6;
     tst=getindex(x);
     if(tst==ndx)
         xout=x;
     end
                    end
                end
            end
        end
    end
end
end

function p1 = iterateC(p0,n)
global discC
for j = 1:n
    p1=discC*p0;
p0 = p1;
end
end
