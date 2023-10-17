% analysis

% set up transition matrix
% hcl_count

% create Markov Chain using transition matrix from hcl_count
sc = max(Cd)+10000;
 discC= diag(ones(nmr,1)) +C/sc;
 size(discC);
 mc = dtmc(discC);
 mc.StateNames = ["Depression" "Recession"];
 isreducible(mc)
% 
%  C(nmr,:)=ones(1,nmr);
% RHS = zeros(nmr,1);
% RHS(nmr) = 1;
% 
% 
%  [L,U] = lu(C)


figure(1)
graphplot(mc,'ColorEdges',true)

% find indexing / binary rep. of isolated states
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

% where flux(i) = sum_j=1^n p1(i)*discC(i,j)
flux = discC'*p1;



asymptotics(mc)

figure(2)
semilogy(p1,'*')




function p1 = iterateC(p0,n)
% discrete Euler simulations
global discC
for j = 1:n
    p1=discC*p0;
p0 = p1;
end
end