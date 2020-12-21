% This function takes out the policy function of the control variables 

function [Theta,theta] = Linear_solution_cntrl(M,m,cntrl,state,r)
[n1,n2] = size(M);
n3 = length(m);
nc = length(cntrl);
ns = length(state);

%checks
if n1 ~= n2
    error('matrix not square');
end

if n3 ~= n1
    error('M and m with different dimensions');
end

if nc+ns ~= n1
    error('control plus states different from total variables');
end

%computations

[V,Lam,U] = eig(M);
L = diag(real(Lam));

conv = abs(L) < 1+r;
expl = abs(L) >= 1+r;

ne = sum(expl);
%check
if ne ~= nc
    warning('no saddle path system');
end

if rank(V) == n1
    mm = V\m;
    
    Ze = (eye(ne)-diag(L(expl)))\mm(expl);
    
    Theta = V(cntrl,conv)/V(state,conv);
    theta = (V(cntrl,expl)-Theta*V(state,expl))*Ze;
elseif rank(V) < n1
    U = U';
    mm = U*m;
    
    Ze = (eye(ne)-diag(L(expl)))\mm(expl);
    
    Theta = -U(expl,cntrl)\U(expl,state);
    theta = U(expl,cntrl)\Ze;
end