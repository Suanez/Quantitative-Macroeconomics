% Final Proyect Quantitaive Macroeoconomics
% Daniel Suañez Ramirez
% In this code we carry out the model, solved in the paper of Finla Project of Quantitative Macroeconomics, the model is
% a New Keynesian Model with a technology shock, we approximate to the
% steady state using a second order Taylor approximation. We follow the
% paper of Lombardo Sutherland, 2007, to develop this project. We
% introduces a shock in the productivity of workers.


clear all
clc
tic

%% Parameters
%Parameters
beta = 0.99; % discoutn rate
sig = 2; %elasticity os subtitution of good
phi = 1.5; %Elasticity of sustitution labour demand
epsi = 1.1; %elasticity of substitution of demand
varphi = 0.75; % probability that the firm change princes
A_star = 1.1; % technology in the new steady state
rho = 0.9; %Shock
var = 0.001;  %variance shock
i_bar = 0.06; % target real interes rate
h_pi =2;  %measure of the percentage deviation of inflation fortheir target rate of inflation
% Here we change our Teylor rule, so that the Central Bank only care the
% inflation



for i=0:1;
h_y =0.8*i;    %measure of thepercentage deviation of the NKM output and the flexible price output of RBC model

%randn('seed',0);
shock = 1; %A stn dev (1), random (2)

%Steady State
R_ss = 1/beta;
A_ss = A_star;
f_pi_ss = @(K) (1+i_bar+h_y*K-R_ss)/(R_ss-h_pi);
f_i_ss = @(K) i_bar+h_pi*f_pi_ss(K)+h_y*K;
f_x_ss = @(K) ((1-(1-varphi)*(1+f_pi_ss(K))^(epsi-1))/varphi)^(1/(1-epsi));
f_w_ss = @(K) (epsi-1)/epsi*f_x_ss(K)*A_ss;
f_N_ss = @(K) f_w_ss(K)^(1/(phi+sig))*A_ss^(-sig/(phi+sig));
f_c_ss = @(K) A_ss*f_N_ss(K);
f_K = @(K) f_c_ss(K)/((epsi-1)/epsi*A_ss^(1+phi))^(1/(sig+phi))-1-K;
KK = fzero(f_K,0);
pi_ss = f_pi_ss(KK);
i_ss = f_i_ss(KK);
x_ss = f_x_ss(KK);
w_ss = f_w_ss(KK);
N_ss = f_N_ss(KK);
c_ss = f_c_ss(KK);
chi1_ss = c_ss/(1-beta*(1-varphi)*(1+pi_ss)^epsi);
chi2_ss = c_ss*w_ss/A_ss/(1-beta*(1-varphi)*(1+pi_ss)^epsi);

% We introudces the first and second order condition matrix
%varibles: c,N,w,R,A,i,pi,x,chi1,chi2
cntrl = [1 2 3 4 6 7 8 9 10];
state = [5];
n_cntrl = length(cntrl);
n_state = length(state);
n = n_cntrl+n_state;
nn_state = n_state*(n_state+1)/2;

F = zeros(n,1);
D1F = zeros(n,n);
D2F = zeros(n,n,n);

G = zeros(n,1);
D1G = zeros(n,n);
D2G = zeros(n,n,n);

G(1) = [phi*log(N_ss)-log(w_ss)+sig*log(c_ss)];
D1G(1,[1 2 3]) = [sig phi -1];

G(2) = [log(A_ss)+log(N_ss)-log(c_ss)];
D1G(2,[1 2 5]) = [-1 1 1];

G(3) = [log(1+pi_ss)+log(R_ss)-log(1+i_ss)];
D1G(3,[4 6 7]) = [1 -1 1];

G(4) = [log(epsi/(epsi-1))+log(chi2_ss)-log(chi1_ss)-log(x_ss)];
D1G(4,[8 9 10]) = [-1 -1 1];

G(5) = [1-varphi*exp((1-epsi)*log(x_ss))-(1-varphi)*exp((epsi-1)*log(1+pi_ss))];
D1G(5,[7 8]) = [-(1-varphi)*exp((epsi-1)*log(1+pi_ss))*(epsi-1) -varphi*exp((1-epsi)*log(x_ss))*(1-epsi)];
D2G(5,7,[7]) = [-(1-varphi)*exp((epsi-1)*log(1+pi_ss))*(epsi-1)^2];
D2G(5,8,[8]) = [-varphi*exp((1-epsi)*log(x_ss))*(1-epsi)^2];

G(6) = [rho*log(A_ss)+(1-rho)*log(A_star)];
D1G(6,[5]) = [rho];

F(6) = [log(A_ss)];
D1F(6,[5]) = [1];

G(7) = [exp(-sig*log(c_ss))*exp(-log(R_ss))];
D1G(7,[1 4]) = [-sig*exp(-sig*log(c_ss))*exp(-log(R_ss)) ...
-exp(-sig*log(c_ss))*exp(-log(R_ss))];
D2G(7,1,[1 4]) = [sig^2*exp(-sig*log(c_ss))*exp(-log(R_ss)) ...
sig*exp(-sig*log(c_ss))*exp(-log(R_ss))];
D2G(7,4,[4]) = [exp(-sig*log(c_ss))*exp(-log(R_ss))];

F(7) = [beta*exp(-sig*log(c_ss))];
D1F(7,[1]) = [beta*(-sig)*exp(-sig*log(c_ss))];
D2F(7,1,[1]) = [beta*sig^2*exp(-sig*log(c_ss))];

G(8) = [exp(log(chi1_ss))-exp(log(c_ss))];
D1G(8,[1 9]) = [-exp(log(c_ss)) exp(log(chi1_ss))];
D2G(8,1,[1]) = [-exp(log(c_ss))];
D2G(8,9,[9]) = [exp(log(chi1_ss))];

F(8) = [beta*(1-varphi)*exp(epsi*log(1+pi_ss)+log(chi1_ss))];
D1F(8,[7 9]) = [beta*(1-varphi)*exp(epsi*log(1+pi_ss)+log(chi1_ss))*epsi ...
beta*(1-varphi)*exp(epsi*log(1+pi_ss)+log(chi1_ss))];
D2F(8,7,[7 9]) = [beta*(1-varphi)*exp(epsi*log(1+pi_ss)+log(chi1_ss))*epsi^2 ...
beta*(1-varphi)*exp(epsi*log(1+pi_ss)+log(chi1_ss))*epsi];
D2F(8,9,[9]) = [beta*(1-varphi)*exp(epsi*log(1+pi_ss)+log(chi1_ss))];

G(9) = [exp(log(chi2_ss))-exp(log(c_ss)+log(w_ss)-log(A_ss))];
D1G(9,[1 3 5 10]) = [-exp(log(c_ss)+log(w_ss)-log(A_ss)) -exp(log(c_ss)+log(w_ss)-log(A_ss)) ...
exp(log(c_ss)+log(w_ss)-log(A_ss)) exp(log(chi2_ss))];
D2G(9,1,[1 3 5]) = [-exp(log(c_ss)+log(w_ss)-log(A_ss)) -exp(log(c_ss)+log(w_ss)-log(A_ss)) ...
exp(log(c_ss)+log(w_ss)-log(A_ss))];
D2G(9,3,[3 5]) = [-exp(log(c_ss)+log(w_ss)-log(A_ss)) exp(log(c_ss)+log(w_ss)-log(A_ss))];
D2G(9,5,[5]) = [-exp(log(c_ss)+log(w_ss)-log(A_ss))];
D2G(9,10,[10]) = [exp(log(chi2_ss))];

F(9) = [beta*(1-varphi)*exp(epsi*log(1+pi_ss)+log(chi2_ss))];
D1F(9,[7 10]) = [beta*(1-varphi)*exp(epsi*log(1+pi_ss)+log(chi2_ss))*epsi ...
beta*(1-varphi)*exp(epsi*log(1+pi_ss)+log(chi2_ss))];
D2F(9,7,[7 10]) = [beta*(1-varphi)*exp(epsi*log(1+pi_ss)+log(chi2_ss))*epsi^2 ...
beta*(1-varphi)*exp(epsi*log(1+pi_ss)+log(chi2_ss))*epsi];
D2F(9,10,[10]) = [beta*(1-varphi)*exp(epsi*log(1+pi_ss)+log(chi2_ss))];

G(10) = 1+i_bar+h_pi*exp(log(1+pi_ss))-h_pi+h_y*(exp(log(c_ss)-1/(phi+sig)*(log((epsi-1)/epsi)+...
(1+phi)*log(A_ss)))-1)-exp(log(1+i_ss));
D1G(10,[1 5 6 7]) = [h_y*(exp(log(c_ss)-1/(phi+sig)*(log((epsi-1)/epsi)+(1+phi)*log(A_ss)))-1) ...
-h_y*(exp(log(c_ss)-1/(phi+sig)*(log((epsi-1)/epsi)+(1+phi)*log(A_ss)))-1)*(1+phi)/(phi+sig) ...
-exp(log(1+i_ss)) h_pi*exp(log(1+pi_ss))];
D2G(10,1,[1 5]) = [h_y*(exp(log(c_ss)-1/(phi+sig)*(log((epsi-1)/epsi)+(1+phi)*log(A_ss)))-1) ...
-h_y*(exp(log(c_ss)-1/(phi+sig)*(log((epsi-1)/epsi)+(1+phi)*log(A_ss)))-1)*(1+phi)/(phi+sig)];
D2G(10,5,[5]) = [h_y*(exp(log(c_ss)-1/(phi+sig)*(log((epsi-1)/epsi)+(1+phi)*log(A_ss)))-1)*(1+phi)^2/(phi+sig)^2];
D2G(10,6,[6]) = [-exp(log(1+i_ss))];
D2G(10,7,[7]) = [h_pi*exp(log(1+pi_ss))];


%varibles: c,N,w,R,A,i,pi,x,chi1,chi2

r = 1/beta-1;
sta = abs(F)<1e-8;
F(sta) = G(sta)/(1+r);
D1F(sta,:) = D1G(sta,:)/(1+r);
D2F(sta,:,:) = D2G(sta,:,:)/(1+r);

Jf = D1F;
Jg = D1G;
g = G-F;
Hf = rect_D2(D2F);
Hg = rect_D2(D2G);

M = Jf\Jg;
m = Jf\g;

[Theta,theta] = Linear_solution_cntrl(M,m,cntrl,state,r-1e-6);
[Omega,omega] = Linear_solution_next_state(M,m,cntrl,state,r-1e-6);

SIG = [var];
chi(cntrl,:) = Theta;
chi(state,:) = eye(n_state);
Gchi = Gamma2(chi);

Gp = Gamma2(ones(1,n));
Hf = Hf.*Gp;
Hg = Hg.*Gp;

Jf2 = [Jf 0.5*Hf*Gchi;zeros(nn_state,n) eye(nn_state)];
Jg2 = [Jg 0.5*Hg*Gchi;zeros(nn_state,n) Gamma2(Omega)];
g2 = [g;vec_var(SIG)];

M2 = Jf2\Jg2;
m2 = Jf2\g2;

state2 = [state n+1:n+nn_state];

[Theta2,theta2] = Linear_solution_cntrl(M2,m2,cntrl,state2,r-1e-6);
[Omega2,omega2] = Linear_solution_next_state(M2,m2,cntrl,state2,r-1e-6);
theta2 = real(theta2);
omega2 = real(omega2);

X0 = [sqrt(var)];
T = 12;

T = 12;
epsilon = zeros(T,1);
epsilon(3) = sqrt(var);
X0 = [0];

for t = 1:T
  XX0 = [X0; vec_var(X0*X0')];
  Y0 = Theta2*XX0;
  variables(t,cntrl) = Y0';
  variables(t,state) = X0';
  X0 = Omega2(1:n_state,:)*XX0+[epsilon(t)];
end

%% We plot our results
%varibles: c,N,w,R,A,i,pi,x,Phi1,Phi2
figure(2)
hold on
subplot(8,1,1)
plot(variables(:,[1]),'LineWidth',2);title('Consumption');grid;
hold on
subplot(8,1,2)
plot(variables(:,[2]),'LineWidth',2);title('Labour Supply');grid;
hold on
subplot(8,1,3)
plot(variables(:,[3]),'LineWidth',2);title('Real Wage');grid;
hold on
subplot(8,1,4)
plot(variables(:,[4]),'LineWidth',2);title('Real Interest Rate');grid;
hold on
subplot(8,1,5)
plot(variables(:,[5]),'LineWidth',2);title('Technology');grid;
hold on
subplot(8,1,6)
plot(variables(:,[6]),'LineWidth',2);title('Nominal Interest Rate');grid;
hold on
subplot(8,1,7)
plot(variables(:,[7]),'LineWidth',2);title('Inflation');grid;

%El hold on lo pongo para que cuando repita no te sobreescriba
toc

end
legend('h_y=0','h_y=0.8')