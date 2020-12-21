% This function takes the policy out of the state's in the following period

function [Omega,omega] = Linear_solution_next_state(M,m,cntrl,state,r)

[Theta,theta] = Linear_solution_cntrl(M,m,cntrl,state,r);

Omega = M(state,state)+M(state,cntrl)*Theta;
omega = m(state)+M(state,cntrl)*theta;