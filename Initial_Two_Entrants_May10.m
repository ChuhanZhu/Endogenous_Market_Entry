function F = Initial_Two_Entrants_May10(a, beta, xi, p, L_bar, gamma_b, delta, lambda1, lambda2, Delta1, Delta2, l_u_1_shadow, l_u_2_shadow,l_u_1_bank,l_u_2_bank,l_I_1_bank,l_I_2_bank,gamma_s)
% Asset sizes of two intermediaries
a1 = a(1);
a2 = a(2);

%production function
f = @(ai, aj)xi*ai^beta/(aj^beta + ai^beta);
f_prime = @(ai, aj)xi*(beta*ai^(beta-1)/(ai^beta + aj^beta) - beta*ai^(2*beta - 1)/((ai^beta + aj^beta)^2));
%value function of becoming a bank
v_bank = @(ai, aj, l_u_i, l_I_i, lambda_i, Delta_i)f(ai, aj)*(p+(1-p)*lambda_i) - ai + ai*L_bar*gamma_b...
                       - (1-p)*ai*delta*(l_u_i - lambda_i)^2  - ai*l_I_i*Delta_i...
                       + ai*l_I_i*(1-p);
%value function of becoming a shadow bank
v_shadow = @(ai, aj, l_u_i, lambda_i)f(ai,aj)*(p+(1-p)*lambda_i)-ai + ai*l_u_i*gamma_s - (1-p)*ai*delta*(l_u_i - lambda_i)^2;
%FOC of shadow bank
FOC_shadow = @(ai, aj, lambda_i)f_prime(ai, aj)*(p+(1-p)*lambda_i) - 1 + lambda_i*gamma_s + (gamma_s)^2/(4*(1-p)*delta);
%FOC of bank
FOC_bank = @(ai, aj, lambda_i, l_u_i, l_I_i,Delta_i)f_prime(ai, aj)*(p+(1-p)*lambda_i) - 1 + gamma_b*L_bar...
                                           - (1-p)*delta*(l_u_i - lambda_i)^2 + l_I_i*(1-p)...
                                           - l_I_i*Delta_i;

%Two FOCs conditional on intermediaries' choices
% i specific variables: lambda1, lambda2, Delta1, Delta2, l_u_1_shadow, l_u_2_shadow,l_u_1_bank,l_u_2_bank,l_I_1_bank,l_I_2_bank
if v_bank(a1, a2, l_u_1_bank, l_I_1_bank, lambda1, Delta1)> v_shadow(a1, a2, l_u_1_shadow, lambda1) && v_bank(a2, a1, l_u_2_bank, l_I_2_bank, lambda2, Delta2)> v_shadow(a2, a1, l_u_2_shadow, lambda2)
    F(1) = FOC_bank(a1, a2, lambda1, l_u_1_bank, l_I_1_bank, Delta1);
    F(2) = FOC_bank(a2, a1, lambda2, l_u_2_bank, l_I_2_bank, Delta2);
elseif v_bank(a1, a2, l_u_1_bank, l_I_1_bank, lambda1, Delta1)> v_shadow(a1, a2, l_u_1_shadow, lambda1) && v_bank(a2, a1, l_u_2_bank, l_I_2_bank, lambda2, Delta2)< v_shadow(a2, a1, l_u_2_shadow, lambda2)
    F(1) = FOC_bank(a1, a2, lambda1, l_u_1_bank, l_I_1_bank, Delta1);
    F(2) = FOC_shadow(a2,a1,lambda2);
elseif v_bank(a1, a2, l_u_1_bank, l_I_1_bank, lambda1, Delta1)< v_shadow(a1, a2, l_u_1_shadow, lambda1) && v_bank(a2, a1, l_u_2_bank, l_I_2_bank, lambda2, Delta2)> v_shadow(a2, a1, l_u_2_shadow, lambda2)
    F(1) = FOC_shadow(a1, a2, lambda1);
    F(2) = FOC_bank(a2, a1, lambda2, l_u_2_bank, l_I_2_bank, Delta2);
else
    F(1) = FOC_shadow(a1, a2, lambda1);
    F(2) = FOC_shadow(a2, a1, lambda2);
end




end