%Let Delta and Lambda be independent.
%Aim1: Make it Stable? : Almost stable
%Aim2: Shadow banks with higher mu and std?
%Intuitively: why shadow banks have higher mu -- selected by higher
%entrycost
%correlated with p, gamma_b, gamma_s

%Solve problem of first error: set two entrants initially

clear;
clc;
warning on;
%% setting paras
tol_val = 1e-25;
options_fsolve = optimoptions('fsolve','Display','off','MaxIter',100000,'MaxFunctionEvaluations', 100000,'OptimalityTolerance',tol_val,'FunctionTolerance',tol_val,'StepTolerance',tol_val);


% p = 0.632;
p = 0.777;  % probability of good state                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
%gamma_s = 0.0089;
gamma_s = 0.0029; 
% gamma_b = 0.0139;
gamma_b = 0.013;
delta = 0.4;
L_bar = 7;

%beta = 0.97;n
beta = 0.4;
xi = 0.01;

mu_Delta = 0.5;
% sigma_Delta = 0.003;
sigma_Delta = 0.03;
mu_lambda = 0.9;
% sigma_lambda = 0.01;
sigma_lambda = 0.3;

% EntryCost_shadow = 1e-6;
% EntryCost_bank = 6e-6;
EntryCost_shadow = 1e-6;
EntryCost_bank = 4e-7;

A = 0; 

maxit = 5000;

lambda_set = [];
choice_set = [];
choice = -1;

a1 = -1;
a2 = -1;
%% Two initial entrants: 
while a1 < 0 || a2 < 0
%Random Productivity types: ensure positive productivity types
lambda1 = -1;  
while lambda1 <= 0
    lambda1 = normrnd(mu_lambda, sigma_lambda);
end
lambda2 = -1;  
while lambda2 <= 0
    lambda2 = normrnd(mu_lambda, sigma_lambda);
end

% lambda1 = normrnd(mu_lamdba, sigma_lambda);
% lambda2 = normrnd(mu_lamdba, sigma_lambda);

%Random deposit risks

Delta1 = normrnd(mu_Delta , sigma_Delta); 
Delta2 = normrnd(mu_Delta , sigma_Delta); 

%solve uninsured leverage of shadow banks in closed form
L_u_shadow_1 = gamma_s/(2*(1-p)*delta) + lambda1;
L_u_shadow_2 = gamma_s/(2*(1-p)*delta) + lambda2;
%solve uninsured leverage and insured leverage of banks in closed form
L_u_banks_1 = (Delta1 - p)/(2*p*delta) + lambda1;
L_i_banks_1 = L_bar - L_u_banks_1;
L_u_banks_2 = (Delta2 - p)/(2*p*delta) + lambda2;
L_i_banks_2 = L_bar - L_u_banks_2;

%solve asset size a_1 and a_2 by solving equations under different
%conditions
ai = ones(2,1);
[ai, errors] = fsolve(@(a)Initial_Two_Entrants_May10(a, beta, xi, p, L_bar, gamma_b, delta, lambda1, lambda2, Delta1, Delta2, L_u_shadow_1, L_u_shadow_2, L_u_banks_1, L_u_banks_2, L_i_banks_1, L_i_banks_2,gamma_s),ai,options_fsolve);
max(abs(errors))
[a1, a2] = deal(ai(1), ai(2));
A = a1^beta + a2^beta;
end

%% Other entrants

for i= 1: maxit-2
% lambda = normrnd(mu_lamdba, sigma_lambda);
lambda = -1;  
while lambda <= 0
    lambda = normrnd(mu_lambda, sigma_lambda);
end
lambda_set(end+1) = lambda;

% while Delta < p-lambda*2*p*delta % ensure that uninsured leverage of banks  >0
Delta = normrnd(mu_Delta , sigma_Delta);  
% end

%For a new A, define a new production function and derivative
f = @(ai)xi*ai^beta/(A + ai^beta);
f_prime = @(ai)xi*(beta*ai^(beta - 1)/(A+ai^beta) - (beta*ai^(2*beta - 1))/((A + ai^beta)^2));

%% Value function of Shadow Banks
%solve uninsured leverage of shadow banks in closed form
L_u_shadow = gamma_s/(2*(1-p)*delta) + lambda;

%asset size
FOC_shadow = @(ai)f_prime(ai)*(p+(1-p)*lambda) - 1 + lambda*gamma_s + (gamma_s)^2/(4*(1-p)*delta);

ai = 1;
[ai_shadow,errors] = fsolve(FOC_shadow,ai,options_fsolve);
max(abs(errors))

%value function of becoming a shadow bank
v_shadow = @(ai, l_u)f(ai)*(p+(1-p)*lambda)-ai + ai*l_u*gamma_s - (1-p)*ai*delta*(l_u - lambda)^2;

%% Value function of Banks
%solve uninsured leverage and insured leverage of banks in closed form
%Note: leverage of banks: < L_bar
L_u_banks = (Delta - p)/(2*p*delta) + lambda;
L_i_banks = L_bar - L_u_banks;

%Asset size
FOC_bank = @(ai)-f_prime(ai)*(p+(1-p)*lambda) + 1 - gamma_b*L_bar + (1-p)*delta*(L_u_banks - lambda)^2 - L_i_banks*(1-p) + L_i_banks*Delta;
ai = 1;
[ai_banks,errors] = fsolve(FOC_bank,ai,options_fsolve);
max(abs(errors))

%value function of becoming a bank
v_bank = @(ai, l_u, l_i)f(ai)*(p+(1-p)*lambda) - ai + ai*L_bar*gamma_b...
                       - (1-p)*ai*delta*(l_u - lambda)^2  - ai*l_i*Delta...
                       + ai*l_i*(1-p);

%% Choice of becoming a bank or shadow bank
if v_bank(ai_banks, L_u_banks, L_i_banks)- EntryCost_bank > v_shadow(ai_shadow, L_u_shadow)- EntryCost_shadow && v_bank(ai_banks, L_u_banks, L_i_banks)- EntryCost_bank >0
    choice = 1;
    ai_final = ai_banks;
elseif v_bank(ai_banks, L_u_banks, L_i_banks)- EntryCost_bank < v_shadow(ai_shadow, L_u_shadow)- EntryCost_shadow && v_shadow(ai_shadow, L_u_shadow)- EntryCost_shadow >0
    choice = 0;
    ai_final = ai_shadow;
else
    choice = 2;   
    ai_final = 0;
end


choice_set(end+1)=choice;

% Since the situation might happen when the former small lambda can not
% enter but the latter one can, we do not break when one lambda can not
% enter. Instead, we stop updating A under that situation.
if choice == 1 || choice == 0 
A = A+ai_final^beta;
end
end
%% Plot
lambda_shadow = lambda_set(choice_set == 0);
lambda_bank   = lambda_set(choice_set == 1);

figure;
hold on;

histogram(lambda_shadow, 'Normalization', 'probability', ...
    'FaceColor', [0.2 0.4 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

histogram(lambda_bank, 'Normalization', 'probability', ...
    'FaceColor', [0.8 0.3 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

% histogram(lambda_shadow, 'Normalization', 'count', ...
%     'FaceColor', [0.2 0.4 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
% histogram(lambda_bank, 'Normalization', 'count', ...
%     'FaceColor', [0.8 0.3 0.3], 'FaceAlpha', 0.5, 'EdgeColor', 'none');

xlabel('Productivity (\lambda)');
ylabel('Frequency');
title('Distribution of Productivity by Bank Type');
legend({'Shadow Banks', 'Banks'});
grid on;
hold off;

%%
v1 = v_bank(ai_banks, L_u_banks, L_i_banks);
v2 = v_shadow(ai_shadow, L_u_shadow);















