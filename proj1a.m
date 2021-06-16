clear;
clc;

fprintf('\n############################################################\n')
fprintf(' Levenberg-Marquardt Non-Linear Least Squares Program ')
fprintf('\n############################################################\n')
fprintf(' Kayode Olumoyin ')
fprintf('\n############################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input the filename of the textfile containing pressure and volume data
% using the virial equation of state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Some matices in this project may be close to singular,
% these matrices may be badly scaled and badly conditioned.
% Their RCOND is close to zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long e % using the highest matlab precision

warning('off','all') % This suppresses the RCOND warnings
prompt = 'Enter the input text file: ';

% Enter any of the following txt files %

% Ar_223.15K_virial_experiment_B2-B6.txt
% Ar_273.15K_virial_theory_B2-B7.txt
% He_148.15K_virial_B2-B4.txt
% He_50K_virial_B2-B5.txt
% air_100K_virial_B2.txt
% ethane_215K_virial_B2-B3.txt

filename = input(prompt,'s');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID = fopen(filename,'r');

header = textscan(fileID,' %s %s %s ',1,'Delimiter',',');
%celldisp(header(1,3))
fprintf('\n')
fprintf('comment = %s\n', string(header(1,3)))


fclose(fileID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab function to detect the input file as a table 
opts = detectImportOptions(filename); 

% set table as a table of double type
file_double = readtable(filename,opts);

% check the variables detected and the datatype detected
% disp([opts.VariableNames' opts.VariableTypes'])

% set table as a table of string type
opts = setvartype(opts,'string');
file_string = readtable(filename,opts);

%preview(filename,opts)

fprintf('\n############################################################\n')

% detect the equation of state from the input file
eos = file_string(1,1);
eos_v = table2array(eos(1,1));
fprintf('\n Equation of state = %s\n', eos_v)

% detect the number of parameters from input file
num_par = file_double(1,2);
num_param = table2array(num_par(1,1));
fprintf('\n Number of parameters = %i\n',num_param)

% detect the initial parameter guesses from input file
guess1 = file_double(2,1:num_param);
guess2 = table2array(guess1);
for k=1:num_param
    fprintf('\n Initial guess, parameter %i = %0.6e',k,guess2(k))
end
fprintf('\n')

% detect temperature coefficient from input file
temp1 = file_double(3,2);
Temp1 = table2array(temp1(1,1));

% detect temperature unit from input file
temp2 = file_string(3,3);
Temp2 = table2array(temp2(1,1));

fprintf('\n Temperature = %i %s \n', Temp1, Temp2)

% detect volume unit from input file
vol_uni = file_string(4,1); 
vol_unit = table2array(vol_uni(1,1));
fprintf('\n Volume unit in input text file = %s \n', vol_unit)

% detect pressure unit from input file
pre_uni = file_string(4,2);
pre_unit = table2array(pre_uni(1,1));
fprintf('\n Pressure unit in input text file = %s \n', pre_unit)

% detect the volume data
vol = file_double(5:end,1);
Vm = table2array(vol(:,1));

% detect the pressure data
pre = file_double(5:end,2);
P = table2array(pre(:,1));

% detect the number of volume and pressure data
data_len = length(Vm);
fprintf('\n Number of input text file = %i \n', data_len)

fprintf('################################################################')

% print first 10 volume data
volume = file_double(5:14,1);
Vm_10 = table2array(volume(:,1));

fprintf('\n first 10 volume element = \n')
disp(Vm_10)

% print first 10 pressure data
pressure = file_double(5:14,2);
P_10 = table2array(pressure(:,1));
fprintf('\n first 10 pressure element = \n')
disp(P_10)

fprintf('################################################################')

% convert Vm to m^3/mol
% new_Vm = zeros(length(Vm));
vol_unit = lower(vol_unit);
if vol_unit == 'm^3/mol'
    new_Vm = Vm;
elseif vol_unit == 'dm^3/mol'
    new_Vm = 0.001*Vm;
elseif vol_unit == 'l/mol'
    new_Vm = 0.001*Vm;
elseif vol_unit == 'cm^3/mol'
    new_Vm = 0.000001*Vm;
end

% convert P to Pa
% new_P = zeros(length(P));
pre_unit = lower(pre_unit);
if pre_unit == 'pa'
    new_P = P;
elseif pre_unit == 'mmhg'
    new_P = 133.322387415*P;
elseif pre_unit == 'atm'
    new_P = 101325*P;
elseif pre_unit == 'torr'
    new_P = 133.3223684211*P;
elseif pre_unit == 'megapa'
    new_P = 1000000*P;
elseif pre_unit == 'kilobar'
    new_P = 100000000*P;
elseif pre_unit == 'bar'
    new_P = 100000*P;
end

% print first 10 SI volume data, m^3/mol
new_Vm_10 = new_Vm(1:10);
fprintf('\n first 10 SI volume data, m^3/mol = \n')
disp(new_Vm_10)

% print first 10 SI pressure data, Pa
new_P_10 = new_P(1:10);
fprintf('\n first 10 SI pressure data, Pa = \n')
disp(new_P_10)

clear opts
fprintf('\n############################################################\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the gas constant and temperature coefficient in the input file
R = 8.31447;                % Pa m^3 K^(-1) mol^(-1)
T = Temp1;                  % K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calling the levenberg-marquardt routine
% lambda (big) correspond to a steepest descent minimizer
% lambda (small) correspond to a Gauss-Newton minimizer

old_param = guess2;
la = 10000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n############################################################\n')
fprintf('\n CYCLE %i, lambda = %0.20e \n', 1, la)

for k=1:num_param
    fprintf('\n old parameter %i = %0.6e',k,old_param(k))
end
fprintf('\n')

score1 = score(new_Vm, new_P, old_param, R, T);
fprintf('\n The old score = %0.20e\n', score1)

bbeta = beta(new_Vm, new_P, old_param, R, T);
fprintf('\n The beta vector = \n')
disp(bbeta)

aalpha = myHessian(new_Vm, old_param, R, T);
fprintf('\n The hessian matrix, alpha = \n')
disp(aalpha)

aalpha_prime = lambda_myHessian(new_Vm, old_param, R, T, la);
fprintf('\n The updated hessian matrix, alpha_prime = \n')
disp(aalpha_prime)

%delta_a = linsolve(aalpha_prime,bbeta);
g = myCholesky(aalpha_prime);
y = LowerTrig(g,bbeta);
delta_a = UpperTrig(g',y);
%delta_a = delta_a';
fprintf('\n The delta_param = \n') 
disp(delta_a)

new_param = old_param + delta_a';
%old_param' + delta_a;
fprintf('\n The new parameters = \n')
disp(new_param')

score2 = score(new_Vm, new_P, new_param, R, T);
fprintf('\n The new score = %0.20e\n', score2)
% fprintf('\n##########################################################\n')

counter = 2;
for i=1:100
    if score2 >= score1
        new_la = la * 10;
        aalpha_prime = lambda_myHessian(new_Vm, old_param, R, T, new_la);
        new_delta_a = linsolve(aalpha_prime,bbeta);
        %new_delta_a = new_delta_a';
        nnew_param = new_param + new_delta_a';
        score2 = score(new_Vm, new_P, nnew_param, R, T);
    else
        new_la = la / 10;
        old_param = new_param;
        
        fprintf('\n')
        fprintf('########################################################')
        fprintf('\n CYCLE %i, lambda = %0.20e \n', counter, new_la)
        
        score1 = score2;
        la = new_la;
        
        bbeta = beta(new_Vm, new_P, old_param, R, T);
        fprintf('\n The beta vector = \n')
        disp(bbeta)

        aalpha = myHessian(new_Vm, old_param, R, T);
        fprintf('\n The hessian matrix, alpha = \n')
        disp(aalpha)

        aalpha_prime = lambda_myHessian(new_Vm, old_param, R, T, la);
        fprintf('\n The updated hessian matrix, alpha_prime = \n')
        disp(aalpha_prime)

        %delta_a = linsolve(aalpha_prime,bbeta);
        g = myCholesky(aalpha_prime);
        y = LowerTrig(g,bbeta);
        delta_a = UpperTrig(g',y);
        %delta_a = delta_a';
        fprintf('\n The delta_param = \n') 
        disp(delta_a)

        new_param = old_param + delta_a';
        %old_param' + delta_a;
        fprintf('\n The new parameters = \n')
        disp(new_param')
        
        fprintf('\n The old score = %0.20e\n', score1)

        score2 = score(new_Vm, new_P, new_param, R, T);
        fprintf('\n The new score = %0.20e\n', score2)
        fprintf('\n########################################################\n')
    end
    counter = counter + 1;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(' ')
fprintf('#################### final statistics ##########################')

fprintf('\n chi square = %0.20e\n', score2)

var = myVariance(new_Vm, new_P, new_param, R, T);
fprintf('\n sample variance = %0.20e\n', var)

fprintf('\n alpha = \n')
disp(aalpha_prime)

Id = eye(length(aalpha_prime)); 
% using cholesky function, see C.F. Van Loan, Intro Sci Comput, pg 261
% %tf = issymmetric(aalpha_prime); % check if A is a symmetric positive definite matrix
g1 = myCholesky(aalpha_prime);
y1 = linsolve(g1,Id);
C = linsolve(g1',y1);
fprintf('\n C = \n')
disp(C)

fprintf('\n The final parameters = \n')
disp(new_param')

se = sqrt(var*C);
fprintf('\n The standard error = \n')
disp(diag(se))

fprintf('\n parameter correlation coefficient, rho = \n')
rho = zeros(num_param);
[rows, cols] = size(se);
for row=1:rows
    for col=1:cols
        if row < col
            rho(row,col) = se(row,col)^2/(se(row,row)*se(col,col));
            fprintf('rho(%i, %i) = %0.20e\n',row, col, rho(row,col));
        end
    end
end

pmean = myMean(new_P);
fprintf('\n pmean = %0.20e\n',pmean)

ssT = ssTotal(new_P);
fprintf('\n SS_total = %0.20e\n',ssT)

RR = 1 - (score2/ssT);
fprintf('\n The coefficient of determination, r^2 = %0.20e\n', RR);

fprintf('\n The correlation coefficient, r = %0.20e\n', sqrt(RR));

ad_R = 1 - ((score2/(length(new_P) - num_param - 1))/(ssT/(length(new_P) - 1)));
fprintf('\n adjusted r^2 = %0.20e\n', ad_R)

r_fac = rFactor(new_Vm, new_P, new_param, R, T);
fprintf('\n R-factor = %0.20e\n', r_fac)

fprintf('\n R-factor in percentage = %0.20e\n', r_fac*100)
fprintf('\n############################################################\n')    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plotting
Volume_fit = new_Vm;
Pressure_fit = f(Volume_fit, new_param, R, T);

figure
plot(Volume_fit,Pressure_fit,'r',new_Vm(1:100:end), new_P(1:100:end),'bo','markersize',4)
xlabel m^3/mol, ylabel Pa
legend('fitted model','data')
title(string(header(1,1)), string(header(1,2)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% eos
function pv = f(x, par, R, T)
new_B = zeros(length(x),1);
num_param = length(par);
for i=1:num_param
    old_B = (R*T*par(i))*(1./x.^(i+1));
    new_B = old_B + new_B;
end    
pv = (R*T)*(1./x) + new_B;
end

% jacobian
function jacc = jac_f(x, par, R, T)
num_param = length(par);
jacc = zeros(length(x),num_param);
for i=1:num_param
    jacc(:,i) = (R*T)*(1./x.^(i+1));
end
end

% score
function loss = score(x, y, par, R, T)
y_calc = f(x, par, R, T);
delta = y - y_calc;
new_delta = delta.^2;
loss = sum(new_delta);
end

% beta
function b = beta(x, y, par, R, T)
y_calc = f(x, par, R, T);
J = jac_f(x, par, R, T);
delta = y - y_calc;
num_param = length(par);
b = zeros(num_param,1);
for i=1:num_param
    b(i) = dot(delta,J(:,i));
end
end

% hessian
function hes = myHessian(x, par, R, T)
J = jac_f(x, par, R, T);
hes = J'*J;
end

% hessian_prime
function hes_prime = lambda_myHessian(x, par, R, T, la)
hes = myHessian(x, par, R, T);
[rows,cols] = size(hes);
for row=1:rows
    for col=1:cols
        if row == col
            hes_prime(row,col) = hes(row,col)*(1+la);
        else
            hes_prime(row,col) = hes(row,col);
        end
    end
end
end

% variance
function v = myVariance(x, y, par, R, T)
num_param = length(par);
K = (y - f(x, par, R, T)).^2;
Sum = sum(K);
v = (1/(length(x)-num_param))*Sum;
end

% mean(y)
function ybar = myMean(y)
K = sum(y);
ybar = (1/length(y))*K;
end

% ssTotal(y)
function sT = ssTotal(y)
K = (y - myMean(y)).^2;
sT = sum(K);
end

% rFactor
function r_fac = rFactor(x, y, par, R, T)
K1 = abs(y - f(x, par, R, T));
Sum1 = sum(K1);
K2 = abs(y);
Sum2 = sum(K2);
r_fac = Sum1/Sum2;
end

% cholesky factorization
function G = myCholesky(A)
[n,n] = size(A); G = zeros(n,n);
for i=1:n
    for j=1:i
        if j==1
            s=A(j,i);
        else
            s=A(j,i)-G(j,1:j-1)*G(i,1:j-1)';
        end
        if j<i
            G(i,j)=s/G(j,j);
        else
            G(i,i)=sqrt(s);
        end
    end
end
end

% forward substitution
function x = LowerTrig(L,b)
n = length(b);
x = zeros(n,1);
for j = 1:n-1
    x(j) = b(j)/L(j,j);
    b(j+1:n) = b(j+1:n) - L(j+1:n,j)*x(j);
end
x(n) = b(n)/L(n,n);
end

% backward substitution
function x = UpperTrig(U,b)
n = length(b);
x = zeros(n,1);
for j = n:-1:2
    x(j) = b(j)/U(j,j);
    b(1:j-1) = b(1:j-1) - x(j)*U(1:j-1,j);
end
x(1) = b(1)/U(1,1);
end
