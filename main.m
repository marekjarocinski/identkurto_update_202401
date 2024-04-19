% Estimate Y = UC, Student-t, unknown v
% where
% Y - observable data, i.i.d.
% U - i.i.d. Student-t WITH UNKNOWN v>1
% C - matrix with the impacts
% We will actually estimate W in YW = U, then recover C = inv(W)
%
% Maximum Likelihood
% Depends on: nlogl_iidstudent_Wzn
% Marek Jarocinski

clear all, close all

% read data
ttab = readtimetable("fomc_surprises_jk.csv",Delimiter=",");
ttab.start.Format = "uuuu-MM-dd HH:mm";
ttab(year(ttab.start)<1991,:) = [];
fprintf('Data from %s to %s, T=%d\n', ttab.start(1), ttab.start(end),size(ttab,1))

ynames = ["MP1","TFUT02","TFUT10","SP500"]

% drop missing values
indmissing = logical(sum(isnan(ttab{:,ynames}),2));
disp("dropping the following observations due to missing values:")
disp(ttab(indmissing,ynames))
ttab(indmissing,:) = [];
fprintf('Data from %s to %s, T=%d\n', ttab.start(1), ttab.start(end),size(ttab,1))

Y = ttab{:,ynames}*100;
ymaturities = [1/12 2 10 NaN];
[T, N] = size(Y);
nv = N; % 1 or N

figure()
for n = 1:N
subplot(2,2,n)
plot(ttab.start, ttab{:,ynames(n)})
title(ynames(n))
end


% Maximum likelihood / posterior mode estimation of W(->C), v
options = optimoptions('fmincon');
options.Display = 'iter';
options.MaxFunctionEvaluations = 1e4*N^2;
options.OptimalityTolerance = 1e-9;
options.Algorithm = 'trust-region-reflective';
options.SpecifyObjectiveGradient = true;

% constraints
A = []; b = []; Aeq = []; beq = []; lb = [repmat(-Inf,N^2,1); repmat(0,N,1)]; ub = []; nonlcon = [];

% starting point
W = inv(chol(cov(Y)));
par0 = [W(:); repmat(log(3),N,1)];

% maximum likelihood
fun = @(par) nlogl_iidstudent_Wzn(Y, reshape(par(1:N^2),N,N), par(N^2+1:end));
%[parmaxlik,fval,exitflag,output,grad,hessian] = fminunc(fun, par0, options);
[parmaxlik,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,par0,A,b,Aeq,beq,lb,ub,nonlcon,options);
W = reshape(parmaxlik(1:N^2), N, N);
C = inv(W);
v = exp(parmaxlik(N^2+1:end));

% Normalize signs and order
% 1. Sign: positive interest rate shock
P = diag(sign(mean(C(:,~isnan(ymaturities)),2)));
% 2. Order by the maximum impact on consecutive variables
idx = nan(1,N);
temp = abs(C);
for n = 1:N
    [~, i] = max(temp(:,n));
    idx(n) = i;
    temp(i,:) = -inf;
end
P = P(idx,:);

%P*C
Wnormalized = W*P';
vnormalized = abs(P)*v;

if 0 % manually fix the ordering
    C = inv(Wnormalized);
    if C(3,4)>C(4,4)
        P = eye(N);
        P = P(:,[1 2 4 3]);
    end
    Wnormalized = Wnormalized*P';
    vnormalized = abs(P)*vnormalized;
end
disp('P'), disp(P)

% Minimize again to get the reordered hessian
par0 = [Wnormalized(:); log(vnormalized)];
[parmaxlik,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,par0,A,b,Aeq,beq,lb,ub,nonlcon,options);

maxl.Y = Y;
maxl.ll = -fval;
maxl.parmaxlik = parmaxlik;
maxl.hessian = hessian;
maxl.v = exp(parmaxlik(N^2+1:end));
maxl.W = reshape(parmaxlik(1:N^2), N, N);
maxl.C = inv(maxl.W);
maxl.U = Y*maxl.W;
maxl.U1s = maxl.U./std(maxl.U);
maxl.C1s = diag(std(maxl.U))*maxl.C;
maxl.vdec = (maxl.C1s.^2)./sum(maxl.C1s.^2);
if N==4, temp = [1 2 3 2]; else, temp = 1:N; end
temp = diag(diag(maxl.C(:,temp)));
maxl.C1bp = temp\maxl.C;
maxl.U1bp = maxl.U*temp;

disp('C'), disp(maxl.C)
disp('C1s'), disp(maxl.C1s)
disp('C1bp'), disp(maxl.C1bp)
disp('vdec'), disp(maxl.vdec)
disp('v'), disp(maxl.v)

% save maxl and shocks
save('maxl','maxl');
%y = round(maxl.Y, 5);
%writetimetable(array2timetable(y, "RowTimes", ttab.start, "VariableNames", ynames), "data.csv");
%u = round(maxl.U, 5);
%writetimetable(array2timetable(u, "RowTimes", ttab.start), "U.csv");
u = round(maxl.U1s, 5);
writetimetable(array2timetable(u, "RowTimes", ttab.start), "U1s.csv");
u = round(maxl.U1bp, 5);
writetimetable(array2timetable(u, "RowTimes", ttab.start), "U1bp.csv");

% Variances
JacobianWC = -kron(maxl.W',maxl.W);
Jacobianzv = diag(maxl.v.^-1);
JacobianWCzv = blkdiag(JacobianWC, Jacobianzv);
asyvarCv = inv(JacobianWCzv'*hessian*JacobianWCzv);
asystdCv = sqrt(diag(asyvarCv));
Cstd = reshape(asystdCv(1:N^2),N,N);
disp([maxl.v asystdCv(N^2+1:end)])

Cm = maxl.C;
Cl = Cm - 2*Cstd;
Cu = Cm + 2*Cstd;

fh = plot_resp(diag(std(maxl.U))*Cm, diag(std(maxl.U))*Cl, diag(std(maxl.U))*Cu, ymaturities, ynames);
exportgraphics(fh, 'C1smaxlik_band.pdf')
fh = plot_resp(temp\Cm, temp\Cl, temp\Cu, ymaturities, ynames);
exportgraphics(fh, 'C1bpmaxlik_band.pdf')




