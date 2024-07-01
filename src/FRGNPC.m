function [VGg,Out] = FRGNPC(Gamma,gamma,opt)
%function [Zexp,VGg,flag,Out] = FRGNPC(Gamma,gamma,opt)
%%% This code solves the FR problem for Gamma rho = gamma, rho psd
%%% Now runs predictor-corrector
%%% Input:   Gamma: cell array holding the constraint data matrices
%%%          gamma : RHS vector for Gamma(rho) = gamma
%%%          opt : struct; option for the algorithm
%%%                tolerFR : stopping tolerance
%%%                SlaterFR : Slater tolerance; mineig(W)>SlaterFR
%%%                iterbnd : max number of iterations
%%%                P : a matrix preventing A*(y) ~= 0; <A*(y),P> = rhspos
%%%                rhspos : positive constant for <A*(y),P> = rhspos
%%%                verbose : print out auxiliary outputs/algorithm progress
%%% Output: Zexp : exposing vector  Hsvec(Zexp) = MGamma y psd,
%%%         VGg : FR vector that yields VRV^T
%%%         flag :(1,0,-1,-2) = ...
%%%               (problem reduced,Slater holds,problem infeas,inconclusve)
%%%
%%%
%%% (primal) min_y {b'*y : A^*(y)>=0, <P,A^*(y)> = rhspos }
%%% (dual)   max_{lambda,W} { rhspos*lambda : A(W)+lambda*A(P)=b, W>= 0 }
%%% There are four blocks in the optimality conditions
%%% dual feasibility  A(W)+lambda*A(P)=b is written with
%%%                      W = What + NW*v  AND
%%%                      lambda = lambdahat + Nlambda*v
%%% primal feasibility <P,A^*(y)> = 1 is written with
%%%                     y = yhat + Ny*u
%%% complementarity  W*A^*(y) - mu*In


%%% NOTE : The code handles both Hermitian and the symmetric variables;
%%%        controlled by varflag


if isfield(opt,'tolerFR'), tolerFR = opt.tolerFR;  else, tolerFR = 1e-12;  end
if isfield(opt,'Slatertoler'), Slatertoler = opt.Slatertoler;  else, Slatertoler = 1e-13;  end
if isfield(opt,'iterbnd'), iterbnd = opt.iterbnd;  else, iterbnd = 5e1;  end
if isfield(opt,'P'), P = opt.P;  else, P = eye(length(Gamma{1}));  end
if isfield(opt,'rhspos'), rhspos = opt.rhspos;
else rhspos = max(floor(length(Gamma{1})/2),1); end
if isfield(opt,'verbose'), verbose = opt.verbose;  else, verbose = 1;  end

if verbose >=2
    fprintf(2,'<starting FRGNPC.m facial reduction> \n')
end

n = length(Gamma{1});
n2 = n^2;
m = length(gamma);   % number of constraints
In = eye(n);
tn = n*(n+1)/2;      % triangular number

if verbose >= 1
    fprintf('dim n = %d, # constraints (m) = %d \n',n,m)
end

% Idendify the space; real or complex
if sum(cellfun(@isreal, Gamma)) == m
    varflag = 'realCase';
    if verbose >=1; fprintf('variable: real symmetric\n'); end
else
    varflag = 'complexCase';
    if verbose >= 1; fprintf('variable: Hermitian\n'); end
end

% Obtain the matrix representation of Gamma()
%  Gamma(rho) = MGamma*HSvec(rho) = gamma
if strcmp(varflag,'realCase')
    MGamma = zeros(m,n*(n+1)/2);
    for ii = 1:length(gamma)
        MGamma(ii,:) = HSvec(Gamma{ii},'realCase');
    end
elseif strcmp(varflag,'complexCase')
    MGamma = zeros(m,n^2);
    for ii = 1:length(gamma)
        MGamma(ii,:) = HSvec(Gamma{ii},'complexCase');
    end
end

if strcmp(varflag,'realCase') && size(MGamma,2) > tn  % if there are too many columns in the real case
    tn1 = n*(n-1)/2;
    MGamma(:,tn1+1:2*tn1) = []; % throw away the complex part
end

% set isomorphic dimension of W
if strcmp(varflag,'realCase')
    isodim = tn;  % dim(W) = tn, when W real symmetric
elseif strcmp(varflag,'complexCase')
    isodim = n2; % dim(W) = n^2, when W Hermitian
end


% preprocessing of the data
MGammaorig = MGamma;
Gammaorig = Gamma;
gammaorig = gamma;


eigGamma = cellfun(@eig,Gamma,'UniformOutput',false);
eigGammamin = cellfun(@min,eigGamma);
eigGammamax = cellfun(@max,eigGamma);
psdind = find(eigGammamin>-eps); % psd data matrices
nsdind = find(eigGammamax<eps);  % nsd data matrices


% Detect infeasibility by looking at the eigenvalues of the data
if any(gamma(nsdind)>=0) && sum(abs(gamma))~=0
    fprintf('The problem is infeasible; see %d-th constraint.\n',nsdind(1))
    fprintf('max eigenvalue is %e and the rhs is %d\n',eigGammamax(nsdind(1)),gamma(nsdind(1)))
    keyboard
    return
end
if any(gamma(psdind)<0)
    fprintf('The problem is infeasible; see %d-th constraint.\n',psdind(1))
    fprintf('min eigenvalue is %e and the rhs is %d\n',eigGammamin(psdind(1)),gamma(psdind(1)))
    keyboard
    return
end


VarNames = {'-log10(mingap)','stepDF','stepPF','sigmaa','objval','pri-feas','dualfeas','cond.Jac','relres(pred)','relres(corr)'};
if verbose >= 2
    fprintf('%-9s','iter#')
    fprintf('%-15s',VarNames{:})
    fprintf('\n')
end

%% Initialization
bestub = inf;
bestlb = -inf;
ubs(1) = inf;
lbs(1) = -inf;

flagsteponeDF = 0;  % number of times step length one taken
flagsteponePF = 0;  % number of times step length one taken

sigmaa = 1;     % centering parameter - changes adaptively
stepPF = 1.2;     % initial primal stepsize - changes adaptively
stepDF = 1.2;     % initial dual stepsize - changes adaptively


tg = zeros(1,m);
for ii = 1:m
    tg(ii) = trace(Gamma{ii});
end

%P = In; % see Latex; better P can be chosed to guarantee strong duality
GI = MGamma*HSvec(P,varflag);


JconstD = [MGamma GI]; % upper constant top part (in writeup) of Jacobian; dual
JconstP = [tg, zeros(1,isodim) ;
    MGamma', -eye(isodim)]; % lower constant part of the Jacobian; primal

%% nullsp. repre dual feasiblity [MGamma, A(I)]*[W;lambda]=gamma
%%% obtain sparse nullspace repre. for JconstD.
rp = 1:m;  % row permutation
[~,cp] = licols(JconstD);   % obtain lin.indep columns of JconstD
if length(cp)~=length(rp)
    fprintf('Error: Check the number of rows and columns in permutation.\n')
    keyboard
end
ccp = setdiff(1:isodim+1,cp);   % complement indices
cp = [cp ccp];              % column permutation
[~,icpD] = sort(cp);        % inverse permutation of cols
JconstD_perm = JconstD(rp,cp); % order of rows/constraints does not change problem
gamma_perm = gamma(rp);  % for RHS

lambdaind = find(cp==isodim+1); % where lambda ind is located after permutation

BMG = sparse(JconstD_perm(:,1:m));
EMG = sparse(JconstD_perm(:,m+1:end));
ND = [-sparse(BMG\EMG); speye(isodim+1-m)];  % sparse nullsp repre for dual feas

JconstD_noperm = JconstD_perm(:,icpD);  % MGammas before column permutations, but with row permutations
ND_noperm = ND(icpD,:); % matrix holding null basis for MGamma_noperm


%% nullsp. repre primal feasiblity trace(A*(y))=rhspos and A^*y-S=0
%%% obtain sparse nullspace repre. for JconstF.
rp = 1:isodim+1;    % row permutation
%[~,indtemp] = max(abs(tg)); % obtain lin.indep columns of JconstP
%indtemp = [indtemp, 1+isodim:m+isodim]; % take advatage of the identity block
[~,cp] = licols(JconstP);   % obtain lin.indep columns of JconstD
if length(cp)~=length(rp)
    keyboard
end
ccp = setdiff(1:isodim+m,cp);   % complement indices
cp = [cp ccp];
[~,icpP] = sort(cp);   % inverse permutation of cols
JconstP_perm = JconstP(rp,cp); % order of rows/constraints does not change problem
BMG = JconstP_perm(:,1:isodim+1);
EMG = JconstP_perm(:,isodim+2:end);
NP = [ -BMG\EMG; speye(m-1) ];  % sparse nullsp repre for primal feas  % #basis elt=m-1
JconstP_noperm = JconstP_perm(:,icpP);  % MGammas before column permutations, but with row permutations
NP_noperm = NP(icpP,:); % matrix holding null basis for MGamma_noperm


% k = size(ND,2) + size(NP,2); % n^2
%% assign NW,Nlambda,Ny  (null basis for W,lambda,y); permutated null-base
% null basis after 'undoing' the permutations
NW_noperm = ND_noperm(1:isodim,:);
Nlambda_noperm = ND_noperm(end,:);
Ny_noperm = NP_noperm(1:m,:);
NS_noperm = NP_noperm(m+1:end,:);


MGamma_perm = JconstD_perm; % MGamma after permutaton
MGamma_perm(:,lambdaind) = []; % remove the lambda column

%% obtain particular solutions (What,lambdahat,yhat)
JacsolD = lsqminnorm(JconstD,gamma);
JacsolP = lsqminnorm(JconstP,[rhspos;zeros(isodim,1)]);

What = HSMat(JacsolD(1:end-1),varflag);
lambdahat = JacsolD(end);
yhat = JacsolP(1:m);
Shat = HSMat(JacsolP(m+1:end),varflag);


%% set initial iterates (v,u,W,lambda,y,S)
v = zeros(size(ND,2),1);
u = zeros(size(NP,2),1);

d = min([eig(What);0.1]);
if d < .1   % ensure pos def
    W = What + (abs(d)+.1)*In;  % ensure pos def start
end
lambda = 0;
y = (MGamma')\HSvec(P,varflag);
S = HSMat(MGamma'*y,varflag);
d = min([eig(S);0.1]);
if d < .1   % ensure pos def
    S = S + (abs(d)+.1)*In;  % ensure pos def start
end
y = (MGamma')\HSvec(S,varflag);

objval = gamma'*y;

gapp = real(trace(W*S));
mingapp = gapp;
muu = gapp/(n);


Fmu{1} = HSMat(NW_noperm*v,varflag) + What - W;
Fmu{2} = Nlambda_noperm*v + lambdahat - lambda;
Fmu{3} = Ny_noperm*u + yhat - y;
Fmu{4} = HSMat(NS_noperm*u,varflag) + Shat - S;
Fmu{5} = W*S - muu*In;


Fmuprvec = [HSvec(Fmu{1},varflag);Fmu{2}; ...
    Fmu{3};HSvec(Fmu{4},varflag);CRvec(Fmu{5},varflag)];
Fmuprnorm = norm(Fmuprvec);

noiter = 0;         % iteration counter
condJac = 0;        % condition number
stalling = 0;       % stalling to quit while loop
earlystop = false;  % early stopping flag

if strcmp(varflag,'realCase')
    JacD = zeros(n2,isodim-m+1);
    JacP = zeros(n2,m-1);
elseif strcmp(varflag,'complexCase')
    JacD = zeros(2*n2,isodim-m+1);
    JacP = zeros(2*n2,m-1);
end
%%%%%%%%%%%%start of main while loop%%%%%%%%%%%%%%%%%%%%%%%%%%%
% keyboard
while noiter<iterbnd && ...
        (mingapp/(abs(objval)+1)> tolerFR) ...
        && (Fmuprnorm /(abs(objval)+1) > tolerFR) && stalling <=  5 ...
        && (~earlystop)
    
    noiter=noiter+1;
    ubs(noiter) = bestub;
    lbs(noiter) = bestlb;
    
    %%%%%%%%%%%%%%%%%%%
    %%%% evaluate Jacobian
    
    % formulate dual feasibility equation
    for ii = 1:isodim-m+1  % check the column size for the real case
        Btemp = HSMat(NW_noperm(:,ii),varflag) * S ;
        JacD(:,ii) = CRvec( Btemp,varflag );
    end
    
    % formulate primal feasibility equation
    for ii = 1:m-1
        Btemp = W * HSMat(NS_noperm(:,ii),varflag) ;
        JacP(:,ii) = CRvec( Btemp,varflag );
    end
    Jac = [JacD,JacP];
    
    JsJac = Jac'*Jac;
    if verbose == 2
        condJac = condest(JsJac);
    end
    d = sqrt(diag(JsJac));   % for diagonal preconditioning
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% predictor step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Jacrhs_pred = -CRvec( W*S + W*Fmu{4} + Fmu{1}*S, varflag ); % RHS
        
    deltas = ((Jac/diag(d))\Jacrhs_pred)./d;  % diag. preconditioned soln
    
    % assign directions deltas=[dv;du]
    dv = deltas(1:isodim-m+1);
    du = deltas(end-m+2:end);
    relresnormpred = norm(JsJac*deltas-Jac'*Jacrhs_pred);
    
    % form the directions (original variables); backsolve
    dW = HSMat(NW_noperm*dv,varflag) + Fmu{1};
    dS = HSMat(NS_noperm*du,varflag) + Fmu{4}; % exposing vector direction
    % dlambda, dy are not used
    %dlambda = Nlambda_noperm*dv + Fmu{2};
    %dy = Ny_noperm*du + Fmu{3};
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Corrector step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Waff = dW;
    Saff = dS;
    
    %obtain step length that keeps the W,S positive definite
    stepDFaff = 1;
    Wp = W + stepDFaff*dW;
    [~,choltrueW] = chol(Wp);
    while choltrueW  > 0 && stepDFaff > 1e-13
        stepDFaff = .93*stepDFaff;  % original backtracking
        Wp = W + stepDFaff*dW;
        [~,choltrueW] = chol(Wp);
    end
    stepPFaff = 1;
    Sp = S + stepPFaff*dS;
    [~,choltrueS] = chol(Sp);
    while  choltrueS > 0 && stepPFaff > 1e-13
        stepPFaff = .93*stepPFaff;  % original backtracking
        Sp = S + stepPFaff*dS;
        [~,choltrueS] = chol(Sp);
    end
    
    stepDFaff = min(1,stepDFaff);
    stepPFaff = min(1,stepPFaff);
    
    
    muu_aff = real(trace((W+stepDFaff*Waff)*(S+stepPFaff*Saff)))/n;
    sigmaa = (muu_aff/muu)^3;
    
    %Jacrhs = -CRvec( Fmu{5} + W*Fmu{4} + Fmu{1}*S, varflag ); % RHS
    CStemp = (W*S + Waff*Saff) - sigmaa*muu*In;
    Jacrhs_corr = -CRvec( CStemp + W*Fmu{4} + Fmu{1}*S, varflag ); % RHS
    
    %d = sqrt(diag(JsJac)); % for diagonal preconditioning (computed above)
    deltas = ((Jac/diag(d))\Jacrhs_corr)./d;  % diag. preconditioned soln
    
    % assign directions deltas=[dv;du]
    dv = deltas(1:isodim-m+1);
    du = deltas(end-m+2:end);
    relresnormcorr = norm(JsJac*deltas-Jac'*Jacrhs_corr);
    
    % form the directions (original variables); backsolve
    dW = HSMat(NW_noperm*dv,varflag) + Fmu{1};
    dS = HSMat(NS_noperm*du,varflag) + Fmu{4}; % exposing vector direction
    % backsolve
    dlambda = Nlambda_noperm*dv + Fmu{2};
    dy = Ny_noperm*du + Fmu{3};
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %obtain step length that keeps the W,S positive definite
    Wp = W + stepDF*dW;
    [~,choltrueW] = chol(Wp);
    while choltrueW  > 0 && stepDF > 1e-13
        stepDF = .8*stepDF;  % original backtracking
        Wp = W + stepDF*dW;
        [~,choltrueW] = chol(Wp);
    end
    Sp = S + stepPF*dS;
    [~,choltrueS] = chol(Sp);
    while  choltrueS > 0 && stepPF > 1e-13
        stepPF = .8*stepPF;  % original backtracking
        Sp = S + stepPF*dS;
        [~,choltrueS] = chol(Sp);
    end
    
    
    %% Step legnth under exactly dual feas. (flag setup)
    if stepDF > (1/.97) && flagsteponeDF == 0  % Newton step = 1 not taken yet
        if verbose == 2
            fprintf('FIRST Newton step = 1 taken for dual feas.\n');
        end
        flagsteponeDF = 1;  % indicates step = 1 taken
        stepDF = 1;  % take a Newton step first time
    else
        if flagsteponeDF == 1 % Newton step taken last iteration
            if verbose == 2
                fprintf('exact df starts \n');
            end
            flagsteponeDF = 2;
        end
        stepDF = .97*stepDF;   % be safe for now????????????
    end
    %stepDF = step;
    
    %% Step legnth under exactly primal feas. (flag setup)
    %step = steptemp;
    if stepPF > (1/.97) && flagsteponePF == 0  % Newton step = 1 not taken yet
        if verbose == 2
            fprintf('FIRST Newton step = 1 taken for priaml feas.\n');
        end
        flagsteponePF = 1;  % indicates step = 1 taken
        stepPF = 1;  % take a Newton step first time
    else
        if flagsteponePF == 1 % Newton step taken last iteration
            if verbose == 2;  fprintf('exact pf starts \n');   end
            flagsteponePF = 2;
        end
        stepPF = .97*stepPF;   % be safe for now????????????
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% iterate/residual update
    v = v + stepDF*dv;
    u = u + stepPF*du;
    % (y,lambda,W,S) are updated below depending on the stepsizes
    %y = y + step*dy;
    %lambda = lambda + step*dlambda;  % update v for simple/exact rho update
    %W = W + step*dW;
    %S = S + stepPF*dS; %S = HMat(MGamma'*y);
    
    if (mingapp/(abs(objval)+1)> 100*tolerFR)
        v(abs(v)<1e-14)=0;
        u(abs(u)<1e-14)=0;
    end
    
    
    %%% dual iterate/residual update
    if flagsteponeDF >= 2   % step = 1 taken at least one iter previously
        W = What + HSMat(NW_noperm*v,varflag); % exact primal feas.
        lambda = lambdahat + Nlambda_noperm*v;
        Fmu{1} = zeros(n);
        Fmu{2} = 0;
        bestlb = max(bestlb,rhspos*lambda);
    else
        W = W + stepDF*dW;
        lambda = lambda + stepDF*dlambda;
        Fmu{1} = HSMat(NW_noperm*v,varflag) + What - W; %MGamma*Hvec(W)-gamma;
        Fmu{2} = Nlambda_noperm*v + lambdahat - lambda ;
    end
    
    %%% primal iterate/residual update
    if flagsteponePF >= 2 % step = 1 taken at least one iter previously
        y = yhat + Ny_noperm*u;
        S = Shat + HSMat(NS_noperm*u,varflag);
        Fmu{3} = zeros(m,1);
        Fmu{4} = zeros(n);
        bestub = min(bestub,gamma'*y);
    else
        y = y + stepPF*dy;
        S = S + stepPF*dS;
        Fmu{3} = Ny_noperm*u + yhat - y;%trace(HMat(MGamma'*y))-1;
        Fmu{4} = HSMat(NS_noperm*u,varflag) + Shat - S;%trace(HMat(MGamma'*y))-1;
    end
    
    Fmu{5} = W*S-muu*In;
    
    gapp = real(trace(W*S));
    mingapp = min(gapp,bestub-bestlb);
    muu = gapp/n;
    
    
    
    %calculate the current objective
    objval = gamma'*y;%rhspos*lambda;
    
    % vectorized residuals
    Fmuprvec = [HSvec(Fmu{1},varflag);Fmu{2};Fmu{3};...
        HSvec(Fmu{4},varflag); CRvec(Fmu{5},varflag)];
    Fmuprnorm = norm(Fmuprvec);
    
    
    
    
    % update stalling condition
    mineigs = min([eig(W);eig(S)]);
    if mineigs < eps  % use eps machine epsilon ?
        stalling = stalling + 1;
        if verbose == 2
            fprintf('min eig of W or S is %g. Stalling count=%d \n', ....
                mineigs,stalling)
        end
    else
        stalling = 0; % restart to 0
    end
    
    if verbose >= 2
        fprintf('%-9i',noiter)
        fprintf('%-15e',-log10(mingapp/(abs(objval)+1)),...
            stepDF,stepPF,sigmaa,...
            objval,norm(Fmu{3})+norm(Fmu{4},'fro'),...
            norm(Fmu{1},'fro')+Fmu{2},...
            condJac,relresnormpred,relresnormcorr);
        fprintf('\n')
    end
    
    % check early stopping condition using the objective values
    if bestlb >1e-8 && flagsteponeDF == 2  % if dual feasibility achieved and lower bound is positive enough
        earlystop = true;  % stop early.
        if verbose == 2
            fprintf('In FRGNPC.m: STOP early - dual feasibility achieved and dual objective %e is positive.\n',...
                bestlb)
        end
    end
    if bestub <-1e-8 && flagsteponePF == 2  % if dual feasibility achieved and lower bound is positive enough
        earlystop = true;  % stop early.
        if verbose == 2
            fprintf('In FRGNPC: STOP early - primal feasibility achieved and primal objective %e is negative.\n',...
                bestub)
        end
    end
    
    objdiff = gamma'*y - rhspos*lambda ; % primal-dual objective difference
    if abs(objdiff) <0.5*1e-8 && rhspos*lambda >1e-8
        earlystop = true;  % stop early.
        if verbose == 2
            fprintf('In FRGNPC: STOP early - gap is small and dual objective value %e is positive.\n',...
                rhspos*lambda)
        end
        bestlb = rhspos*lambda; % not a true obj value! this is to give flag = 0 below
    end
    
    
    
    % update muu and sigmaa depending on stepsize
    minstep = min(stepDF,stepPF);
    if minstep < .001
        stepDF = .9; stepPF = .9;
    elseif minstep < .8
        stepDF = 1.0; stepPF = 1.0;
    else
        stepDF = 1.2; stepPF = 1.2;
    end
    
    
    gapps(noiter) = gapp;
    mingapps(noiter) = mingapp;
    stepDFs(noiter) = stepDF;
    stepPFs(noiter) = stepPF;
    sigmaas(noiter) = sigmaa;
    objvals(noiter) = objval;
    condJacs(noiter) = condJac;
    relresnormpreds(noiter) = relresnormpred;
    relresnormcorrs(noiter) = relresnormcorr;
end          % end of main while loop
% keyboard


% flag meaning
% -2 : inconclusive
% -1 : infeasible
% 0 : Slater holds
% 1 : reducible (obtain a nontrivial V)

if bestub < -1e-9
    fprintf('The problem is infeasible \n')
    flag = -1;
    Zexp = HSMat(MGamma'*y,varflag);  % this is a certificate of infeasibility
    VGg = NaN;
elseif bestlb > 1e-8
    flag = 0;  % Slater considion holds
    Zexp = zeros(n,n);
    VGg = eye(n);
elseif min(eig(W+lambda*P)) > Slatertoler && flagsteponeDF == 2
    flag = 0;  % Slater considion holds
    Zexp = zeros(n,n);
    VGg = eye(n);
elseif abs( objval ) < 1e-11    % pick a proper tolerance
    flag = 1;
    Zexp = HSMat(MGamma'*y,varflag);      % exposing vector
    Zexp = Zexp - min(eig(Zexp))*In; % this is the due to the trace constraint;
    
    [UZ,DZ] = eig(Zexp);
    %rankZ = length(find(DZ > 1e-13 ));  % question
    rankZ = rank(DZ,n*eps(norm(DZ)));  % numerical rank
    rankZ = max(1,rankZ); % prevent from not choosing any column from exposing vector
    
    feaspoint = W+lambda*In;  % comes from the dual constraint A(W)+lambda*A(In)=b
    feaspoint  = (feaspoint+feaspoint')/2;  % symmetrize
    % (Safeguard) The if statement below is to avoid the numerical rank error arising
    % in rankZ. It determines the rankZ by determiing the rank of the feasible point feaspoint 
    if norm(MGamma*HSvec(feaspoint,varflag) - gamma) / norm(gamma) < tolerFR    % if feaspoint is feasible 
        [~,Dfeaspoint] = eig(feaspoint,'vector');
        rankfeaspoint = length( find(Dfeaspoint>tolerFR) );  % determine the rank of feaspoint
        rankZ = min(rankZ,n-rankfeaspoint);  % determine the rank of the exposing vector
    end
    
    VGg = UZ(:,1:n-rankZ);     % facial vector
    
else
    flag = -2;  % inconclusive
    Zexp = zeros(n,n);
    VGg = eye(n);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(lbs) < noiter
    lbs(length(lbs)+1:noiter) = lbs(end);
end
if length(ubs) < noiter
    ubs(length(ubs)+1:noiter) = ubs(end);
end

if verbose >= 1
    fprintf('%-9s','iter#')
    fprintf('%-15s',VarNames{:})
    fprintf('\n')
    fprintf('%-9i',noiter)
    fprintf('%-15e',-log10(mingapp/(abs(objval)+1)),stepPF,stepDF,sigmaa,...
        objval,...
        norm(Fmu{3})+norm(Fmu{4},'fro'),...
        norm(HSMat(NW_noperm*v,varflag) + What - W,'fro')+Nlambda_noperm*v + lambdahat - lambda,...
        condJac,...
        relresnormpred,relresnormcorr);
    fprintf('\n')
end


if verbose >= 1
    % print the reason for stopping
    if noiter < iterbnd && stalling <= 5
        if earlystop == false
            fprintf('in FRGNPC.m: Toler achieved - algorithm ends\n')
        elseif earlystop == true
            fprintf('in FRGNPC.m: Early stopping \n')
        end
    end
    if noiter >= iterbnd  && stalling <= 5
        fprintf('in FRGNPC.m: ERROR;  exceeded iter bnd: %i  \n',...
            noiter);
    end
    if stalling > 5
        fprintf('in FRGNPC.m: stopping due to stalling \n ')
    end
    
    % print out the conclusion
    if flag == -1
        fprintf('Conclusion: Problem is infeasible! \n')
    elseif flag == 0
        fprintf('Conclusion: Slater condition holds \n')
    elseif flag == 1
        fprintf('Conclusion: Slater condition fails \n')
        fprintf('            facial range vector V dimension: %d-by-%d \n',size(VGg))
        fprintf('            variable reduction: from %d to %d \n',n,size(VGg,2))
    elseif flag == -2
        fprintf('Conclusion: Inconclusive \n')
    end
    if verbose >= 2
        fprintf('mineig, cond#: W,Z, %g %g; %g %g \n',...
            min(eig(W)),condest(W),...
            min(eig(S)),condest(S))
    end
end

% important output
Out.flag = flag;
Out.Zexp = Zexp;
Out.W = W;
Out.y = y;
Out.lambda = lambda;
Out.S = S;
Out.V = VGg;

% auxiliary output
Out.mingapps = mingapps;
Out.muu = muu;
Out.objval = objval;
Out.ubs = ubs;
Out.lbs = lbs;
Out.noiters = noiter;
Out.gapps = gapps;
Out.stepPFs = stepPFs;
Out.stepDFs = stepDFs;
Out.sigmaas = sigmaas;
Out.objvals = objvals;
Out.condJacs = condJacs;


end   % end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
