
%% Test SOS matrices -- SDD matrices

n = 4;                    % number of variables in each polynomial
degree = 2;               % degree of each polynomial


    Nt = 5:5:50;                   % Dimension of matrices
    
    Bound = zeros(length(Nt),4);
    TimeSolver = zeros(length(Nt),4);
    
for index = 1:length(Nt)
    
    N = Nt(index);
    
    %% generate data
    x      = sdpvar(n,1);          % polynomial variables
    Lower  = sdpvar(1,1);          % decision variables
    mBasis = monolist(x,degree/2,1);
    N1 = length(mBasis);

    % for SPOTLESS
    y      = msspoly('y',n);
    z      = msspoly('z',N);
    basis  = monomials(y,1:degree/2);
    Ps     = 0;

    P = sdpvar(N,N);
    for i = 1:N
        for j = i: N
            if i == j
                 tmp = 0.2 + rand(N1,1);
                 %tmp = randi([1,10],N1,1);
    %             tmp = tmp*tmp';     % PSD matrix
                P(i,i) = (tmp.*mBasis)'*(tmp.*mBasis);
                Ps = Ps + (tmp.*basis)'*(tmp.*basis)*z(i)^2;
            else
                tmp = rand(N1,1);
                 P(i,j) = sum(tmp.*mBasis);
                 %P(i,j) = 1;
                 P(j,i) = P(i,j);

                 Ps = Ps + 2*sum(tmp.*basis)*z(i)*z(j);
                 %Ps = Ps + 2*z(i)*z(j);
            end
        end
    end



    %% Method 1: Matrix-original
    opts = sdpsettings('solver','mosek');
    %opts = sdpsettings('solver','sedumi');
    F = sos(P+Lower*eye(N));
    [sol,v,Q] = solvesos(F,Lower,opts);

    %% Method 2: sdd on SOS matrices, converted to scalar 
    z1      = sdpvar(N,1);
    Lower1  = sdpvar(1,1);
    P1      = z1'*P*z1+Lower1*z1'*z1;

    % constructed supported basis 
    combos = combntns(1:N,2);
    I      = num2cell(combos,2);
    monoms = cellfun(@(C)monolist([x;z1(C)],(degree+2)/2), I, 'UniformOutput', 0);
    [sol1,v1,Q1]   = solvesos(sos(P1),[Lower1],opts,[],{monoms});


%% Method 3: solution via Spotless - DSOS/SDSOS
    
    
    
    prog   = spotsosprog;
    prog   = prog.withIndeterminate(y);
    prog   = prog.withIndeterminate(z);  
    

    % New free variable gamma
    [prog,gamma] = prog.newFree(1);
    % DSOS constraint
    prog = prog.withDSOS(Ps + gamma*(z'*z)); % Only line that changes between DSOS,SDSOS,SOS programs

    % MOSEK options
    options = spot_sdp_default_options();
    options.solveroptions.verbose = 1;
    options.verbose = 1;
    options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean up
    options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)

    % Solve program
    sol2 = prog.minimize(gamma, @spot_mosek, options);
    if strcmp(sol2.status,'STATUS_PRIMAL_INFEASIBLE')
        opt_dsos = NaN;
    else
        opt_dsos = double(sol2.eval(gamma));
    end
    %% SDSOS
    prog = spotsosprog;
    prog = prog.withIndeterminate(y);
    prog = prog.withIndeterminate(z);
    [prog,gamma1] = prog.newFree(1);
    % DSOS constraint
    prog = prog.withSDSOS(Ps + gamma1*(z'*z)); % Only line that changes between DSOS,SDSOS,SOS programs
    % Solve program
    sol3 = prog.minimize(gamma1, @spot_mosek, options);
    opt_sdsos = double(sol3.eval(gamma1));
    

    %% Summary
    Bound(index,:) = [value(Lower),value(Lower1),opt_dsos,opt_sdsos];

    TimeSolver(index,:) = [sol.solvertime,sol1.solvertime,sol2.info.wtime,sol3.info.wtime];
    
    save Result0801_v1 Bound TimeSolver Nt n degree
    
end
    