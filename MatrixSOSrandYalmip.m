
%% Random test

clc;clear;

Num = 30:10:50;                       % matrix dimension
ResultLower = zeros(length(Num),3);   % Matrix - orignal; Matrix decomposition; Scalar-csp
ResultTime  = zeros(length(Num),3);

for Index = 1:length(Num)
    
    N = Num(Index);
    n = 2;                    % number of variables in each polynomial
    degree = 2;               % degree of each polynomial

    %% generate data
    x      = sdpvar(n,1);          % polynomial variables
    Lower  = sdpvar(1,1);          % decision variables
    mBasis = monolist(x,degree/2);

    N1  = length(mBasis);
    tmp = rand(N1);
    tmp = tmp*tmp';     % PSD matrix
    P   = mBasis'*tmp*mBasis;
    for i = 2:N
        tmp = rand(N1);
        tmp = tmp*tmp';     % PSD matrix
        P = blkdiag(P,mBasis'*tmp*mBasis);
    end

    for i = 2:N
        P(1,i) = rand(1,N1)*mBasis;
        P(i,1) = P(1,i);
    end

    %% Method 1: Matrix-original
    opts = sdpsettings('solver','mosek');
    %opts = sdpsettings('solver','sedumi');
    F = sos(P+Lower*eye(N));
    [sol,v,Q] = solvesos(F,Lower,opts);

    %% method 2: based on decomposition
    Lower1 = sdpvar(1,1);
    Pi = cell(N-1,1);
    for i = 1:N-1
        Pi{i}      = sdpvar(2);
        Pi{i}(1,2) = P(1,i+1);Pi{i}(2,1) = Pi{i}(1,2); Pi{i}(2,2) = P(i+1,i+1) + Lower1;
        Pi{i}(1,1) = polynomial(x,degree);
    end

    % constraint
    F   = [];
    tmp = 0;
    for i = 1:N-1
        tmp = tmp + Pi{i}(1,1);
        F   = [F, sos(Pi{i})];
    end
    F = [F, coefficients(tmp-P(1,1)-Lower1,mBasis) == 0];
    [sol1,v,Q] = solvesos(F,Lower1,opts);

    %% method 3: converted into scalar polynomials
    z      = sdpvar(N,1);
    Lower3 = sdpvar(1,1);
    P3     = z'*P*z+Lower3*z'*z;
    F      = sos(P3);
    opts.sos.csp = 1;                    % correlative sparsity pattern
    [sol3,v3,Q3] = solvesos(F,Lower3,opts);

    ResultLower(Index,:) = [value(Lower),value(Lower1),value(Lower3)];
    ResultTime(Index,:)= [sol.solvertime,sol1.solvertime,sol3.solvertime];

end

