function [H,err,pD] = H_als(p1,p2)
    max_iter = 50;
    tol = 1e-20;
    [Nch,Npx] = size(p1);
    vind = (sum(p1>0 ,1)==Nch & sum(p2>0 ,1)==Nch);
    P = p1;
    Q = p2;
    N = P;
    pD = speye(Npx);
    errs = Inf(max_iter+1,1); 
    H=eye(3);

    % solve the homography using ALS
    n_it = 1; d_err = Inf;
    while ( n_it-1<max_iter && d_err>tol )
        n_it = n_it+1; 

        Dr = (1:Npx*Nch)'; 
        Dc = repmat(1:Npx,[Nch,1]);
        Dc = Dc(:); 
        A = sparse(Dr,Dc,N(:),Npx*Nch,Npx);
        B = Q(:);
        A1 = A'*A;
        D = A1\(A'*B);
        D = spdiags(D,0,Npx,Npx);

        P_d = P*D;
        M = Q(:,vind)/P_d(:,vind); 
        N = M*P;

        NDiff = (N*D-Q).^2; 
        errs(n_it) = mean(mean(NDiff(:,vind))); 
        d_err = errs(n_it-1) - errs(n_it);
    end

    H = M;
    err = errs(n_it);
    pD = D;
end