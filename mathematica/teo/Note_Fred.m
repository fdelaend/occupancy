niter = 1000;
deltad = zeros(1,niter);
deltaoffd = zeros(1,niter);

N = 1000;
stdval = 0.001;
sigma = stdval^2;
for ii = 1:niter
    A = stdval*randn(N,N);
    A = A-diag(diag(A))+eye(N);
    X = inv(A);
    Xd = diag(X);
    Xoffd = X- diag(Xd);
    
    deltad(ii) = mean(Xd)-(1+sigma*(N-2))/((1-sigma)*(1+sigma*(N-1)));
    deltaoffd(ii) = sum(Xoffd(:))/(N*(N-1))-sigma/((1-sigma)*(1+sigma*(N-1)));
end

figure
histogram(deltad)
xlabel('$\varepsilon_{diag}$','Interpreter','latex')
set(gca,'FontSize',24)
set(gca,'FontName','Times')
set(gca,'FontAngle','italic')

figure
histogram(deltaoffd)
xlabel('$\varepsilon_{off}$','Interpreter','latex')
set(gca,'FontSize',24)
set(gca,'FontName','Times')
set(gca,'FontAngle','italic')
