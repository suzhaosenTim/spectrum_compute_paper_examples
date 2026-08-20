clear


cparams = [];

cparams.eps = 1e-6;
cparams.nover = 0;


rad = 0.1;
delta = 0.25;                        % separation
ctr = 0;


nch = 1:1:8;

rel_err_max = zeros(1, length(nch));
rel_err_min = zeros(1, length(nch));
npts = zeros(1, length(nch));
cond_number = zeros(1, length(nch));

for ii = 1: length(nch)

    start = tic; 
    chnkrorg = chunkerfuncuni(@(t) rad.*circle(t),nch(ii), cparams);
    t1 = toc(start);

    fprintf('%5.2e s : time to build geo\n',t1)
    
    ctr1 = [ctr + delta/2; 0];
    ctr2 = [-ctr - delta/2; 0];

    chnkr = chnkrorg.move([], ctr1 );

    chnkr1 = chnkrorg.move([], ctr2);

    chnkr_int = [chnkr, chnkr1];

    chnkrs = merge(chnkr_int);


    npts(ii)= chnkrs.npt;


    fkern = @(s,t) chnk.lap2d.kern(s,t,'sprime');


    opts = []; 
    opts.sing = 'log';

    start = tic;
    sysmat = chunkermat(chnkrs, fkern, opts);
    t1 = toc(start);
    fprintf('%5.2e s : time to discretize the integral operator\n',t1)


    lhs = 0.5.*eye(chnkrs.npt) + sysmat  ;     % left hand side

    lambda = eig(lhs);

    iid = lambda>1e-5;

    lambda_max = max(lambda);

    lambda_min = min(lambda(iid));





% Compute two families max eigenvalues via formula in Bonnetier paper
% (slightly more general)

    a = sqrt((ctr + (delta/2)-rad)*((ctr + (delta/2)-rad) + 2*rad));

    rhodist = (ctr + (delta/2) -a)/rad ;


    s_para_lower =  0.5 - rhodist^2/2;
    s_para_higher = 0.5 + rhodist^2/2;

    
    max_err_upp = norm(lambda_max - s_para_higher)/norm(s_para_higher);
    max_err_low = norm(lambda_min - s_para_lower)/norm(s_para_lower);

    rel_err_max(ii) = max_err_upp;
    rel_err_min(ii) = max_err_low;

end

fig = figure;
clf

t = tiledlayout(fig,1,2, ...
    'TileSpacing','compact', ...
    'Padding','loose');


nexttile
loglog(npts(1:8), rel_err_max(1:8), '.-', 'LineWidth',1,'MarkerSize',12)
hold on 
loglog(npts(1:8), rel_err_min(1:8), '.-', 'LineWidth',1,'MarkerSize',12)
hold on 
loglog(npts(1:8), (10^22).*npts(1:8).^(-14), '--k', 'LineWidth',1)
 
legend({'Max eigenvalue relative error', 'Min eigenvalue relative error', 'n^{-14}'})
xlabel('N', 'FontSize',12)
ylabel('Relative error','FontSize',12)
title('eigenvalue computation', 'FontSize',12)
axis square


nexttile
plot(chnkrs,'-k', 'LineWidth',1)
hold on 
plot(chnkrs, 'xk', 'LineWidth',1)
title('geometry', 'FontSize',12)
axis equal
pbaspect([1 1 1]) 

set(gcf,'Position', [900  500  900   500])


exportgraphics(gcf,'convergence_figure.pdf','ContentType','vector')