clear




cparams = [];

cparams.eps = 1e-6;
cparams.nover = 0;


rad = 0.1;
delta = 0.25;                        % separation
ctr = 0;

% two circle case
start = tic; 
chnkrorg = chunkerfunc(@(t) rad.*circle(t),cparams);
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)
    
ctr1 = [ctr + delta/2; 0];
ctr2 = [-ctr - delta/2; 0];

chnkr = chnkrorg.move([], ctr1 );

chnkr1 = chnkrorg.move([], ctr2);

chnkr_int = [chnkr, chnkr1];

chnkrs = merge(chnkr_int);
chnkrs = refine(chnkrs, struct("nover", 1));

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

lambda_corr = lambda(iid);                  % throw out the known 0 eigenvalue 

lambda_sing = lambda_corr(1:2:end);

[n, ~] = size(lambda_sing);

% Compute two families max eigenvalues via formula in Bonnetier paper
% (slightly more general)

a = sqrt((ctr + (delta/2)-rad)*((ctr + (delta/2)-rad) + 2*rad));

rhodist = (ctr + (delta/2) -a)/rad ;


d2 = 1:n;

s_para_lower =  0.5 - rhodist.^(2.*d2)/2;
s_para_higher = 0.5 + rhodist.^(2.*d2)/2;

lambda_sing = sort(lambda_sing);

iid2 = (abs(lambda_sing - 0.5)>1e-5);               % throw out the trivial 1/2 accumulation point
lambda_sing = lambda_sing(iid2);

tru_eig = sort([s_para_lower, s_para_higher]);
iid3 = (abs(tru_eig - 0.5)>1e-5);
tru_eig = tru_eig(iid3);

rel_err1 = abs(lambda_sing - tru_eig.')./(tru_eig.');
N = length(rel_err1);


% ellipse case
a = 2.75;
b = 1.25;

chnkr3 = chunkerfunc(@(t)ellipse(t, a,b), cparams);         % an ellipse



start = tic;
sysmat2 = chunkermat(chnkr3, fkern, opts);
t1 = toc(start);
fprintf('%5.2e s : time to discretize the integral operator\n',t1)


lhs_e = 0.5.*eye(chnkr3.npt) + sysmat2  ;      % left hand side

lambda_e = eig(lhs_e);

iid = lambda_e>1e-5;

lambda_corr_e = lambda_e(iid);                  % throw out the known 0 eigenvalue 
d2 = length(lambda_e)/2;



eig_true_e_upp = 0.5 + 0.5.*((a-b)/(a+b)).^(1:d2);

eig_true_e_low = 0.5 - 0.5.*((a-b)/(a+b)).^(1:d2);

eig_true_e = sort([eig_true_e_low, eig_true_e_upp]);
lambda_corr_e = sort(lambda_corr_e);

iid4 = (abs(lambda_corr_e - 0.5)>1e-5);
iid5 = (abs(eig_true_e- 0.5)>1e-5);

eig_true_e = eig_true_e(iid5); 
lambda_corr_e = lambda_corr_e(iid4);

rel_err2 = abs(lambda_corr_e - eig_true_e.')./(eig_true_e.');
N2 = length(rel_err2);





fig = figure;
clf

t = tiledlayout(fig,1,2, ...
    'TileSpacing','compact', ...
    'Padding','loose');


nexttile
scatter(1:N, rel_err1, "*")

xlabel('n-th non-trivial eigenvalue', 'FontSize',12)
ylabel('Relative error','FontSize',12)
ylim([1e-16, 1])
yscale log
%title('Relative error(two circles)', 'FontSize',12)
axis square


nexttile
scatter(1:N2, rel_err2, "*")
xlabel('n-th non-trivial eigenvalue', 'FontSize',12)
ylabel('Relative error','FontSize',12)
ylim([1e-16, 1])
%title('Relative error(ellipse)', 'FontSize',12)
axis square

yscale log
pbaspect([1 1 1]) 


%set(gcf,'Position', [900  500  900   500])
set(gcf,'Position',[541 592 927 317])
fontname(gcf, 'Helvetica')

exportgraphics(gcf,'eigvals error.pdf','ContentType','vector')
