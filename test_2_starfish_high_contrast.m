clear




cparams = [];
cparams.eps = 1.0e-6;
cparams.nover = 0;
% 

narms =6;
amp = 0.25;

start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp), cparams); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

chnkr1 = chnkr.move([], [3;1]);

chnkr_int = [chnkr, chnkr1];

chnkrs = merge(chnkr_int);



fkern = @(s,t) chnk.lap2d.kern(s,t,'sprime');


opts = []; 
opts.sing = 'log';

start = tic;
sysmat = chunkermat(chnkrs, fkern, opts);
t1 = toc(start);
fprintf('%5.2e s : time to discretize the integral operator\n',t1)


    
rem = [];
for i = 1: length(chnkr_int)
      rem = [rem, onesmat(chnkr_int(i))];
end

ninc = length(chnkr_int);
    nchs = zeros(1,ninc);                   % indices to extract the corressponding corrector
    for i=1:ninc
        nchs(i) = chnkr_int(i).nch;
    end

    ngl = chnkr_int(1).k;

    r2 = cumsum(nchs.*ngl);

    r1 = [1, r2(1:end-1) + 1];

    A = {};
    for k = 1:length(r1)
        A{k} = rem(:, r1(k):r2(k));
    end

    stb = blkdiag(A{:});

    lhs = 0.5.*eye(chnkrs.npt) + sysmat + stb ;     % left hand side

    src = [3;-1];
    [~,grad] = chnk.lap2d.green(src,chnkrs.r,true);

    nx = chnkrs.n(1,:);
    ny = chnkrs.n(2,:);

    rhs = grad(:, :, 1).*(nx.') + grad(:, :, 2).*(ny.');

    start = tic; 
    rho = lhs\rhs;
    t1 = toc(start);
    fprintf('%5.2e s : time for dense solve\n',t1)    


x1 = linspace(-1.5, 4.5 ,300);
x2 = linspace(-1.5, 4.5, 300);
[xx,yy] = meshgrid(x1,x2);
targets = [xx(:).'; yy(:).'];


% identify points in computational domain
in1 = chunkerinterior(chnkr, targets);
in2 = chunkerinterior(chnkr1, targets);

skern = @(s,t) chnk.lap2d.kern(s, t, 's');

start4 = tic;
solin_1 = chunkerkerneval(chnkrs,skern,rho,targets(:,in1));
solin_2 = chunkerkerneval(chnkrs,skern,rho,targets(:,in2));
t4 = toc(start4);
fprintf('%5.2e s : time to evaluate\n',t4)



val1 = chnk.lap2d.green(src, [0;0]);
val2 = chnk.lap2d.green(src, [3;1]);

solin_t_1 = chunkerkerneval(chnkrs,skern,rho,[0;0]);

solin_t_2 = chunkerkerneval(chnkrs,skern,rho,[3;1]);

sol_N =solin_1 - (solin_t_1 - val1);
sol_N2 = solin_2 - (solin_t_2 - val2);                  % manual subtraction 

[~, nt] = size(targets);
true_sol = zeros(nt,1);
utarg = zeros(nt, 1);
[val1, ~] = chnk.lap2d.green(src, targets(:, in1));
[val2, ~] = chnk.lap2d.green(src, targets(:, in2));
true_sol(in1) = val1;
utarg(in1) = sol_N ;


true_sol(in2) = val2;
utarg(in2) = sol_N2 ;
uerr = (utarg - true_sol)./true_sol;
uerr = reshape(uerr,size(xx));



rhs2 = -chnkrs.n(2,:).';                                % test for the high contrast computation

start = tic; 
rho2 = lhs\rhs2;
t1 = toc(start);
fprintf('%5.2e s : time for dense solve\n',t1)


in = chunkerinterior(chnkrs, {x1, x2});
out = ~in;



kern_grad = @(s,t) chnk.lap2d.kern(s,t, 'sgrad');               % gradient of the single layer




tstart = tic;


uugrad = nan(2,size(xx(:),1));
curr = nan(2,size(xx(:),1));



E_solin = chunkerkerneval(chnkrs,kern_grad,rho2,targets(:,in));
E_solout = chunkerkerneval(chnkrs,kern_grad,rho2,targets(:,out)) ;

uugrad(:,in(:)) = (reshape(E_solin,2,[]) + [0;1]);
uugrad(:,out(:)) = (reshape(E_solout,2,[]) + [0;1]);
telec = toc(tstart);
fprintf('%5.2e s : Evaluate the electric field \n',telec)

strE = abs(uugrad(:, in(:)));

strE_out = abs(uugrad(:, out(:)));

str_plot = nan(size(xx));

str_plot(in) = sqrt(strE(1,:).^2 + strE(2,:).^2); 
str_plot(out) = sqrt(strE_out(1,:).^2 + strE_out(2,:).^2); 

str_plot_1 = nan(size(xx));

str_plot_1(out) = sqrt(strE_out(1,:).^2 + strE_out(2,:).^2); 

fig = figure;
clf

t = tiledlayout(fig,1,3, ...
    'TileSpacing','compact', ...
    'Padding','loose');

nexttile

h = pcolor(xx,yy,str_plot_1);
set(h,'EdgeColor','None'); hold on;
plot(chnkrs,'k','LineWidth',1);
xlabel('x-axis', 'FontSize',12)
ylabel('y-axis', 'FontSize',12)
title('Exterior electric field ', 'FontSize',12)
colorbar
axis equal
pbaspect([1 1 1]) 


nexttile
h = pcolor(xx,yy,log10(str_plot));
set(h,'EdgeColor','None'); hold on;
plot(chnkrs,'k','LineWidth',1);
xlabel('x-axis', 'FontSize',12)
ylabel('y-axis', 'FontSize',12)
title('log10 Electric field ', 'FontSize',12)
colorbar
axis equal
pbaspect([1 1 1]) 


nexttile
h = pcolor(xx,yy,log10(abs(uerr)));
set(h,'EdgeColor','None'); hold on;
plot(chnkrs,'k','LineWidth',1);
xlabel('x-axis', 'FontSize',12)
ylabel('y-axis', 'FontSize',12)
title('Field Error (Relative)', 'FontSize',12)
colorbar
axis equal
pbaspect([1 1 1]) 


set(fig,'Position',[100 100 1300 500])

exportgraphics(gcf,'high contrast starfish.pdf','ContentType','vector')


