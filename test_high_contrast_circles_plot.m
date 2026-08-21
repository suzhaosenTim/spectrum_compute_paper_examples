clear

% the style of genetrating geometry is from the demo_mult_connect.m

iseed = 8675309;
rng(iseed,'twister');

cparams = [];
cparams.eps = 1.0e-6;
cparams.nover = 0;

chnkr_int = [];
L = 1.5;                % domain size (square box (2x 2))
ntry = 1000;          
rad = 0.1;          
ninc = 16;              % number of inclusions

ctrs = zeros(2,0);    

for i = 1:ninc
    placed = false;    
    chnkr_i = chunkerfunc(@(t) rad .* circle(t), cparams);
    chnkr_i = refine(chnkr_i, struct('nover',1));

    for j = 1:ntry
        tmp = (L - 2*rad)*rand(2,1);
        if isempty(ctrs)
            dist_ok = true;
        else
             dists = vecnorm(tmp - ctrs);
             dist_ok = all(dists > 2*rad); 
        end

        if dist_ok
           chnkr_i = chnkr_i.move([],tmp );
           chnkr_int = [chnkr_int, chnkr_i];
           ctrs(:,end+1) = tmp;
           placed = true;
           break
        end
    end

    if ~placed
        error('Could not place circle without touching, enlarge the wall L');
    end

end


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
lhs = 0.5.*eye(chnkrs.npt) + sysmat + stb ;





x1 = linspace(-0.2, 1.4 ,500);
[xx,yy] = meshgrid(x1,x1);
targets = [xx(:).'; yy(:).'];


   

rhs2 = -chnkrs.n(2,:).';

start = tic; 
rho2 = lhs\rhs2;
t1 = toc(start);
fprintf('%5.2e s : time for dense solve\n',t1)




% identify points in computational domain
in = chunkerinterior(chnkrs,{x1,x1});
out = ~in;


kern_grad = @(s,t) chnk.lap2d.kern(s,t, 'sgrad');               % gradient of the single layer




tstart = tic;


uugrad = nan(2,size(xx(:),1));
curr = nan(2,size(xx(:),1));

E_solin = chunkerkerneval(chnkrs,kern_grad,rho2,targets(:,in));
E_solout = chunkerkerneval(chnkrs,kern_grad,rho2,targets(:,out)) ;

uugrad(:,in(:)) = (reshape(E_solin,2,[]) + [0;1]);
uugrad(:,out(:)) = (reshape(E_solout,2,[]) + [0;1]);
tstream = toc(tstart);
fprintf('%5.2e s : Evaluate the electric field \n',tstream)

% 
% 
% 
% 
% 
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

h = pcolor(xx,yy,(str_plot_1));
set(h,'EdgeColor','None'); hold on;
plot(chnkrs,'k-','LineWidth',1);
xlabel('x-axis', 'FontSize',12)
ylabel('y-axis', 'FontSize',12)
title('Exterior electric field ', 'FontSize',12)
colorbar

axis equal
pbaspect([1 1 1]) 


nexttile
h = pcolor(xx,yy,log10(str_plot));
set(h,'EdgeColor','None'); hold on;
plot(chnkrs,'k-','LineWidth',1);
xlabel('x-axis', 'FontSize',12)
ylabel('y-axis', 'FontSize',12)
title('log10 Electric field ', 'FontSize',12)
colorbar
axis equal
pbaspect([1 1 1]) 



src = [3;-1];                                                   % analytic solution test
[~,grad] = chnk.lap2d.green(src,chnkrs.r,true);

nx = chnkrs.n(1,:);
ny = chnkrs.n(2,:);

rhs = grad(:, :, 1).*(nx.') + grad(:, :, 2).*(ny.');

start = tic; 
rho = lhs\rhs;
t1 = toc(start);

fprintf('%5.2e s : time for dense solve\n',t1)    

skern = @(s,t) chnk.lap2d.kern(s, t, 's');



[~, nt] = size(targets);
true_sol = zeros(nt,1);
utarg = zeros(nt, 1);


nexttile
 h = pcolor(xx,yy,nan(size(xx)));
 set(h,'EdgeColor','None'); 
 hold on;
 plot(chnkrs,'k-','LineWidth',1);
 xlabel('x-axis', 'FontSize',12)
 ylabel('y-axis', 'FontSize',12)
 title('Field Error (Relative)', 'FontSize',12)
 colorbar
axis equal
pbaspect([1 1 1]) 
    
for j = 1:length(chnkr_int)
    in = chunkerinterior(chnkr_int(j), targets);
    start4 = tic;
    solin_1 = chunkerkerneval(chnkrs,skern,rho,targets(:,in));

    val1 = chnk.lap2d.green(src, ctrs(:,j));

    solin_t_1 = chunkerkerneval(chnkrs,skern,rho,ctrs(:,j));

    sol_subtract = solin_1 - (solin_t_1 - val1);    % subtract off the constant



    true_sol(in) = chnk.lap2d.green(src, targets(:, in));

    utarg(in) = sol_subtract ;
    uerr = (utarg - true_sol)./true_sol;
    uerr = reshape(uerr,size(xx));

    h.CData = log10(abs(uerr));
end


set(fig,'Position',[100 100 1300 500])

exportgraphics(gcf,'high contrast circles.pdf','ContentType','vector')