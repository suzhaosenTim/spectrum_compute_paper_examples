clear



cparams = [];

cparams.eps = 1e-6;
cparams.nover = 0;
cparams.maxchunklen = 0.5;



rad = 0.1;           
nx  = 4;             % number of circles in x-direction
ny  = 4;             % number of circles in y-direction (ninc = nx * ny)
                     
spacing = 0.2556;   


Lx = (nx-1)*spacing;
Ly = (ny-1)*spacing ;



chnkrorg = chunkerfunc(@(t) rad.*circle(t), cparams);


chnkr_int = [];

n_pts = [];

for x = 1:nx
    chnkr_i = chunkerfunc(@(t) rad.*circle(t), cparams);
    for y = 1:ny

        cx = -Lx/2 + (x-1)*spacing;
        cy = -Ly/2 + (y-1)*spacing;
        
        chnkr_i = chnkrorg.move([], [cx; cy]);

        chnkr_int = [chnkr_int, chnkr_i];
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


lhs = 0.5.*eye(chnkrs.npt) + sysmat  ;     % left hand side

[rho, lambda_max] = eigs(lhs, 1, 'largestreal');


L = max(abs(chnkrs.r),[],"all");

x1 = linspace(-L-0.15, L + 0.15 ,500);
[xx,yy] = meshgrid(x1,x1);
targets = [xx(:).'; yy(:).'];                      % generate some targets to evaluate the eigenfunctions

in = chunkerinterior(chnkrs, targets);              % identify interior and exterior (there is a hacky way to do this via FMM)
out = ~in;

eval_kern = @(s,t) chnk.lap2d.kern(s,t,'sgrad');  % gradient of the single layer

start2 = tic;
gamma_eig_in= chunkerkerneval(chnkrs, eval_kern, rho, targets(:,in));
gamma_eig_out = chunkerkerneval(chnkrs, eval_kern, rho, targets(:,out));
t3 = toc(start2);
fprintf('%5.2e s : time to plot \n',t3)

gamma_eigvec = nan(2,size(xx(:),1));

gamma_eigvec(:,in(:)) = reshape(gamma_eig_in, 2, []);

gamma_eigvec(:,out(:)) = reshape(gamma_eig_out, 2,[]);

strE = abs(gamma_eigvec(:, in(:)));

strE_out = abs(gamma_eigvec(:, out(:)));

str_plot = nan(size(xx));

str_plot(in) = sqrt(strE(1,:).^2 + strE(2,:).^2); 
str_plot(out) = sqrt(strE_out(1,:).^2 + strE_out(2,:).^2); 

fig = figure;
clf

t = tiledlayout(fig,1,2, ...
    'TileSpacing','compact', ...
    'Padding','loose');




nexttile
h = pcolor(xx, yy, (str_plot));
set(h,'EdgeColor','none');
hold on
plot(chnkrs,'k')
xlabel('x-axis', 'FontSize',12)
ylabel('y-axis', 'FontSize',12)
title('Gamma operator eigenfunction ', 'FontSize',12)
axis equal
colorbar
pbaspect([1 1 1]) 




nexttile
h = pcolor(xx, yy, log10(str_plot));
set(h,'EdgeColor','none');
hold on
plot(chnkrs,'k')
xlabel('x-axis', 'FontSize',12)
ylabel('y-axis', 'FontSize',12)
title('Gamma operator eigenfunction in log scale ', 'FontSize',12)
axis equal
colorbar
pbaspect([1 1 1]) 




set(gcf,'Position', [900  500  900   500])

exportgraphics(gcf,'null_vector_ordered.pdf','ContentType','vector',...
    'Resolution',600)