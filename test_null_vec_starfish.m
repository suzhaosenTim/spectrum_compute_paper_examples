clear

% the genetrating geometry is from the demo_mult_connect.m


iseed = 185309;
rng(iseed,'twister');


cparams = [];
cparams.eps = 1.0e-6;
cparams.nover = 0;
                              

narms =5;
amp = 0.25;
scale = .3;
rad = scale*(amp+1);
ntry = 1000;

% make interior boundaries with random locations
chnkr_int = [];
L = 7;
theta = 2*pi*rand();
ctrs = L*rand()*[cos(theta);sin(theta)];
n_pts = [];
ninc = 80;
for i  = 1:ninc
    % interior boundary is rotated starfish with a random number of arms
    phi = 2*pi*rand();
    narms = randi([3,6]);
    chnkr_i = chunkerfunc(@(t) starfish(t,narms,amp,ctrs(:,i),phi,scale), ...
        cparams); 
    

    % track number of points
    n_pts(i) = chnkr_i.npt;
    chnkr_int = [chnkr_int,chnkr_i];

    % try to find location for next chunker
    for j = 1:ntry
        theta = 2*pi*rand();
        tmp = L*rand()*[cos(theta);sin(theta)];
        rmin = min(vecnorm(tmp - ctrs));
        if (rmin > rad * 3); break; end
    end
    if j == ntry; error('Could not place next boundary'); end
    ctrs = [ctrs,tmp];
end
ctrs = ctrs(:,1:end-1);

chnkrs = merge(chnkr_int);

fprintf('Geometry generated\n')
figure(3); clf;
plot(chnkrs)
quiver(chnkrs)
axis equal




fkern = @(s,t) chnk.lap2d.kern(s,t,'sprime');

opts = []; 
opts.sing = 'log';

start = tic;
sysmat = chunkermat(chnkrs, fkern, opts);
t1 = toc(start);
fprintf('%5.2e s : time to discretize the integral operator\n',t1)


lhs = 0.5.*eye(chnkrs.npt) + sysmat  ;     % left hand side

[rho, lambda_max] = eigs(lhs, 1, 'largestreal');

r2 = chnkrs.r(1,:).';
msure_wts = chunkerintegral(chnkrs, rho.*r2);


L = max(abs(chnkrs.r),[],"all");

x1 = linspace(-L-0.15, L + 0.15 ,500);

[xx,yy] = meshgrid(x1,x1);
targets = [xx(:).'; yy(:).'];                       % generate some targets to evaluate the eigenfunctions

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

exportgraphics(gcf,'null_vector_starfish.pdf','ContentType','vector')