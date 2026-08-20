clear

% the genetrating geometry is from the demo_mult_connect.m

iseed = 8675309;
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
L = 9;
theta = 2*pi*rand();
ctrs = L*rand()*[cos(theta);sin(theta)];
n_pts = [];
ninc = 100;
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


opts = []; 
opts.sing = 'log';

fkern = @(s,t) chnk.lap2d.kern(s,t,'sprime');
start = tic;
sysmat = chunkermat(chnkrs, fkern, opts);
t1 = toc(start);
fprintf('%5.2e s : time to discretize the integral operator\n',t1)




rem = {};
for i = 1: length(chnkr_int)
      rem{i} = onesmat(chnkr_int(i));
end


stb = blkdiag(rem{:});

lhs = 0.5.*eye(chnkrs.npt) + sysmat + stb ;              % change the stabilizer construction due to different # of discretization ptss

rhs = -chnkrs.n(2,:).';


start = tic; 
rho = gmres(lhs, rhs, [], 1e-12, 200);
t1 = toc(start);
fprintf('%5.2e s : time for gmres\n',t1)    


L = max(abs(chnkrs.r),[],"all");

x1 = linspace(-L-0.15, L + 0.15 ,800);
[xx,yy] = meshgrid(x1,x1);
targets = [xx(:).'; yy(:).'];

% identify points in computational domain
in = chunkerinterior(chnkrs,targets);
out = ~in;


kern_grad = @(s,t) chnk.lap2d.kern(s,t, 'sgrad');               % gradient of the single layer




tstart = tic;


uugrad = nan(2,size(xx(:),1));
curr = nan(2,size(xx(:),1));

E_solin = chunkerkerneval(chnkrs,kern_grad,rho,targets(:,in));
E_solout = chunkerkerneval(chnkrs,kern_grad,rho,targets(:,out)) ;

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


h = pcolor(xx,yy,(str_plot_1));
set(h,'EdgeColor','None'); hold on;
plot(chnkrs,'k-','LineWidth',1);
xlabel('x-axis', 'FontSize',12)
ylabel('y-axis', 'FontSize',12)
title('Exterior electric field ', 'FontSize',12)
colorbar

axis equal
pbaspect([1 1 1]) 

set(fig,'Position',[100 100 1300 500])

exportgraphics(gcf,'high contrast starfish large.pdf','ContentType','image',...
    'Resolution',600)
