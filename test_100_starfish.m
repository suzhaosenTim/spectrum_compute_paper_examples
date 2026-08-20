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
L = 10;
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

lhs = 0.5.*eye(chnkrs.npt) +sysmat + stb ;


src = [-4;-6];
[~,grad] = chnk.lap2d.green(src,chnkrs.r,true);

nx = chnkrs.n(1,:);
ny = chnkrs.n(2,:);

rhs = grad(:, :, 1).*(nx.') + grad(:, :, 2).*(ny.');


start = tic; 
rho = gmres(lhs, rhs, [], 1e-12, 200);
t1 = toc(start);
fprintf('%5.2e s : time for gmres\n',t1)    






L = max(abs(chnkrs.r),[],"all");

x1 = linspace(-L-0.15, L + 0.15 ,500);
[xx,yy] = meshgrid(x1,x1);
targets = [xx(:).'; yy(:).'];


skern = @(s,t) chnk.lap2d.kern(s, t, 's');



[~, nt] = size(targets);
true_sol = zeros(nt,1);
utarg = zeros(nt, 1);


nexttile
h = pcolor(xx,yy,nan(size(xx)));
set(h,'EdgeColor','None'); 
hold on;
plot(chnkrs,'k-','LineWidth',2);
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
