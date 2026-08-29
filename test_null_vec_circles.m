clear

% the style of genetrating geometry is from the demo_mult_connect.m

iseed = 8675309;
rng(iseed,'twister');

cparams = [];
cparams.eps = 1.0e-6;
cparams.nover = 0;

chnkr_int = [];
L = 1.5;                % domain size 
ntry = 1000;          
rad = 0.1;          
ninc = 16;              % number of inclusions

ctrs = [];    

for i = 1:ninc
    placed = false;    
    chnkr_i = chunkerfunc(@(t) rad .* circle(t), cparams);

    for j = 1:ntry
        tmp = (L - 2*rad) * rand(2,1) ;
        if isempty(ctrs)
            dist_ok = true;
        else
             dists = vecnorm(tmp - ctrs);
             dist_ok = all(dists > 2*rad); 
        end

        if dist_ok
           chnkr_i = chnkr_i.move([],tmp);
           chnkr_int = [chnkr_int, chnkr_i];
           ctrs =[ctrs, tmp];
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


lhs = 0.5.*eye(chnkrs.npt) + sysmat  ;     % left hand side

[rho, lambda_max] = eigs(lhs, 1, 'largestreal');

r2 = chnkrs.r(1,:).';
msure_wts = chunkerintegral(chnkrs, rho.*r2);






x1 = linspace(-0.2, 1.4 ,500);               % generate some targets to evaluate the eigenfunctions
[xx,yy] = meshgrid(x1,x1);
targets = [xx(:).'; yy(:).'];
                   

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
title('log scale Gamma operator eigenfunction', 'FontSize',12)
axis equal
colorbar
pbaspect([1 1 1]) 




set(gcf,'Position', [900  500  900   500])

exportgraphics(gcf,'null_vector_unordered.pdf','ContentType','vector',...
    'Resolution',600)
