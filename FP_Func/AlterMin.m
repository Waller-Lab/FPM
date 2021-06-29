function [O, P, err, scale, Ns] = AlterMin( I, No, Ns, opts)
%AlterMin Implements alternating minimization sequentially on a stack of
%measurement I (n1 x n2 x nz). It consists of 2 loop. The main loop update
%the reconstruction results O and P. the inner loop applies projectors/minimizers
%P1 and P2 on each image I and steps through the entire dataset
%   Outputs:
%   O: reconsturcted high-resolution complex object
%   P: reconstructed complex pupil function
%   err: errors at each iteration
%   scale: LED brightness map
%   Ns: estimated LED positions using local search algorithms
%
%   Inputs:
% Measurements data
%   I: intensity measurements by different LEDs
% Reconstruction parameters
%   No = [Ny_obj,Nx_obj]: size of the reconstructed image
% Illumination coding parameters
%   Ns = [Nsy,Nsx]: centers of corresponding lpf regions for
%   the illumination pattern

% Iteration parameters: opts
%   tol: maximum change of error allowed in two consecutive iterations
    %   maxIter: maximum iterations 
    %   minIter: minimum iterations
    %   monotone (1, default): if monotone, error has to monotonically dropping
    %   when iters>minIter
%   display: display results (0: no (default) 1: yes)
    %   saveIterResult: save results at each step as images (0: no (default) 1: yes)
    %   mode: display in 'real' space or 'fourier' space.
    %   out_dir: saving directory
%   O0, P0: initial guesses for O and P
    %   OP_alpha: regularization parameter for O
    %   OP_beta: regularization parameter for P
%   scale: LED brightness map
%   H0: known portion of the aberration function, 
        % e.g. sample with a known defocus induce a quardratic aberration
        % function can be defined here
%   poscalibrate: flag for LED position correction using
    % '0': no correction
    % 'sa': simulated annealing method
        % calbratetol: parameter in controlling error tolence in sa
    % 'ga': genetic algorithm
    % caution: takes consierably much longer time to compute a single iteration
%   F, Ft: operators of Fourier transform and inverse

% Last modified on 10/07/2017
% by Lei Tian, lei_tian@alum.mit.edu

%% derived constants
% size of measurement
[Nmy,Nmx,Nimg] = size(I);
Np = [Nmy,Nmx];
% r0 defines # of LEDs lit up in each pattern
[r0,~,~] = size(Ns);
cen0 = round((No+1)/2);
row = @(x) x(:).';

%% options for the algorithms
if nargin<4
    % default values
    opts.tol = 1;
    opts.maxIter = 50;
    opts.minIter = 3;
    opts.monotone = 1;
    opts.display = 0;
    opts.saveIterResult = 0;
    opts.out_dir = [];
    opts.O0 = Ft(sqrt(I(:,:,1)))/r0;
    opts.O0 = padarray(opts.O0,(No-Np)/2);
    opts.P0 = ones(Np);
    opts.OP_alpha = 1;
    opts.OP_beta = 1;
    opts.mode = 'real';
    opts.scale = ones(Nled,1);
    opts.H0 = ones(Np);
    opts.poscalibrate = 0;
    opts.calbratetol = 1e-1;
    opts.F = @(x) fftshift(fft2(x));
    opts.Ft = @(x) ifft2(ifftshift(x));
else
    if ~isfield(opts,'tol')
        opts.tol = 1;
    end
    if ~isfield(opts,'maxIter')
        opts.maxIter = 50;
    end
    if ~isfield(opts,'minIter')
        opts.minIter = 3;
    end
    if ~isfield(opts,'monotone')
        opts.monotone = 1;
    end
    if ~isfield(opts,'display')
        opts.display = 0;
    end
    if ~isfield(opts,'saveIterResult')
        opts.saveIterResult = 0;
    end
    if ~isfield(opts,'out_dir')
        opts.out_dir = ['IterResults'];
        if opts.saveIterResult
            mkdir(opts.out_dir);
        end
    end
    if ~isfield(opts,'O0')
        opts.O0 = Ft(sqrt(I(:,:,1)))/r0;
        opts.O0 = padarray(opts.O0,(No-Np)/2);
    end
    if ~isfield(opts,'P0')
        opts.P0 = ones(Np);
    end
    if ~isfield(opts, 'StepSize')
        opts.StepSize = 1;
    end
    if ~isfield(opts,'OP_alpha')
        opts.OP_alpha = 1;
    end
    if ~isfield(opts,'OP_beta')
        opts.OP_beta = 1;
    end
    if ~isfield(opts,'mode')
        opts.mode = 'real';
    end
    if ~isfield(opts,'Ps')
        opts.Ps = 1;
    end
    if ~isfield(opts,'iters')
        opts.iters = 10;
    end
    if ~isfield(opts,'scale')
        opts.scale = ones(Nled,1);
    end

    if ~isfield(opts,'H0')
        opts.H0 = ones(Np);
    end
    if ~isfield(opts,'poscalibrate')
        opts.poscalibrate = 0;
    end
    if ~isfield(opts,'F')
        opts.F = @(x) fftshift(fft2(ifftshift(x)));
    end
    if ~isfield(opts,'Ft')
        opts.Ft = @(x) fftshift(ifft2(ifftshift(x)));
    end
    if ~isfield(opts,'calbratetol')
        opts.calbratetol = 1e-1;
    end
end

H0 = opts.H0;
F = opts.F;
Ft = opts.Ft;

%% Operators
Ps = opts.Ps;
% operator to crop region of O from proper location at the O plane
downsamp = @(x,cen) x(cen(1)-floor(Np(1)/2):cen(1)-floor(Np(1)/2)+Np(1)-1,...
    cen(2)-floor(Np(2)/2):cen(2)-floor(Np(2)/2)+Np(2)-1);

T0 = clock;

fprintf('| iter |  rmse    |\n');
for j=1:20, fprintf('-'); end
fprintf('\n');



%% initialization in FT domain
P = opts.P0; opts.P0 = 0;
O = opts.O0; opts.O0 = 0;
err1 = inf;
err2 = 50;
err = [];
iter = 0;
scale = opts.scale;
scale = reshape(scale,r0,Nimg);

if opts.display
    if strcmp(opts.mode,'real')
        o = O;
    elseif strcmp(opts.mode,'fourier')
        o = Ft(O);
    end
    f1 = figure(88);
    subplot(221); imagesc(abs(o)); axis image; colormap gray; colorbar;
    title('ampl(o)');
    subplot(222); imagesc(angle(o)); axis image; colormap gray; colorbar;
    title('phase(o)');
    subplot(223); imagesc(abs(P)); axis image; colormap gray; colorbar;
    title('ampl(P)');
    subplot(224); imagesc(angle(P).*abs(P)); axis image; colorbar;
    title('phase(P)');
    drawnow;
end

if opts.saveIterResult
    export_fig(f1,[opts.out_dir,'\R_',num2str(iter),'.png'],'-m4');
end


%% main algorithm starts here
% stopping criteria: when relative change in error falls below some value,
% can change this value to speed up the process by using a larger value but
% will trading off the reconstruction accuracy
% error is defined by the difference b/w the measurement and the estimated
% images in each iteration
fprintf('| %2d   | %.2e |\n',iter,err1);

sp0 = max(row(abs(Ns(:,1,:)-Ns(:,2,:))));

while abs(err1-err2)>opts.tol&&iter<opts.maxIter
%     psistack = zeros(64,64,293);
    err1 = err2;
    err2 = 0;
    iter = iter+1;   
    for m = 1:Nimg
        % initilize psi for correponing image, ROI determined by cen
        Psi0 = zeros(Np(1),Np(2),r0);
        Psi_scale = zeros(Np(1),Np(2),r0);
        cen = zeros(2,r0);
        scale0 = zeros(r0,1);
        for p = 1:r0
            cen(:,p) = cen0-row(Ns(p,m,:));
            scale0(p) = scale(p,m);
            Psi0(:,:,p) = downsamp(O,cen(:,p)).*P.*H0;
            Psi_scale(:,:,p) = sqrt(scale0(p))*Psi0(:,:,p);
        end
        % measured intensity
        I_mea = I(:,:,m);
        % compute field in real space
        psi0 = Ft(Psi_scale);
%         psistack(:,:,m) = psi0;
        % estimated intensity w/o correction term 
        I_est = sum(abs(psi0).^2,3);
        Psi = Proj_Fourier_v2(psi0, I_mea, I_est, scale0, F);
        % projection 2
        dPsi = Psi-Psi0;
        Omax = abs(O(cen0(1),cen0(2)));
        if r0 == 1
            P2 = @(O,P,dpsi,Omax,cen)...
                GDUpdate_Multiplication_rank1(O,P,dpsi,Omax,cen,Ps,...
                opts.OP_alpha,opts.OP_beta, opts.StepSize);
        else
            P2 = @(O,P,dpsi,Omax,cen)...
                GDUpdate_Multiplication_rank_r(O,P,dpsi,Omax,cen,Ps,...
                opts.OP_alpha,opts.OP_beta);
        end
        [O,P] = P2(O,P,dPsi./repmat(H0,[1,1,r0]),Omax,cen);
        
        %% position correction
        poscost = @(ss) sum(sum((abs(Ft(downsamp(O,ss).*P.*H0)).^2-I_mea).^2));
        if strcmp(opts.poscalibrate,'sa')
            optsanneal = saoptimset('Display','off','TolFun',opts.calbratetol);
            cen_correct = round(simulannealbnd(poscost,...
                cen(:,1),cen(:,1)-sp0/3,cen(:,1)+sp0/3,optsanneal));
            Ns(:,m,:) = cen0-cen_correct';
        elseif strcmp(opts.poscalibrate,'ga')
            optsanneal = saoptimset('Display','off');
            cen_correct = ga(poscost,2,[],[],[],[],...
                cen(:,1)-sp0/3,cen(:,1)+sp0/3,[],[1,2],optsanneal);
            Ns(:,m,:) = cen0-cen_correct;
        end
        % compute the total difference to determine stopping criterion
        err2 = err2+sqrt(sum(sum((I_mea-I_est).^2)));
        
    end
    if strcmp(opts.display,'full')
        if strcmp(opts.mode,'real')
            o = O;
        elseif strcmp(opts.mode,'fourier')
            o = Ft(O);
        end
        f1 = figure(88);
        subplot(221); imagesc(abs(o)); axis image; colormap gray; colorbar;
        title('ampl(o)');
        subplot(222); imagesc(angle(o)); axis image; colormap gray; colorbar;
        title('phase(o)');
        subplot(223); imagesc(abs(P)); axis image; colormap gray; colorbar;
        title('ampl(P)');
        subplot(224); imagesc(angle(P).*abs(P)); axis image; colorbar;
        title('phase(P)');
        drawnow;
    end
    %end
    
    %% compute error
    % record the error and can check the convergence later.
    err = [err,err2];
    
    if strcmp(opts.display,'iter')
        if strcmp(opts.mode,'real')
            o = O;
        elseif strcmp(opts.mode,'fourier')
            o = Ft(O);
        end
        f1 = figure(88);
        subplot(221); imagesc(abs(o)); axis image; colormap gray; colorbar;
        title('ampl(o)');
        subplot(222); imagesc(angle(o)); axis image; colormap gray; colorbar;
        title('phase(o)');
        subplot(223); imagesc(abs(P)); axis image; colormap gray; colorbar;
        title('ampl(P)');
        subplot(224); imagesc(angle(P).*abs(P)); axis image; colorbar;
        title('phase(P)');
        drawnow;
    end
    
    fprintf('| %2d   | %.2e |\n',iter,err2);
    
    if opts.saveIterResult
        export_fig(f1,[opts.out_dir,'\R_',num2str(iter),'.png'],'-m4');
        %     saveas(f2,[opts.out_dir,'\Ph_',num2str(iter),'.png']);
    end
    
    if opts.monotone&&iter>opts.minIter
        if err2>err1
            break;
        end
    end
    
    
end
if strcmp(opts.mode,'fourier')
    O = Ft(O);
end

fprintf('elapsed time: %.0f seconds\n',etime(clock,T0));

end

