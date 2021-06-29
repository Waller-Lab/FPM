%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main file to implement Fourier Ptychography reconstruction algorithm
% ref
% Lei Tian, et.al, Biomedical Optics Express 5, 2376-2389 (2014).
%
% last modified on 10/07/2015
% by Lei Tian, lei_tian@alum.mit.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% To do list for the user: (marked by 'TODO#')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) specify where the data located in 'filedir'
% 2) specify where you want to store the results in 'out_dir'
% 3) Find coordinates for estimating background levels and input into 'bk1'
% and 'bk2'.
% 4) specify a threshold value, above which the estimated background value
% is signal rather than noise.
% 5) make sure the LED index (used for taking the images) are properly defined
% in 'lit'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Reconstruction library locates here
%clear all;
addpath(['./FP_Func']);
addpath('natsortfiles');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO 1: specify the file directory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Multiplex image directory
filedir = ['../data/'];
% Generate the image list, in 'tif' image format (depending on your image format)
imglist = dir([filedir,'*.tif']);
nstart = [100, 100];
N = natsortfiles({imglist.name});%sorting the images in 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO 2: specify output folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%out_dir = ['/home/'];
%mkdir(out_dir);

%% define # of LEDs used to capture each image
numlit = 1;
% raw image size
n1 = 2160; n2 = 2560;


%% read in all images into the memory first
fprintf(['loading the images...\n']);
tic;
Nimg = length(imglist);
Iall = zeros(n1,n2,Nimg,'uint16');
Ibk = zeros(Nimg,1);
for m = 1:Nimg

    %     fn = [filedir,imglist(m).name]; %removed for sorting purpose
    fn = [filedir,N{m}];  %added for sorting purpose
    disp(fn);
    % all image data
    I = double(imread(fn));
    Iall(:,:,m) = I  ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% TODO 3: specify background region coordinates
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % backgr4ound noise esimtation 
    bk1 = mean2(double(Iall(1:100,1:100,m)));
    bk2 = mean2(double(Iall(492:600,2380:2520,m)));

    Ibk(m) = mean([bk1,bk2]);
    % brightfield background processing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO 4: if Ibk is larger than some threshold, it is not noise! (chnage 
    % this value correpondingly)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if Ibk(m)>300*(2^8-1)/(2^16-1)&& m>1
%         Ibk(m) = Ibk(m-1);
%     end
    if Ibk(m)>300 
        Ibk(m) = Ibk(m-1);
    end
end
fprintf(['\nfinish loading images\n']);
toc;

%% define processing ROI
Np = 344;
% nstart = [1,1];


%% read system parameters
% USAF_Parameter();
U2OS_Parameter();
%% load in data: read in the patch from the memory
Imea = double(Iall(nstart(1):nstart(1)+Np-1,nstart(2):nstart(2)+Np-1,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO 5: Define the LED index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% in this example, the LEDs are scanned sequentially within illumination NA
% defined in the systemsetup();
% In other cases, make the correponding changes. In the end, all we need is
% the matrix 'lit' that specifies the LED index corresponding image.

ledidx = 1:Nled;
ledidx = reshape(ledidx,numlit,Nimg);
lit = Litidx(ledidx);
lit = reshape(lit,numlit,Nimg);
% reorder LED indices based on illumination NA
[dis_lit2,idx_led] = sort(reshape(illumination_na_used,1,Nled));

Nsh_lit = zeros(numlit,Nimg);
Nsv_lit = zeros(numlit,Nimg);

for m = 1:Nimg
    % corresponding index of spatial freq for the LEDs are lit
    lit0 = lit(:,m);
    Nsh_lit(:,m) = idx_u(lit0);
    Nsv_lit(:,m) = idx_v(lit0);
end

% reorder the LED indices and intensity measurements according the previous
% dis_lit
Ns = [];
Ns(:,:,1) = Nsv_lit;
Ns(:,:,2) = Nsh_lit;

Imea_reorder = Imea(:,:,idx_led);
Ibk_reorder = Ibk(idx_led);
% gain = a(idx_led);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processing the data to DENOISING is IMPORTANT
% background subtraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ithresh_reorder = Imea_reorder;
for m = 1:Nimg
    Itmp = Ithresh_reorder(:,:,m);
    Itmp = Itmp-Ibk_reorder(m);
%     Itmp = awgn(Itmp,0,'measured');
    Itmp(Itmp<0) = 0;
    Ithresh_reorder(:,:,m) = Itmp;
    
end

%load('C:\Users\Muneeb\Desktop\ajmal fpm\laura\multiplexed fpm\Ns_cal289.mat\');
% Ns_reorder = Ns(:,idx_led,:);
% Ithresh_reorder = Ithresh_reorder(:,:,1:89);
Ns_reorder = Ns(:,idx_led,:);
clear Imea
%% reconstruction algorithm
% select the index of images that will be used in the processing
Nused = 293;
idx_used = 1:Nused;
I = Ithresh_reorder(:,:,idx_used);
Ns2 = Ns_reorder(:,idx_used,:);

%% reconstruction algorithm options: opts
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
opts.tol = 1;
opts.maxIter = 10; 
opts.minIter = 2;
opts.monotone = 1;
% 'full', display every subroutin,
% 'iter', display only results from outer loop
% 0, no display
opts.display = 'full';%0;%'iter';
upsamp = @(x) padarray(x,[(N_obj-Np)/2,(N_obj-Np)/2]);
opts.O0 = F(sqrt(I(:,:,1)));
opts.O0 = upsamp(opts.O0);
opts.P0 = w_NA;
opts.Ps = w_NA;
opts.iters = 1;
opts.mode = 'fourier';
opts.scale = ones(Nused,1);
opts.OP_alpha = 1;
opts.OP_beta = 1e3 ;
opts.poscalibrate =0;
opts.calbratetol = 1e-1;
opts.F = F;
opts.Ft = Ft;
opts.StepSize = 0.1;
%% algorithm starts
[O,P,err_pc,c,Ns_cal] = AlterMin(I,[N_obj,N_obj],round(Ns2),opts);

%% save results
fn = ['RandLit-',num2str(numlit),'-',num2str(Nused)];
%save([out_dir,'\',fn],'O','P','err_pc','c','Ns_cal');

%f1 = figure; imagesc(-angle(O),[-.6,1]); axis image; colormap gray; axis off

fprintf('processing completes\n');

%I = mat2gray(real(O));
%figure(2);imshow(I);

figure(2);imshow(angle(O),[]);
figure(1);imshow(abs(O),[])
% figure(3);imagesc(-angle(O));colormap gray;
