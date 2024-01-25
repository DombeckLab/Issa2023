function [coeff,Bfull,Bhist,spk_est] = findbn(t,y,tau,neu,max_freq,exp_tau,plt,lam_est,f_est,max_iter)
% FINDBN estimtes y as a sum of a baseline function and a signal
%
% [coeff,Bfull,Bhist,spk_est] = FINDBN(t,y,tau,neu)
% [coeff,Bfull,Bhist,spk_est] = findbn(t,y,tau,neu,max_freq,exp_tau,plt,lam_est,f_est,max_iter)
%
% The full model is y = spk_est*h + B * coeff where h is a kernel defined
% by time constants tau, B is the basis functions, and coeff is the
% associated constants. Thus spk_est*h (convolution operation) is the
% estimate of F(t) and B*coeff (matrix multiplication) is the estimate of
% the baseline. Algorithm iterates between solving for coeff (using lsqlin)
% and updating spk_est (using oasisAR2).
%
% The full set of basis function is a sum of a constant function,
% sinusoids, exponentials, and a user-provided neuropil signal. The
% parameters of these functions are determined by inputs neu, max_freq, and
% exp_tau. Optimization of weights utilizes lsqlin and uses certain
% constraints.
%
% t: time vector
% y: signal
% tau: time constant of fluorescence signal (decay and rise time, seconds)
% neu: neuropil signal (optional). This gets incorporated into the basis
% max_freq: maximum sinusoid frequency in Hz (optional). 0 = none
% exp_tau: decay constants for exponentials (seconds). 0 = none
% plt: whether to show and update plot during fitting
% lam_est: scales sparsity penalty used by oasis (lower => more firing)
% f_est: initial estimate of firing f (spk_est*h)
% max_iter: max iterations to run (default: 25)
%
% coeff: N x 1, coefficients of basis functions
% Bfull: T x N, basis functions
% Bhist: N x max_iter, stores coeff for each iteration
% spk_est: T x 1, estimate of deconvolved signal
%
% t,y,neu are T-length vectors. row vectors are rescaled to column vectors
% N: # of basis functions, returned in Bfull. First column is constant
% function (1s), followed by sin/cos functions and then exponentials. N-1
% column is neuropil and last column (N) is estimated firing f_est (spk_est
% convolved with the kernel). coeff is the corresponding scaling factors.
%
% Example
%   [coeff,Bfull,Bhist,spk_est] = FINDBN(t,y,[.6 .1],neu,.002,[150 1000], ...
%       true,[.2 .2],[],25);
%   or
%   [coeff,Bfull,Bhist,spk_est] = FINDBN(t,y); % other inputs are optional
%
% The following variables may be of use:
%   Bnpil = Bfull(:,end-1)*coeff(end-1); % scaled neuropil
%   Bline = Bfull(:,1:end-2)*coeff(1:end-2); % scaled baseline
%   F_est = (Bfull(:,end)*coeff(end))./Bline; % F-F0 estimate
%   dF = (y-Bnpil-Bline)./Bline; % dF/F0 estimate
%   z_norm = spk_est./Bline; % scale the spike estimate

%% Author: John Issa, Northwestern University, 2019 (updated 2023)
%% References: Issa et.al., Nat Neuro 2024, Lateral entorhinal cortex subpopulations represent experiential epochs surrounding reward
%% available at https://www.nature.com/articles/s41593-023-01557-4
%%%%%% on March 13, 2022, copied oasisAR2 from CaImAn-MATLAB-master to replace oasisAR2 from zhoupc

%% preliminaries
if ~exist('tau','var') || isempty(tau), tau = [.6 .1]; end % fall and rise time, in seconds
if ~exist('neu','var') || isempty(neu), neu = []; else neu = double(neu); c = max(neu(:)); neu = neu/c; end % scale down neuropil
if ~exist('max_freq','var') || isempty(max_freq), max_freq = .001; end % maximum sinusoid frequency in Hz (.01)
if ~exist('exp_tau','var') || isempty(exp_tau), exp_tau = [250 1000]; end % decay constants for exponentials
if ~exist('plt','var') || isempty(plt), plt = 1; end % whether to plot results
if ~exist('lam_est','var') || isempty(lam_est), lam_est = [1 1]; end % scale of lam and smin
if ~exist('max_iter','var') || isempty(max_iter), max_iter = 25; end % max iterations to run
y = double(y(:)); t = double(t(:)); T = size(y,1); exp_tau = exp_tau(:)';

if ~exist('f_est','var') || isempty(f_est), f_est = zeros(T,1); end % initial estimate of f
f_est = double(f_est);

t0 = t(1);
t = t-t0;

cr = [.25 1.5]; % range for neuropil ratio. default: [.25 1.5]

dt = median(diff(t)); dr = exp(-dt./tau); g = [sum(dr) -prod(dr)]; % exp2ar(tau/dt)
min_freq = 1/(T*dt); % minimum frequency in Hz

%% initialize baseline parameters
if max_freq == 0
    fr = []; Bf = [];
else
    fr = min_freq:min_freq:max_freq; % fr = [6:-1:1]*.0025; % 6:-1:1
    Bf = [cos(2*pi*t*fr) sin(2*pi*t*fr)];
end
if exp_tau == 0
    Be = []; exp_tau = [];
else
    Be = [exp(-t./exp_tau) flipud(exp(-t./exp_tau))];
end

Bset = [ones(size(t)) Bf Be neu]; % basis set

Nf = 2*length(fr); % number of sinusoids
Ne = 2*length(exp_tau); % number of exponentials
Nn = size(neu,2); % number of neuropil components
NB = size(Bset,2); % Nc+Nf+Ne+Nn
i_neu = NB+1 - (Nn:-1:1); % index of neuropil

Bcoeff = zeros(NB,1);

Bhist = NaN(NB+1,max_iter); % stores Bcoeff for each iteration

%% set-up parameters for lsqlin (constrained linear least squares)
options = optimoptions(@lsqlin,'display','none');
lb = zeros(NB+1,1); lb(1) = min(y)/20;
lb(1+(1:Nf)) = -Inf; % let sinusoids be negative

ub = Inf(NB+1,1);
A = [-.1  ones(1,Nf/2) -ones(1,Nf/2) zeros(1,Ne) zeros(1,Nn) 0; ...
     -.1 -ones(1,Nf/2)  ones(1,Nf/2) zeros(1,Ne) zeros(1,Nn) 0];
Amax = [0; 0];
if isempty(neu)
    % do nothing
else
    A = [A; -1 zeros(1,Nf) zeros(1,Ne) .25*ones(1,Nn) 0]; Amax = [Amax; 0];
    lb(i_neu) = cr(1)*c; ub(i_neu) = cr(2)*c;
end

%% iterate between estimating basis function weights (est_basis) and deconvolving (est_spk)
est_basis
for i = 1:max_iter
    est_spk
    est_basis
    Bhist(:,i) = coeff;

    if plt
        plot(t+t0,y); hold on; plot(t+t0,Bset*Bcoeff,'linewidth',2); plot(t+t0,Bset*Bcoeff+f_est*coeff(end),'linewidth',1); hold off; box off
    end
    fprintf('.'); if mod(i,50)==0, fprintf('\n'); end
    if plt, drawnow; end

    % stop if solution has converged
    if i > 1
        if all(abs(Bhist(:,i)-Bhist(:,i-1))/abs(Bhist(:,i)) < 1e-3)
            fprintf('terminating early. %d of %d',i,max_iter)
            Bhist = Bhist(:,1:i);
            break
        end
    end
end
fprintf('\n');

if ~isempty(neu)
    coeff(i_neu) = coeff(i_neu)/c;
    Bfull(:,i_neu) = Bfull(:,i_neu)*c;
end

%%

    function est_spk
        %% estimate spikes
        y0 = y - Bset*Bcoeff;
        
        % try to estimate a semi-normalized lam and smin for oasis
        % [~,~,sig] = findbaseline(y0,1); % old, not used
        % [sig] = est_sig(y0);
        sig = std(y0-smoothdata(y0,1,'movmean',3)); % fast estimate of noise
        lam = 2*(sig/sqrt(abs(skewness(y0))))^1.25;
        smin = sig;
        
        lam = lam*lam_est(1); smin = smin*lam_est(2);

        if length(tau) == 2
            [f_est, spk_est] = oasisAR2(y0,g,lam,smin);

            % if oasisAR2 fails, try something else
            if min(y - f_est) < -sig*20
                disp(['alternate: ' num2str([min(y0 - f_est)   -sig*20])])
                [f_est, spk_est] = thresholded_oasisAR2(y0, g, sig);
            end
        else
            [f_est, spk_est] = oasisAR1(y0,g(1),lam,smin);
        end
    end

    function est_basis
        %% estimate basis set
        Bfull = [Bset f_est];
        coeff = lsqlin(Bfull,y,A,Amax,[],[],lb,ub,[],options);
        if isempty(coeff), coeff = zeros(NB+1,1); end
        Bcoeff = coeff(1:NB);
    end

    % function [sig] = est_sig(x)
    %     % estimate noise
    %     x2 = sort(x);
    %     [mu,sig] = normfit(x2(1:round(length(x)*.05))); % find lowest 5% of data
    %     x = x(find(x<=mu+8*sig));
    % 
    %     for j = 1:50
    %         [mu,sig]=normfit(x);
    %         x = x(find(x<=mu+1.5*sig));
    %     end
    % 
    %     if sig==0 % for the odd case where nothing is found
    %         [mu,sig]=normfit(x2);
    %         x = x2(find(x2<=mu+1.5*sig));
    %     end
    %     sig = std(x);
    % end


end
