function [u, u_hat, omega] = MVMD(signal, alpha, tau, K, DC, init, tol)
% Multivariate Variational Mode Decomposition
%
% The function MVMD applies the "Multivariate Variational Mode Decomposition (MVMD)" algorithm to multivariate or multichannel data sets. 
% We have verified this code through simulations involving synthetic and real world data sets containing 2-16 channels. 
% However, there is no reason that it shouldn't work for data with more than 16 channels.
% 
%
% Input and Parameters:
% ---------------------
% signal  - input multivariate signal that needs to be decomposed
% alpha   - the parameter that defines the bandwidth of extracted modes (low value of alpha yields higher bandwidth)
% tau     - time-step of the dual ascent ( pick 0 for noise-slack )
% K       - the number of modes to be recovered
% DC      - true if the first mode is put and kept at DC (0-freq)
% init    - 0 = all omegas start at 0
%         - 1 = all omegas start uniformly distributed
%         - 2 = all omegas initialized randomly
% tol     - tolerance value for convergence of ADMM
%
%
% Output:
% ----------------------
% u       - the collection of decomposed modes
% u_hat   - spectra of the modes
% omega   - estimated mode center-frequencies
%
%
% Syntax:
% -----------------------
% [u, u_hat, omega] = MVMD(X, alpha, tau, K, DC, init, tol)
%   returns:
%			 a 3D matrix 'u(K,L,C)' containing K multivariate modes, each with 'C' number of channels and length 'L', that are 
%            computed by applying the MVMD algorithm on the C-variate signal (time-series) X of length L.
%    		 - To extract a particular channel 'c' corresponding to all extracted modes, you can use u_c = u(:,:,c).
%			 - To extract a particular mode 'k' corresponding to all channels, you can use u_k = u(k,:,:).
%			 - To extract a particular mode 'k' corresponding to the channel 'c', you can use u_kc = u(k,:,c).

%			 3D matrix 'u_hat(K,L,C)' containing K multivariate modes, each with 'C' number of channels and length 'L', that  
%            are computed by applying the MVMD algorithm on the C-variate signal (time-series) X of length L.

%			 2D matrix 'omega(N,K)' estimated mode center frequencies

% Usage:
% -----------------------
% 	Example 1: Mode Alignment on Synthetic Data
% 	T = 1000; t = (1:T)/T;
% 	f_channel1 = (cos(2*pi*2*t)) + (1/16*(cos(2*pi*36*t))); % Channel 1 contains 2 Hz and 36Hz tones
% 	f_channel2 = (1/4*(cos(2*pi*24*t))) + (1/16*(cos(2*pi*36*t))); % Channel 2 contains 24 Hz and 36Hz tones
% 	f = [f_channel1;f_channel2]; % Making a bivariate signal
% 	[u, u_hat, omega] = MVMD(f, 2000, 0, 3, 0, 1, 1e-7); 

% 	Example 2: Real World Data (EEG Data)
% 	load('EEG_data.mat');
% 	[u, u_hat, omega] = MVMD(data, 2000, 0, 6, 0, 1, 1e-7);


% 	Authors: Naveed ur Rehman and Hania Aftab
% 	Contact Email: naveed.rehman@comsats.edu.pk
%
% 	Acknowledgments: The MVMD code has been developed by modifying the univariate variational mode decomposition code that has 
%                 been made public at the following link:
%                 https://www.mathworks.com/matlabcentral/fileexchange/44765-variational-mode-decomposition
%                 by K. Dragomiretskiy, D. Zosso.
%                 
%
% 	Please cite the following papers if you use this code in your work:
%   -----------------------------------------------------------------
% 
%  [1] N. Rehman, H. Aftab, Multivariate Variational Mode Decomposition, arXiv:1907.04509, 2019. 
%  [2] K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Transactions on Signal Processing, vol. 62, pp. 531-544, 2014. 


%---------- Check for Input Signal              
% Check for getting number of channels from input signal
[x, y] = size(signal);
if x > y
	C = y;
	signal = signal';
else
	C = x;
end

%---------- Preparations
for c = 1:C
    signal_c = signal(c,:); 
    save_T = length(signal_c);
    fs = 1/save_T;
	
    % Mirroring
    T = save_T;
    f_mirror(1:T/2) = signal_c(T/2:-1:1);
    f_mirror(T/2+1:3*T/2) = signal_c;
    f_mirror(3*T/2+1:2*T) = signal_c(T:-1:T/2+1);
    f(c,:) = f_mirror;
	
    % Time Domain 0 to T (of mirrored signal)
    T = length(f);
    t = (1:T)/T;
	
    % frequencies
    freqs = t-0.5-1/T;
	
	% Construct and center f_hat
    f_hat(c,:) = fftshift((fft(f(c,:))));
    f_hat_plus(c,:) = f_hat(c,:);
    f_hat_plus(c,1:T/2) = 0;
end

%------------ Initialization
% Maximum number of iterations 
N = 500;
% For future generalizations: individual alpha for each mode
Alpha = alpha*ones(1,K);
% matrix keeping track of every iterant 
u_hat_plus = zeros(N, length(freqs), K, C);
omega_plus = zeros(N, K);

% initialize omegas uniformly
switch init
	case 1
        for i = 1:K
            omega_plus(1,i) = (0.5/K)*(i-1);
        end
    case 2
            omega_plus(1,:) = sort(exp(log(fs) + (log(0.5)-log(fs))*rand(1,K)));
    otherwise
            omega_plus(1,:) = 0;
end

% if DC mode imposed, set its omega to 0
if DC
    omega_plus(1,1) = 0;
end

% start with empty dual variables
lambda_hat = zeros(N, length(freqs), C); 

% other inits
uDiff = tol+eps; % update step
n = 1; % loop counter
sum_uk = zeros(C, length(freqs)); % accumulator

%--------------- Algorithm of MVMD

while ( uDiff > tol &&  n < N ) % not converged and below iterations limit
    % update first mode
    k = 1;
    for c = 1:C
		% update first mode accumulator
        sum_uk(c,:) = u_hat_plus(n,:,K,c) + sum_uk(c,:) - u_hat_plus(n,:,1,c);
        % update spectrum of first mode through Wiener filter of residuals
        u_hat_plus(n+1,:,k,c) = (f_hat_plus(c,:) - sum_uk(c,:) - lambda_hat(n,:,c)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
	end	
	
	% update first omega if not held at 0
    if ~DC
		temp1 = 0;
		temp2 = 0;
		for c = 1:C
            numerator = freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k,c)).^2)';
            denominator = sum(abs(u_hat_plus(n+1,T/2+1:T,k,c)).^2);
            temp1 = numerator + temp1;
            temp2 = denominator + temp2;
		end
		%center frequency
        omega_plus(n+1,k) = temp1/temp2;
    end
		
	% update of any other mode
	for k=2:K
		for c = 1:C
			% accumulator
			sum_uk(c,:) = u_hat_plus(n+1,:,k-1,c) + sum_uk(c,:) - u_hat_plus(n,:,k,c);
			% mode spectrum
			u_hat_plus(n+1,:,k,c) = (f_hat_plus(c,:) - sum_uk(c,:) - lambda_hat(n,:,c)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
		end
		% center frequencies
		temp1=0;
		temp2=0;
		for c = 1:C
			numerator = freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k,c)).^2)';
			denominator = sum(abs(u_hat_plus(n+1,T/2+1:T,k,c)).^2);
			temp1 = numerator + temp1;
			temp2 = denominator + temp2;
		end
		omega_plus(n+1,k) = temp1/temp2; 
	end
    
	% Dual ascent
	for c = 1:C
		lambda_hat(n+1,:,c) = lambda_hat(n,:,c) + tau*(sum(u_hat_plus(n+1,:,:,c),3) - f_hat_plus(c,:));
    end
	
	% loop counter
	n = n+1;
    
	% converged yet?
	uDiff = eps;
	for i=1:K
		for c = 1:C
			uDiff = uDiff + 1/T*(u_hat_plus(n,:,i,c)-u_hat_plus(n-1,:,i,c))*conj((u_hat_plus(n,:,i,c)-u_hat_plus(n-1,:,i,c)))';
		end
	end
	uDiff = abs(uDiff);
end

%------ Post-processing and cleanup
% discard empty space if converged early
N = min(N,n);
omega = omega_plus(1:N,:);

% Signal reconstruction
u_hat = zeros(T, K, C);
for c = 1:C
	u_hat((T/2+1):T,:,c) = squeeze(u_hat_plus(N,(T/2+1):T,:,c));
	u_hat((T/2+1):-1:2,:,c) = squeeze(conj(u_hat_plus(N,(T/2+1):T,:,c)));
	u_hat(1,:,c) = conj(u_hat(end,:,c));
end

u = zeros(K,length(t),C);
for k = 1:K
	for c = 1:C
		u(k,:,c)=real(ifft(ifftshift(u_hat(:,k,c))));
	end
end

% remove mirror part
u = u(:,T/4+1:3*T/4,:);

% recompute spectrum
clear u_hat;
for k = 1:K
	for c = 1:C
		u_hat(:,k,c)=fftshift(fft(u(k,:,c)))';
	end
end

u_hat = permute(u_hat, [2 1 3]);
end