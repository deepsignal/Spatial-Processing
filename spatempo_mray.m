
%% ----- Spatio_Temporal processing for multi-ray ------- %%

clear all
%close all


ray_proc_funs

% Generate OFDM 
%data = [zeros(5,1) ; ones(26,1) ; 1 ; ones(27,1) ; zeros(5,1)] ;
data = [1 ; ones(27,1) ; zeros(5,1) ; zeros(5,1) ; ones(26,1) ] ;


data_time = ifft(data,64);
data_pack = repmat(data_time,100,1);
data_pack_rx = resample(data_pack,16,1,50);


% Now to delay and also project on angular manifold, each ray separately..... 
phis = [80 70] ; %-90]  ; %-60 40];
dels = [16 16+8] ; % 64 64+16]; 
m = 2 ; %length(phis) ; %3;
d = 0.5 ; % antenna spacing
j=sqrt(-1);

ang_mani=exp(-2*pi*1i*d*[0:m-1].'*sin([phis]*pi/180));

% ---- Delay Params ------ %
data_pack_rx_mp_stack = [];
for i=1:length(dels)
    data_pack_rx_mp_stack(i,:) = data_pack_rx(1+dels(i):1024*4+dels(i)).' ;
end

data_rx_antenna = ang_mani*data_pack_rx_mp_stack;
figure ; plot(real(data_rx_antenna(1,:))) ; hold on ; plot(imag(data_rx_antenna(1,:)),'r')


% Reciever processing and symbol equalization
data_rx_dsp = resample(data_rx_antenna.',1,16);

symbs_per_ant = length(data_rx_dsp)/64 ; 
H_symb_marray = [];

% reshape each antenna rx, into individual symbols, and multiply with reference
% to generate the channel response
for i=1:size(data_rx_dsp,2)
    H_symb_marray(:,:,i) = fft(reshape(data_rx_dsp(1:end,i),[64 4]),64,1).*repmat(data,1,symbs_per_ant);
end


figure ; plot(fftshift(real(H_symb_marray(:,1,1)))); hold all
plot(fftshift(real(H_symb_marray(:,2,1))));
plot(fftshift(real(H_symb_marray(:,3,1))));


figure ; plot(fftshift(real(H_symb_marray(:,:,1)))); hold all

%equalized_symbs1 = (fft(data_rx_dsp(1:64,1))).*data ;  % conjugate data in reality, but only real in this case
%equalized_symbs2 = (fft(data_rx_dsp(65:128,1))).*data ; 
%equalized_symbs3 = (fft(data_rx_dsp(1:64,2))).*data ;
%equalized_symbs4 = (fft(data_rx_dsp(65:128,2))).*data ;


%% -- Pmusic based detection -- %%

% dimensions of m-array are [carriers : time : antennas] 
Xm = [fftshift(H_symb_marray(:,1,1))] ; % ; fftshift(equalized_symbs3); fftshift(equalized_symbs4)] ;
Xmr = corrmtx(Xm,40,'modified');
XmR = Xmr'*Xmr;

% % average over antennas
% if 1
%     XmR=zeros(40+1,40+1);
%     for i=1:2
%        Xm = [fftshift(H_symb_marray(:,1,i))] ; 
%        Xmr = corrmtx(Xm,40,'modified');
%        XmR = XmR + Xmr'*Xmr; 
%     end
% end

[P,f]=pmusic(XmR,2,1024*2,'corr');
figure(786) ; plot(log10(fftshift(abs(P)./max(abs(P))))); hold all ; grid on
title('Built-in music with correlation mat')

[P,f]=peig(XmR,2,1024*2,'corr');
figure(787) ; plot(log10(fftshift(abs(P)./max(abs(P))))); hold all ; grid on



% MUSIC with sub-array smoothing
num_subarrays = 40 ;
corr_lags = 64-num_subarrays ;
Xm_subarray_fwd = []; Xm_subarray_bwd = [];

for i=1:num_subarrays
    Xm_subarray_fwd(:,i)= [H_symb_marray(i:i+64-num_subarrays,1,1)] ; % ; H_symb_marray(i:i+64-num_subarrays,2,1)];
    Xm_subarray_fwd1(:,i)= [H_symb_marray(i:i+64-num_subarrays,1,2)] ; % ; H_symb_marray(i:i+64-num_subarrays,2,1)];
    %Xm_subarray_fwd2(:,i)= [H_symb_marray(i:i+64-num_subarrays,1,3)] ; % 
    Xm_subarray_bwd(:,i)= conj(H_symb_marray(i:i+64-num_subarrays,1,1));
end

Xm_subarray_fwd_bwd = [Xm_subarray_fwd ; Xm_subarray_bwd ];

R_subarray = Xm_subarray_fwd_bwd*Xm_subarray_fwd_bwd' ;
R_subarray1 = Xm_subarray_fwd1*Xm_subarray_fwd1' ;
%R_subarray2 = Xm_subarray_fwd2*Xm_subarray_fwd2' ;

%R_subarray = Xm_subarray_fwd_bwd*Xm_subarray_fwd_bwd' ;

R_subarray = R_subarray ; %+ R_subarray1 ; %+ R_subarray2;
% to use the forward-backward approach, uncomment the next line
%R_subarray=(R_subarray+fliplr(eye(length(R_subarray)))*R_subarray.'*fliplr(eye(length(R_subarray))))/2;

[P,f]=pmusic(R_subarray,2,1024*2,'corr');
figure(789) ; plot(20*log10(fftshift(abs(P)./max(abs(P))))); hold all ; grid on
title('Pmusic with subarray formation')
music_result = P;

%% --- Tensor decomposition --- %


Y = transpose([fftshift(H_symb_marray(:,1,1))  fftshift(H_symb_marray(:,1,2))]) ; % a symbol from two antennas

% form smoothed spatemp matrix

spatemp_mat = [];
subc_elements = 54 ;
subarrays = 64-subc_elements  ; % 4 elements across each symbol subcarriers

for i=1:subarrays
    spatemp_mat(:,i) = reshape(transpose(Y(:,i:i+subc_elements-1)),2*subc_elements,1);
end

% form spatio-temporal basis

fft_sz = 1024 ;
A = dftmtx(fft_sz);

fft_sz_step = 200 ; %fft_sz_stp/2;
A_xarray = A(1:subc_elements,1:fft_sz_step); % selecting three only to make a six element array

m = 2   ; L = 180 ; d = 0.5 ;
steering_mat = [];
k=1;



for i = 1 : L,
   %steer(i,:)=exp(-2*pi*1i*d*sin(-pi/2 + pi*(i-1)/L)*[0:m-1].');
   a_tau_theta_1 = exp(-2*pi*1i*d*sin(-pi/2 + pi*(i-1)/L)*[0]).*A_xarray;
   a_tau_theta_2 = exp(-2*pi*1i*d*sin(-pi/2 + pi*(i-1)/L)*[m-1]).*A_xarray;
   steering_mat(1:2*subc_elements,(i-1)*fft_sz_step+1:(i)*fft_sz_step) =  [a_tau_theta_1 ; a_tau_theta_2];
end


% corr matrix 
R_spatemp = spatemp_mat*spatemp_mat' ;

% pseudospectrum
%[V_Ravg D_Ravg] = eig(R_spatemp);
[V_Ravg D_Ravg U_Ravg] = svd(R_spatemp);
spatemp_projec = (sum(abs(transpose(steering_mat)*(V_Ravg(:,end-3:-1:end-80))),2)) ; % project on noise eigens
spatemp_spec = log10(1./(spatemp_projec)) ; %log10(1./abs(spatemp_projec)) ; 
spa_temp_image = reshape(spatemp_spec,fft_sz_step,180);

figure ; imagesc(spa_temp_image) ; colorbar

spa_temp_image_seg = sum(spa_temp_image([34:60],:).^1,1);
figure ; plot(fftshift(spa_temp_image_seg))


angular_spec = (fftshift(spa_temp_image_seg));
angular_spec_mean = mean(angular_spec);
angular_spec = angular_spec-angular_spec_mean; 
angles_below_mean = find(angular_spec<0);
angular_spec(angles_below_mean)=0;
angular_spec_disp =angular_spec ; 

figure(1034) ;  h1=polar(([-89:90].*pi/180),angular_spec_disp./max(angular_spec_disp)) ; set(h1,'linewidth',4)


%% ---- Doing steering mat (spatio-temporal) fit the sparse way ----- %%

%steering_mat
%spatemp_mat

fft_sz = 1024;
%A = dftmtx(fft_sz);
%A = (A(1:16,1:200)) ; %
                         
%Y= h_raw_allq_cell{sel_path}(label1_start:label1_start+5,27:42)';
%Y= h_raw_allq_cell{sel_path}(3500:3500+5,27:42)';

D_eigs = diag(D_Ravg); 
Y =  V_Ravg(:,1:15)*diag(D_eigs(1:15)); %  spatemp_mat  ; % V_Ravg(:,1:10) ;
A = steering_mat ;

% noise variance setup
[m,n] = size(A);
r = size(Y,2);
X= zeros(n,r);

sigma =norm(Y)/1e5;  s_max =4;

% Algorithm step
supp = [];
supp_c = setdiff(1:n,supp);
res = Y;
omp_iter = 1;

while length(supp) < s_max % norm(res(:)) > sigma %length(supp) < s_max  %&&  norm(res(:)) > sigma
    supp_c = setdiff(1:n,supp);
    obs = sum(abs(res'*A(:,supp_c)).^2,1);
    [mval,midx] = max(obs);
    supp = union(supp,supp_c(midx));
    [tmpU,R] = qr(A(:,supp),0);
    res= Y - 1.*(tmpU*tmpU'*Y);
    res_norm(omp_iter) = norm(res'*res);
    omp_iter = omp_iter+1;
end

supp = sort(supp,'ascend');
X(supp, :) = pinv(A(:,supp))*Y;

figure(710) ; imp_resp = sum(X(1:end,:),2) ; stem(abs(imp_resp)./max(abs(imp_resp))) ; hold all
title('Batch pursuits with spatiotemporal arrays')



spa_temp_image = reshape(imp_resp,fft_sz_step,180);
figure ; imagesc(abs(spa_temp_image)) ; colormap hot



