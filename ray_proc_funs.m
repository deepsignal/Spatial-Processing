


%% For ray processing using multi-antenna, multi-carrier data. Extracting 
%s super-resolved time delays and AoAs from channel data


function message = ray_proc_funs 

    assignin('base','STO_sanitize',@STO_sanitize)
    assignin('base','H2h_subspace',@H2h_subspace)
	assignin('base','spat_spec_SS',@spat_spec_SS)
    assignin('base','simulpursuits',@simulpursuits)
    
    
    message = 'Done importing functions to workspace' ;
    
end


function STO_sanitized_out = STO_sanitize(ch_resp_mat)
% will process complex channel response to output one
% with STO removed to zero
    if 1
        H_pair1  = ch_resp_mat ;
        phi_ant1 = unwrap(angle(H_pair1));    
        nf_del = 1:length(phi_ant1);
        phi_ant =  [-phi_ant1];
        nf_del = nf_del(:) ;  phi_ant = phi_ant(:);
        A_phi = [nf_del ones(length(phi_ant),1)];
        b_phi = phi_ant ; 
        phi_line = A_phi\b_phi; % solving for phase slope and intercept
        phi_rec = A_phi*phi_line;
        phi_ant_rec = phi_ant1 + transpose(phi_rec) ;
        STO_sanitized_out = abs(H_pair1).*exp(1i*phi_ant_rec);
        %figure(255+1) ; plot(phi_ant_rec(1:16)) ; hold on ; % plot(phi_ant_rec(17:32),'r')
        %fprintf('Phase slop computed is : Slope = %f Linear \n',phi_line(1))
    else
        disp('This has multiple antennas')
    end
end


function h_abs = H2h_subspace(Hcmplx,corr_len,fft_sz,spikes)
    % Testing with explicit corr_mat & explicit fourier proj
    %corr_len ->2nd arg
    %fft_sz -> 3rd arg
    Xm = corrmtx(Hcmplx,corr_len,'mod');
    XmR = Xm'*Xm;
    %XmR_in = (XmR' + XmR)/2 ; 
    %[~,d,v] = svd(XmR,0);
    [E,D] = eig((XmR+XmR')/2); % Ensure Hermitian
    [eigenvals,indx] = sort(diag(D),'descend');
    eigenvects = E(:,indx);
    %N =length(XmR);  % FFT matrix formation starts
    %Nf = fft_sz ; FN=[];
    %for n=1:Nf
    %    FN(n,:)=exp(-2*pi*1i*((n*N/Nf)-1)*(0:(N-1))/N)/sqrt(N);
    %end
    %pmusic_ex=FN*v(:,8);
    Pgoet=0;
    for j=spikes+1:corr_len %M-p Loop operating on noise vectors
       h = fft(eigenvects(:,j),fft_sz);
       Pgoet=Pgoet + abs(h.^2);
    end;
    %figure ; plot(20*log10(fftshift(1./Pgoet))); hold all
    h_abs = ((1./Pgoet));

end




function h_abs = LSp2SW(A,b)
% Squared LS stage wise regreassion: LS within the atom finding stage
% and LS on the overall atoms to determine the residual
        % Apply LS-MP
    thrLSp2SW=1e-4;
    r=b;
    SS=[];
    szA =  size(A);
    m = szA(2);
    while r'*r>thrLSp2SW,
        Z=zeros(m,1);
        for jj=1:1:m
            SStemp=[SS,jj]; 
            rtemp=b-A(:,SStemp)*pinv(A(:,SStemp))*b;
            Z(jj)=rtemp'*rtemp;
        end;
        posZ=find(Z==min(Z),1);
        SS=sort([SS,posZ(1)]);
        r=b-A(:,SS)*pinv(A(:,SS))*b;    
    end;
    h_abs=zeros(m,1);
    h_abs(SS)=pinv(A(:,SS))*b;

end


function [theta_spec_Ravg eff_packs] = spat_spec_SS(q_deck,h_raw_allq,ant_pair)
% Spatial spectrum generation using q-deck and h_raw_allq
% using sub-space decomposition method 

    %ant_pair = [1 2]; 
    [q_valid_inds1 v_bools1] = find(q_deck{ant_pair(1)}(:)>10);
    [q_valid_inds2 v_bools2] = find(q_deck{ant_pair(2)}(:)>10);

    C = intersect(q_valid_inds1,q_valid_inds2) ;
    A1 = h_raw_allq{ant_pair(1)}(C,25:25+16) ;
    A2 = h_raw_allq{ant_pair(2)}(C,25:25+16) ;

     for i=1:length(C)
        H_pair1_sto_comp_A1(i,:) = STO_sanitize(A1(i,:));
        H_pair1_sto_comp_A2(i,:) = STO_sanitize(A2(i,:));
     end

    % generate temporo_spatial basis
    A = dftmtx(128);
    A_xarray = A(1:3,:); % selecting three only to make a six element array
    m = 2   ;
    L = 180 ;
    d = 0.5 ;
    steer = zeros(180,2);
    a_tau_theta = zeros(6,180*128);
    k=1;

    for i = 1 : L,
       steer(i,:)=exp(-2*pi*1i*d*sin(-pi/2 + pi*(i-1)/L)*[0:m-1].');
       a_tau_theta(1:3,k:k+128-1) = exp(-2*pi*1i*d*sin(-pi/2 + pi*(i-1)/L)*[0]).*A_xarray;
       a_tau_theta(4:6,k:k+128-1) = exp(-2*pi*1i*d*sin(-pi/2 + pi*(i-1)/L)*[m-1]).*A_xarray;
       %phi(i)=real(a'*R*a);
       %phi(i)=1/(a'*inv(R)*a);
       k=k+128;
    end

    theta_spec_hist =[];
    R_tau_ang_avg = zeros(6,6);

    for p=1:length(C)

        H_concat = [H_pair1_sto_comp_A1(p,:) ; H_pair1_sto_comp_A2(p,:)];
        H_concat_mat = zeros(6,13);

        for i=1:16-3
            H_concat_mat(:,i) = reshape(transpose(H_concat(:,i:i+3-1)),6,1);
        end

        R_tau_ang = H_concat_mat*H_concat_mat';
        [V_R D_R] = eig(R_tau_ang);

        R_tau_ang_avg = R_tau_ang_avg + R_tau_ang ;

        tau_theta_spec = abs(a_tau_theta'*V_R(:,1:3));
        tau_theta_spec = sum(tau_theta_spec,2);
        %figure(641) ; plot(tau_theta_spec) ; hold all

        theta_spec = reshape(tau_theta_spec,128,180);
        theta_spec_avg = sum(theta_spec,1);
        %figure(741) ; plot(fftshift(theta_spec_avg)); hold all

        theta_spec_hist(:,p)=theta_spec_avg;
    end
    
    R_tau_ang_avg = R_tau_ang_avg ./ length(C);
    [V_Ravg D_Ravg] = eig(R_tau_ang_avg);
    tau_theta_spec_Ravg = log10(1./abs(a_tau_theta'*V_Ravg(:,1:2))); % we have to go inverse on noise eigens
    %tau_theta_spec_Ravg = abs(a_tau_theta'*V_Ravg(:,end:-1:end-3)); % straight on signal eigens
    tau_theta_spec_Ravg = sum(tau_theta_spec_Ravg,2);
    
    theta_spec_Ravg = reshape(tau_theta_spec_Ravg,128,180);
    theta_spec_Ravg = sum(theta_spec_Ravg,1);
    
    eff_packs = length(C);

end



%% --- Simultaneous pursuits ---- %%

function [pulses psupport] = simulpursuits(Y,fft_sz,basispan,numbasis,numspikes)

    %fft_sz = 1024;
    %basispan = 16 ; numbasis = 200 ;

    A = dftmtx(fft_sz);
    A = (A(1:basispan,1:numbasis)) ; %(fft_sz/2)-1)); % 64 time samples of 1024 different frequency shifts
                              % we are all samples here
    %Y = [h_raw_allq_cell{sel_path}(label1_start:label1_start+5,27:42)'];

    realY = real(Y) - repmat(mean(real(Y)),16,1) ;
    imagY = imag(Y) - repmat(mean(imag(Y)),16,1) ;

    Y = realY + 1i.*imagY ;


    % noise variance setup
    [m,n] = size(A);
    r = size(Y,2);
    X= zeros(n,r);

    sigma =norm(Y)/1e5;  s_max = numspikes;

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

    %figure(710) ; stem(abs(sum(X(1:end,:),2))) ; hold all
    %title('Batch pursuits result')

    pulses = sum(abs(X),2) ;
    psupport = supp ;
    
end









