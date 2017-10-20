

% -- Motion processing function calls to generate motion markers -- %



function message = mot_proc_funs 

    assignin('base','vecspace_motion',@vecspace_motion)
    
    message = 'Done importing motion functions to workspace' ;
    
end



function final_motion_out = vecspace_motion(h_raw_allq_cell,ch_q_cell,st_montime,paths2use)

% catch the history for final combining, in actual implementation, happening inline 
var_vec_hist  = zeros(paths2use,length(h_raw_allq_cell{1}));
conv_vec_hist = zeros(paths2use,length(h_raw_allq_cell{1}));

% Motion vector generation
for m = 1:2

    test_path = m ;

    % Prinicpal extraction using progression approximation subspace.

    %test_path = 6;
    ch_sel = [ h_raw_allq_cell{test_path}(:,25:25+16-1) ] ; % h_raw_cell{test_path+2}(:,25:25+16-1) ; h_raw_cell{test_path-1}(:,25:25+16-1) ]  ; 

    lambda = 0.95 ;
    R = 0.01.*eye(16);
    mean_eig = 0 ;
    var_eig = 0.00;
    lambda2 = 0.99 ;
    R_eig =[] ; var_eig_h=[]; mean_eig_h=[];

    W     = eye(16); 
    p_wts = 0.01.*ones(16,1);
    %beta = [0.99 0.99 0.99 0.99 0.90 0.99 0.99 0.99 0.99 0.99 0.99 0.99 0.99] ;
    beta = 0.99.*(ones(1,13));
    eig_sorted = []; pre_proj = []; med_queue = zeros(1,10);
    Reig_buff = zeros(1,30) ; 

    nbound = 6 ; % the sigma variation beyond which to trigger an event
    rej_flag = zeros(1,length(ch_sel));
    proj_vec_dminus1s = 1e-4.*ones(1,10);
    proj_vec_dminus1f = 1e-4.*ones(1,10);
    
    var_vec = 1e-4.*ones(1,10); 
    vpcktcnt = 0 ;
    conv_count = 0;
    conv_flag = 0;
    

    for k = 2:length(ch_sel)

        %projpre = ch_sel(k,:)- (ch_sel(k,:)*W(:,1:3)*transpose(W(:,1:3)));
        ch_sel1 = transpose(ch_sel(k,:))./(max(abs(ch_sel(k,:)))*1);
        ch_sel2 = ch_sel1 - (W(:,1)'*ch_sel1).*W(:,1);
        ch_sel3 = ch_sel2 - (W(:,2)'*ch_sel2).*W(:,2);
        ch_sel4 = ch_sel3 - (W(:,3)'*ch_sel3).*W(:,3);
        ch_sel5 = ch_sel4 - (W(:,4)'*ch_sel4).*W(:,4);
        ch_sel6 = ch_sel5 - (W(:,5)'*ch_sel5).*W(:,5);
        ch_sel7 = ch_sel6 - (W(:,6)'*ch_sel6).*W(:,6);
        ch_sel8 = ch_sel7 - (W(:,7)'*ch_sel7).*W(:,7);
        ch_sel9 = ch_sel8 - (W(:,8)'*ch_sel8).*W(:,8);
        ch_sel10 = ch_sel9 - (W(:,9)'*ch_sel9).*W(:,9); 

        pre_proj(k,1)= (abs(W(:,1)'*ch_sel1)^2); 
        pre_proj(k,2)= (abs(W(:,2)'*ch_sel2)^2);
        pre_proj(k,3)= (abs(W(:,3)'*ch_sel3)^2);
        pre_proj(k,4)= (abs(W(:,4)'*ch_sel4)^2);
        pre_proj(k,5)= (abs(W(:,5)'*ch_sel5)^2);
        pre_proj(k,6)= (abs(W(:,6)'*ch_sel6)^2);
        pre_proj(k,7)= (abs(W(:,7)'*ch_sel7)^2);
        pre_proj(k,8)= (abs(W(:,8)'*ch_sel8)^2);
        pre_proj(k,9)= (abs(W(:,9)'*ch_sel9)^2);
        pre_proj(k,10)= (abs(W(:,10)'*ch_sel10)^2);

        int_index = 3; % specified how many layers under interference to look
        
        %pre_proj(k,2)= (abs(W(:,5)'*ch_sel5)^2); %p_wts(5) ; %(abs(W(:,5)'*ch_sel5)^2);
        if ch_q_cell{m}(k)<14
            mean_level(k)  = proj_vec_dminus1s(int_index);
            disturb_metr(k)= proj_vec_dminus1f(int_index);
            mean_levelone(k)    = proj_vec_dminus1s(1);
            disturb_levelone(k) = pre_proj(k,1);
            var_vec_hist(m,k) = var_vec_hist(m,k-1);
            conv_vec_hist(m,k) = conv_vec_hist(m,k-1);
            continue;
        end
        
        % increment the valid packet count
        vpcktcnt = vpcktcnt + 1 ;
   
        proj_vec = [pre_proj(k,1:10)];

        if vpcktcnt>100
    %         % exit in instances of gross interference
    %         if k>500 && pre_proj(k,1)<(0.8*(proj_vec_dminus1s(1)))
    %            rej_flag(k)=1;
    %            break; 
    %         end



            % Fast adaptation common to the entire vector 
            proj_lambda_f = 0.95;
            proj_vec_dminus1f = proj_lambda_f*proj_vec_dminus1f + (1-proj_lambda_f).*proj_vec ;

            th_set = 3 ; % determines how often we lock the slow filter 

            % Slow adaptation is differential across test and track objectives 
            proj_lambda_s = 0.999;
            proj_vec_dminus1s(1) = proj_lambda_s*proj_vec_dminus1s(1) + (1-proj_lambda_s).*proj_vec(1) ;
            % similarly others in the series below

            if (proj_vec_dminus1f(int_index) < th_set*proj_vec_dminus1s(int_index)) || (vpcktcnt<300) % to allow the adaptation of slow one
                proj_lambda_s = 0.995 ; %0.995;
                proj_vec_dminus1s(int_index) = proj_lambda_s*proj_vec_dminus1s(int_index) + (1-proj_lambda_s).*proj_vec(int_index) ;
            else
                proj_vec_dminus1s(int_index) = proj_vec_dminus1s(int_index) ; 
            end

            % can aggregate two or more int_indexes here below
            mean_level(k)  = proj_vec_dminus1s(int_index);
            disturb_metr(k)= proj_vec_dminus1f(int_index);


            % ++ Accumulating the primary projection metrics ++ %
            mean_levelone(k)    = proj_vec_dminus1s(1);
            disturb_levelone(k) = pre_proj(k,1); 
            
            
            % ++ Accumulate the primal variance in system to determine quality ++ %
            var_lambda = 0.999;
            e_var = abs(disturb_levelone(k)- mean_levelone(k));
            var_vec(1) = var_lambda*var_vec(1) + (1-var_lambda)*e_var ; 
            
            % Store path quality (simulation only) 
            var_vec_hist(m,k) = var_vec(1) ; 
            
            % ++ Convergence detector ++ %
            e_sf = abs(proj_vec_dminus1f(int_index)-proj_vec_dminus1s(int_index)) ;
            if (e_sf < 0.1*proj_vec_dminus1s(int_index)) && (proj_vec_dminus1s(int_index)<0.1) && vpcktcnt>300
                conv_count = conv_count+1 ;
            else
                conv_count = 0 ;
            end
            if conv_count>20
                conv_flag = 1;
            end
            conv_vec_hist(m,k) = conv_flag ;

        end
        

            % The adaptive computing of interference directions
            v = ch_sel1;    
            for j=1:13
                if vpcktcnt>500 && pre_proj(k,1)<(0.8*(proj_vec_dminus1s(1)))
                    rej_flag(k)=1;
                    break; 
                end
                % Amari notation
                y = W(:,j)'*v;
                p_wts(j)=(beta(j).*p_wts(j))+((1)*(abs(y)^2));
                e_p = v - W(:,j)*y;
                if j==6 % debug only condition
                    e1(k) = abs(e_p'*e_p); %pre_proj(k)=(abs(y)^2) ; %p_wts(j);
                    W1(:,k) = W(:,6);
                end
                W(:,j) = W(:,j) + (e_p*(y'/p_wts(j)));
                v = v - W(:,j)*y ; 
            end

            e_res(k) = e_p'*e_p; 
            p_comp(k,:) = p_wts;   

    end

start_off=10;
% figure(645) ; plot(disturb_metr(start_off:end)); hold on ; plot(mean_level(start_off:end),'r'); hold off
figure(646) ; plot(disturb_levelone) ; hold on ; plot(mean_levelone,'r') ; hold off
motion_out_ppath(m,:)  = disturb_metr(start_off:end)./mean_level(start_off:end);
path_metric(m) = var_vec(1);
%close all 

end

% Generate final motion output : Combine top three paths (in practice wud be inline)
sel_inds = 1:length(motion_out_ppath);
inv_var_hist = (1./var_vec_hist(:,sel_inds)) ;
[vals inds] = sort(inv_var_hist,'descend');

for k=1:length(motion_out_ppath)
    [vals inds] = sort(inv_var_hist(:,k),'descend');
    eff_wts = inv_var_hist(inds(1:2),k);
    motion_out_normed = tanh(motion_out_ppath./100).*100 ;
    eff_motion = motion_out_normed(inds(1:2),k);
    wtd_motion(k) = sum(eff_motion)/2 ;
end

% figure ;plot(wtd_motion);
% title('Weighted motion output')
% 
% 
% figure ; plot(median(motion_out_ppath))
% title('Median of aggregated motion paths')
% 
% sorted_motion_out = sort(motion_out_ppath);
% figure ; plot(sorted_motion_out(end-1,:))
% title('Max minus one of aggregated motion paths')
% 
% figure ; plot(tanh(sorted_motion_out(end-1,:)./100).*100)
% title('Max minus one with HyperTan')
% 
% sorted_motion_out = sort(motion_out_ppath);
% figure ; plot(sorted_motion_out(end-2,:))
% title('Max minus two of aggregated motion paths')
% 
% 
% nl_normalized = tanh(motion_out_ppath./100).*100 ;
% max_ratio_motion = nl_normalized.*repmat(transpose((1./(path_metric))),1,length(motion_out_ppath));
% max_ratio_motion = max_ratio_motion./(sum(1./path_metric));
% figure ; plot(sum(max_ratio_motion))
% title('Maximal ratio combining of motion paths - theoretical only')
% 
% % convergence plot
% figure ; plot(conv_vec_hist') ; axis([0 length(conv_vec_hist) -1 +2])
% title('Convergence of different paths')

sorted_motion_out = sort(motion_out_ppath);
final_motion_out = sorted_motion_out(end-1,:)>3 ;

% figure ; plot(st_montime(10:end),final_motion_out) ; hold on
% h_raw_allq = h_raw_allq_cell{test_path};
% plot(st_montime(1:end),abs(h_raw_allq(:,29))','r');
% axis([0 st_montime(end) 0 3])


end