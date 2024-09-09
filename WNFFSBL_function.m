%% Consider an M elements vertical array, and a shallow water environment in ocean.
%% rr is the  range grid points, dp is the depth grid points, the number of the total grid points is N.
%% Using  the true source  position and source parameters to generate the recieved acoustic data Y(:,:,:) in frequency domain by kraken, Y(:,:,:)  is a (M,L,F) size matrix
%% Generate dictionary Green_all(:,:,:,), which is a (M,N,F) size  matrix.

% then use the following codes to get input:
% ii=0;
% for i=1:length(rr)
%     for j=1:length(dp)
%         ii=ii+1;
%         g_all(:,ii,fi)=[Green_all(j,:,i,fi)];
%     end
% end
% glsum=sqrt(sum(abs(g_all(:,:,fi)).^2,1));
% for i=1:size(glsum,2) %归一化，使得范数为1
% g_all(:,i,fi)=g_all(:,i,fi)/glsum(i);
% end
% Aa=g_all;
% Target_value= ((sum(abs(g_all(:,:,1)'*Yf(:,:,1)).^2,2)))./(sum(abs(g_all(:,:,1)).^2,1).');
function [SPENFFSP,SPENFFSBL,SP_target,count,rroff_grid,dpoff_grid]=MNFFSBL(Y,Aa,rr,dp,Target_value,K)
%% Input { Y:data, Aa(M*N*F matrix): Overcomplete dictionary, rr:range gridpoints, dp:depth gridpoints, Target_value}


[M,L,F]=size(Y);
N1=length(rr);
N2=length(dp);
N=N1*N2;
a_=10^(-6);
b_=10^(-4);
b=10^(-6);
es=10^(-12);
    rroff_grid=rr;
    dpoff_grid=dp;
source_power=ones(N,1);
[max_value,SP_target]=max(Target_value);
source_power(SP_target) =  (max_value);
SP_E=source_power(SP_target);
A_=[];B_=[];C_=[];Z=zeros(N,F);Q=zeros(N,L,F);W=zeros(1,F);
q=zeros(N,L,F);
for ff=1:F
Aa_(:,:,ff)=Aa(:,SP_target,ff);
gamma=diag(SP_E);
sigma_y(:,:,ff)=eye(M)+Aa_(:,:,ff)*gamma*Aa_(:,:,ff)';
sigma_y_inv(:,:,ff)=inv(sigma_y(:,:,ff));
    Mu_x(:,:,ff)=gamma*Aa_(:,:,ff)'*sigma_y_inv(:,:,ff)*Y(:,:,ff);
    sigma_x(:,:,ff)=inv(Aa_(:,:,ff)'*Aa_(:,:,ff)+gamma^(-1));
%% SS,QQ,method2:
Aa_sigma_y_inv=Aa(:,:,ff)'*sigma_y_inv(:,:,ff);
for i=1:N
    Z(i,ff)=  Aa_sigma_y_inv(i,:)*Aa(:,i,ff);
end
Q(:,:,ff)=(Aa(:,:,ff)'*sigma_y_inv(:,:,ff)*Y(:,:,ff));
W(ff)=real(trace(Y(:,:,ff)'*sigma_y_inv(:,:,ff)*Y(:,:,ff))+b);
end
%% Keep track of the positions selected during the iterations
selected = SP_target;
deleted = [];
max_it =300;
co=0;add=0;res=0;

for count = 1:max_it
    source_power_old=source_power;
    co=co+1;
    for ff=1:F
    %% 相当于功率为0时初始化z,q,g：
    z(:,ff)=real(Z(:,ff));
    q(:,:,ff)=Q(:,:,ff);
    w(:,ff)=ones(N,1)*W(ff);
    %% 只对利用的原子进行赋值
    z(SP_target,ff)=real(Z(SP_target,ff)./(1-SP_E.*Z(SP_target,ff)));
    for l=1:L
        q(SP_target,l,ff)=Q(SP_target,l,ff)./(1-SP_E.*Z(SP_target,ff));
    end
    q_2(:,ff)=sum(abs(q(:,:,ff)).^2,2);
    sum_Q2=sum(abs(Q(:,:,ff)).^2,2);
    w(SP_target,ff)=W(ff)+SP_E.*sum_Q2(SP_target)./(1-SP_E.*real(Z(SP_target,ff)));
    
    A_(:,ff)=b_/F*z(:,ff).*q_2(:,ff)./w(:,ff)-b_/F*abs(z(:,ff)).^2;
    B_(:,ff)=-L*abs(z(:,ff)).^2+L*z(:,ff).*q_2(:,ff)./w(:,ff)-2*b_/F*z(:,ff)+b_/F*q_2(:,ff)./w(:,ff);
    C_(:,ff)=(M*L+a_)*q_2(:,ff)./w(:,ff)-b_/F-L*z(:,ff);
    end
 %%   
      Asum=sum(A_,2);
      Bsum=sum(B_,2);
      Csum=sum(C_,2);
    delta=  Bsum.^2-4*Asum.*Csum;
    delta_d0= find(delta<0);
    source_power=(- Bsum-sqrt(Bsum.^2-4*Asum.*Csum ))./(2*Asum );
    ml = -inf*ones(N,1);  %initialize ml function
    ml_ = -inf*ones(N,F);  %initialize ml function
    %       ML=ML_;     % ML_ : new  , ML : old
    iglower0 = find(source_power<0);
    source_power(iglower0)=es;
    C_A=sum(C_,2)./ sum(A_,2);
    CA_positive= find(C_A>0);
    ig0 = find(source_power>es);
    
    %% indices for re-estimation, find the index which satisfy the re-estimation condition:
    [id_re,~,which] = intersect(ig0,SP_target);  % return the same entry in two arrays
    
    if ~isempty(id_re)
        SP_T = source_power(id_re);
        for ff=1:F
        ml_(id_re,ff) =-(M*L+a_)*log(1-sum((abs(q(id_re,:,ff))).^2,2)./w(id_re,ff)./ (SP_T.^(-1)+z(id_re,ff) ) ) -L*log(1+SP_T.*z(id_re,ff))-b_/F*SP_T-...
            (-(M*L+a_)*log(1-sum((abs(q(id_re,:,ff))).^2,2)./w(id_re,ff)./ (SP_E(which).^(-1)+z(id_re,ff) ) ) -L*log(1+SP_E(which).*z(id_re,ff))-b_/F*SP_E(which) ) ;
        end
        ml=sum(ml_,2);
    end
    
    %% indices for adding,find the index which satisfy the adding condition:
    id_ad = setdiff(ig0,id_re);  % 返回前一个矩阵里有而后一个矩阵里没有的元素
    if ~isempty(id_ad)
        SP_T = source_power(id_ad);
        for ff=1:F
        ml_(id_ad,ff) =-(M*L+a_)*log( 1-sum((abs(q(id_ad,:,ff))).^2,2)./w(id_ad,ff)./ (SP_T.^(-1)+z(id_ad,ff) ) ) -L*log(1+SP_T.*z(id_ad,ff))-b_/F*SP_T;
        end
        ml=sum(ml_,2);
        which = intersect(deleted,id_ad);
        ml(which) = -inf;
    end
    is0 = setdiff([1:N],ig0); % 找出delta<0的基的下标
    %% indices for deleting
    [id_del,foo,which] = intersect(is0,SP_target);
    if ~isempty(id_del)
        
        if length(SP_target) == 1
            ml(id_del) = -inf;
        else
             for ff=1:F
            ml_(id_del,ff) =((M*L+a_)*log( 1-sum((abs(q(id_del,:,ff))).^2,2)./w(id_del,ff) )./(SP_E(which).^(-1)+z(id_del,ff)) ) +L*log(1+SP_E(which).*z(id_del,ff))+b_/F*SP_E(which);
             end
             ml=sum(ml_,2);
        end
        
    end
    
    [ML(count),idx] = max(real((ml)));
    [ ml_sort,ml_index]= sort(real(ml),'descend');
    
    %   SP_nozero=ml_index(1:length(SP_nozero));
    %% check convergence
    tol=10^(-7);
    if count >=2
        if  norm(ML(count)-ML(count-1)) < norm(ML(count)-ML(1))*tol
%        if norm((source_power-source_power_old),2)/norm(source_power_old,2)< 0.001
            break;
        end
    end
    %% update alphas
    % Choose the basis which results in the largest increase in the
    % likelihood
    
    which = find(SP_target==idx);%判断已有的基里对应最大化ML的下标

    

         if source_power(idx) > es
            
        if ~isempty(which) % reestimate a basis
            res=res+1;
            %% 附录快速方法
             for ff=1:F
            SP_T=source_power(idx);
            Sigii = real(sigma_x(which,which,ff));
            Sigi = ( sigma_x(:,which,ff) );
            mui = Mu_x(which,:,ff);
            delta_sp=SP_T^(-1)-SP_E(which)^(-1);
            ki = (delta_sp)/(1+Sigii*(delta_sp) );
            for i=1:length(mui)
                Mu_x(:,i,ff) = Mu_x(:,i,ff)-ki*mui(i)*Sigi;
            end
            sigma_x(:,:,ff) = sigma_x(:,:,ff)-ki*Sigi*Sigi';
            comm = Aa(:,:,ff)'*(Aa_(:,:,ff)*Sigi);
            Z(:,ff) = Z(:,ff) + ki*abs(comm).^2;
            for i=1:length(mui)
                Q(:,i,ff) = Q(:,i,ff) + ki*mui(i)*comm; %
            end

            W(ff)=W(ff)+ki*sum(abs(Y(:,:,ff)'*(Aa_(:,:,ff)*Sigi)).^2);
            
             end
             SP_E(which) = SP_T;
        else
            sigma_x_=sigma_x;
          Mu__x=Mu_x;
          A_a_=Aa_;
            %% add a basis
            SP_T=source_power(idx);
            add=add+1;
         for ff=1:F
            
            Aaii = Aa(:,idx,ff); Sigii = 1/(SP_T^(-1)+real(Z(idx,ff)));
            mui = Sigii*Q(idx,:,ff);
            comm1 = (sigma_x_(:,:,ff))*(A_a_(:,:,ff)'*Aaii);
            ei = Aaii-A_a_(:,:,ff)*comm1;
            off = -Sigii*comm1;
            sigma_x__ = [sigma_x_(:,:,ff)+Sigii*comm1*comm1', off; off', Sigii];
            if ff==1
            sigma_x=zeros(size(sigma_x__,1),size(sigma_x__,2),F);
            end
            sigma_x(:,:,ff)=sigma_x__;
            Mu_x_=Mu__x(:,:,ff);
            for i=1:length(mui)
                Mu_x_(:,i)=   Mu__x(:,i,ff)-mui(i) *comm1;
            end
                if ff==1
            Mu_x=zeros(size(Mu_x_,1)+1,size(Mu_x_,2),F);
                end
            Mu_x(:,:,ff) = [Mu_x_; mui];
            comm2 = Aa(:,:,ff)'*ei;
            Z(:,ff) = Z(:,ff)- Sigii*abs(comm2).^2;
            for i=1:length(mui)
                Q(:,i,ff)  = Q(:,i,ff)  - mui(i)*comm2;
            end
            W(ff)=W(ff)-Sigii*sum(abs(Y(:,:,ff)'*ei).^2);

            Aa__ = [A_a_(:,:,ff),Aaii];
                if ff==1
            Aa_=zeros(size(Aa__,1),size(Aa__,2),F );
                end
            Aa_(:,:,ff)=Aa__;
         end
             SP_target = [SP_target;idx];
            SP_E = [SP_E;SP_T];
        end
    else
        if ~isempty(which) && length(SP_target) > 1 % delete a basis
             for ff=1:F
            %% 
            Sigii = real(sigma_x(which,which,ff)); mui = Mu_x(which,:,ff); Sigi = (sigma_x(:,which,ff));
            sigma_x(:,:,ff) = sigma_x(:,:,ff)-Sigi*Sigi'/Sigii; 

            for i=1:length(mui)
                Mu_x(:,i,ff)=   Mu_x(:,i,ff)-mui(i)/Sigii*Sigi;
            end
            
            comm = Aa(:,:,ff)'*(Aa_(:,:,ff)*Sigi);
            Z(:,ff) = Z(:,ff) + abs(comm).^2/Sigii;  % 一定要加上abs模值
            for i=1:length(mui)
                Q(:,i,ff)  = Q(:,i,ff)  + mui(i)/Sigii*comm;
            end
            W(ff)=W(ff)+Sigii^(-1)*sum(abs(Y(:,:,ff)'*(Aa_(:,:,ff)*Sigi)).^2);

            
             end
             deleted = [deleted idx];
               SP_E(which) = [];
            source_power(SP_target(which))=0;
             SP_target(which) = [];
             sigma_x(:,which,:) = []; Mu_x(which,:,:) = [];
            sigma_x(which,:,:) = [];
             Aa_(:,which,:) = [];
        elseif ~isempty(which) && length(SP_target) == 1
            % Something is wrong, trying to delete the only coefficient
            % that has been added.
            break;
      
        end
            
        end
        

    selected =union(selected,idx);
  
         end

%     SPENFFSP=zeros(N1*N2,1);
%     SPENFFSP_=zeros(N1*N2,1);
        SPENFFSP=ones(N1*N2,1)*es;
    SPENFFSP_=ones(N1*N2,1)*es;
    [SP_target,spid]=sort(SP_target);
    SP_E=SP_E(spid);
    SPENFFSP(SP_target)=real(SP_E)/max(real(SP_E) );
    SPENFFSP_(SP_target)=real(SP_E);
   source_power= SPENFFSP_;
    source_id=[];
    for i=1:N1
        SPENFFSBL(:,i)=SPENFFSP((i-1)*N2+1:(i-1)*N2+N2);
    end
    
%%  refined DOA original:
     SP_id=imregionalmax(SPENFFSBL);
    [depth_id,range_id]= find(SP_id~=0);
%     K_e= min(length(depth_id),M-1); % if K is not known
   K_e= min(length(depth_id),K); % if K is  known
    Power_select=[];
    for i=1:length(range_id)
    Power_select(i)=SPENFFSBL(depth_id(i),range_id(i));
    end
    [~,peaks]=sort(Power_select,'descend');
    Position_id=[range_id( peaks(1:K_e)),depth_id( peaks(1:K_e))];  % Position_es第一列是距离，第二列是深度。
    [~,ascend_rid]=sort(Position_id(:,1),'ascend');
    Position_id=[Position_id(ascend_rid,:)]; % 按照距离升序排列行，以便于后边RMSE的计算。
    


for i=1:size(Position_id,1)
%     %%
    ri=0;step=0.2;
   
    Aa_k=[];qs=[];
     if (Position_id(i,1))>1&&((Position_id(i,1)-1)*N2+Position_id(i,2) +N2)<=N %避免索引出错
for ff=1:F
   sigma_y(:,:,ff)=eye(M)+Aa_(:,:,ff)*diag(source_power(SP_target))*Aa_(:,:,ff)';
        Aa_k= [Aa(:, (Position_id(i,1)-1)*N2+Position_id(i,2) ,ff)];         % exclude the angle index
        gamma_k= [source_power((Position_id(i,1)-1)*N2+Position_id(i,2) )];
        gamma_k=diag(gamma_k);
        sigma_yk=sigma_y(:,:,ff)-Aa_k* gamma_k*Aa_k';    % 
        sigma_yk_inv(:,:,ff)=inv( sigma_yk);
        
         gk(ff)=real(trace(Y(:,:,ff)'*sigma_yk_inv(:,:,ff)*Y(:,:,ff))+b);
end
     end
     
if (Position_id(i,1))>1&&((Position_id(i,1)-1)*N2+Position_id(i,2) +N2)<=N %避免索引出错
 sp_sumQ1=source_power((Position_id(i,1)-1)*N2+Position_id(i,2) +size(SPENFFSBL,1) ) +source_power((Position_id(i,1)-1)*N2+Position_id(i,2) -1 ) +source_power((Position_id(i,1)-1)*N2+Position_id(i,2) +size(SPENFFSBL,1)-1 ) ;
 sp_sumQ2=source_power((Position_id(i,1)-1)*N2+Position_id(i,2) -size(SPENFFSBL,1) ) +source_power((Position_id(i,1)-1)*N2+Position_id(i,2) -1 ) +source_power((Position_id(i,1)-1)*N2+Position_id(i,2) -size(SPENFFSBL,1)-1 ) ;
 sp_sumQ3=source_power((Position_id(i,1)-1)*N2+Position_id(i,2) -size(SPENFFSBL,1) ) +source_power((Position_id(i,1)-1)*N2+Position_id(i,2) +1 ) +source_power((Position_id(i,1)-1)*N2+Position_id(i,2) -size(SPENFFSBL,1)+1 ) ;
 sp_sumQ4=source_power((Position_id(i,1)-1)*N2+Position_id(i,2) +size(SPENFFSBL,1) ) +source_power((Position_id(i,1)-1)*N2+Position_id(i,2) +1 ) +source_power((Position_id(i,1)-1)*N2+Position_id(i,2) +size(SPENFFSBL,1)+1 ) ;
spr_array=[sp_sumQ1,sp_sumQ2,sp_sumQ3,sp_sumQ4];
   [spQid]= max(spr_array);
    
 if  spQid==sp_sumQ1
     spre=[];L_off=[];
     
         for l1=0.1:step:0.9
               ri=ri+1;di=0;
        for l2=0.1:step:0.9
           di=di+1;
        aa=[];qk=[];zk=[];
        for ff=1:F
         aa(:,ff)=Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) ,ff)+l1*( Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) +size(SPENFFSBL,1),ff)-Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) ,ff) )+...
             l2*( Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) +1,ff)-Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) ,ff));
         qs(:,:,ff)=aa(:,ff)'*sigma_yk_inv(:,:,ff)*Y(:,:,ff);
        qk_sum(:,ff)= sum(abs(qs(:,:,ff).^2),2); zk(ff)=real(aa(:,ff)'*sigma_yk_inv(:,:,ff)*aa(:,ff));
       A_k(:,ff)=b_/F*zk(ff).*qk_sum(:,ff)./gk(ff)-b_/F*abs(zk(ff)).^2;
       B_k(:,ff)=-L*abs(zk(ff)).^2+L*zk(ff).*qk_sum(:,ff)./gk(ff)-2*b_/F*zk(ff)+b_/F*qk_sum(:,ff)./gk(:,ff);
       C_k(:,ff)=(M*L+a_)*qk_sum(:,ff)./gk(ff)-b_/F-L*zk(ff);
        end
         spre= max( (- sum(B_k,2)-sqrt(sum(B_k,2).^2-4*sum(A_k,2).*sum(C_k,2)))./(2*sum(A_k,2)), (- sum(B_k,2)+sqrt(sum(B_k,2).^2-4*sum(A_k,2).*sum(C_k,2)))./(2*sum(A_k,2))) ; %
         for ff=1:F
         L_off(ff) =-(M*L+a_)*log( 1-qk_sum(:,ff)./gk(ff)./ (spre.^(-1)+zk(ff) ) ) -L*log(1+spre.*zk(ff))-b_/F*spre;
         end
        L_offgrid(ri,di)=real(sum(L_off));   
    end
    end
     [r_max(i),d_max(i)]=find((L_offgrid)==max(max(L_offgrid)));
    ll1=0.1:step:0.9;
    ll2=0.1:step:0.9;
    r_offgrid(i)=rr(Position_id(i,1))+0.05*ll1(r_max(i));
    d_offgrid(i)=dp(Position_id(i,2))-5*ll2(d_max(i));
    rroff_grid(Position_id(i,1))=r_offgrid(i);
    dpoff_grid(Position_id(i,2))= d_offgrid(i);
    
 elseif spQid==sp_sumQ2
         for l1=0.1:step:0.9
             ri=ri+1;di=0;
        for l2=0.1:step:0.9
        di=di+1;
        aa=[];qk=[];zk=[];
        for ff=1:F
         aa(:,ff)=Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) ,ff)+l1*(Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) -size(SPENFFSBL,1),ff)  -Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) ,ff))+...
             l2*(  Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) -1,ff)-Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) ,ff));
         qs(:,:,ff)=aa(:,ff)'*sigma_yk_inv(:,:,ff)*Y(:,:,ff);
        qk_sum(:,ff)= sum(abs(qs(:,:,ff).^2),2); zk(ff)=real(aa(:,ff)'*sigma_yk_inv(:,:,ff)*aa(:,ff));
       A_k(:,ff)=b_/F*zk(ff).*qk_sum(:,ff)./gk(ff)-b_/F*abs(zk(ff)).^2;
       B_k(:,ff)=-L*abs(zk(ff)).^2+L*zk(ff).*qk_sum(:,ff)./gk(ff)-2*b_/F*zk(ff)+b_/F*qk_sum(:,ff)./gk(:,ff);
       C_k(:,ff)=(M*L+a_)*qk_sum(:,ff)./gk(ff)-b_/F-L*zk(ff);
        end
spre= max( (- sum(B_k,2)-sqrt(sum(B_k,2).^2-4*sum(A_k,2).*sum(C_k,2)))./(2*sum(A_k,2)), (- sum(B_k,2)+sqrt(sum(B_k,2).^2-4*sum(A_k,2).*sum(C_k,2)))./(2*sum(A_k,2))) ; %         for ff=1:F
         for ff=1:F
L_off(ff) =-(M*L+a_)*log( 1-qk_sum(:,ff)./gk(ff)./ (spre.^(-1)+zk(ff) ) ) -L*log(1+spre.*zk(ff))-b_/F*spre;
         end
        L_offgrid(ri,di)= real(sum(L_off));     
    end
         end
         [r_max(i),d_max(i)]=find((L_offgrid)==max(max(L_offgrid)));
    ll1=0.1:step:0.9;
    ll2=0.1:step:0.9;
    r_offgrid(i)=rr(Position_id(i,1))-0.05*ll1(r_max(i));
    d_offgrid(i)=dp(Position_id(i,2))-5*ll2(d_max(i));
    rroff_grid(Position_id(i,1))=r_offgrid(i);
    dpoff_grid(Position_id(i,2))= d_offgrid(i);
     elseif spQid==sp_sumQ3
         for l1=0.1:step:0.9
             ri=ri+1;di=0;
        for l2=0.1:step:0.9
        di=di+1;
        aa=[];qk=[];zk=[];
        for ff=1:F
         aa(:,ff)=Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) ,ff)+l1*( Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) -size(SPENFFSBL,1),ff)-Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) ,ff) )+...
             l2*( Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) +1,ff)-Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) ,ff));
         qs(:,:,ff)=aa(:,ff)'*sigma_yk_inv(:,:,ff)*Y(:,:,ff);
        qk_sum(:,ff)= sum(abs(qs(:,:,ff).^2),2); zk(ff)=real(aa(:,ff)'*sigma_yk_inv(:,:,ff)*aa(:,ff));
       A_k(:,ff)=b_/F*zk(ff).*qk_sum(:,ff)./gk(ff)-b_/F*abs(zk(ff)).^2;
       B_k(:,ff)=-L*abs(zk(ff)).^2+L*zk(ff).*qk_sum(:,ff)./gk(ff)-2*b_/F*zk(ff)+b_/F*qk_sum(:,ff)./gk(:,ff);
       C_k(:,ff)=(M*L+a_)*qk_sum(:,ff)./gk(ff)-b_/F-L*zk(ff);
        end
spre= max( (- sum(B_k,2)-sqrt(sum(B_k,2).^2-4*sum(A_k,2).*sum(C_k,2)))./(2*sum(A_k,2)), (- sum(B_k,2)+sqrt(sum(B_k,2).^2-4*sum(A_k,2).*sum(C_k,2)))./(2*sum(A_k,2))) ; %         for ff=1:F
        for ff=1:F
L_off(ff) =-(M*L+a_)*log( 1-qk_sum(:,ff)./gk(ff)./ (spre.^(-1)+zk(ff) ) ) -L*log(1+spre.*zk(ff))-b_/F*spre;
         end
        L_offgrid(ri,di)= real(sum(L_off));    
    end
         end
              [r_max(i),d_max(i)]=find((L_offgrid)==max(max(L_offgrid)));

    ll1=0.1:step:0.9;
    ll2=0.1:step:0.9;
    r_offgrid(i)=rr(Position_id(i,1))-0.05*ll1(r_max(i));
    d_offgrid(i)=dp(Position_id(i,2))+5*ll2(d_max(i));
    rroff_grid(Position_id(i,1))=r_offgrid(i);
    dpoff_grid(Position_id(i,2))= d_offgrid(i);
 else
         for l1=0:step:0.9
               ri=ri+1;di=0;
        for l2=0:step:0.9
      di=di+1;
        aa=[];qk=[];zk=[];
        for ff=1:F
         aa(:,ff)=Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) ,ff)+l1*( Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) +size(SPENFFSBL,1),ff)-Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) ,ff) )+...
             l2*( Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) +1,ff)-Aa(:,(Position_id(i,1)-1)*N2+Position_id(i,2) ,ff));
         qs(:,:,ff)=aa(:,ff)'*sigma_yk_inv(:,:,ff)*Y(:,:,ff);
        qk_sum(:,ff)= sum(abs(qs(:,:,ff).^2),2); zk(ff)=real(aa(:,ff)'*sigma_yk_inv(:,:,ff)*aa(:,ff));
       A_k(:,ff)=b_/F*zk(ff).*qk_sum(:,ff)./gk(ff)-b_/F*abs(zk(ff)).^2;
       B_k(:,ff)=-L*abs(zk(ff)).^2+L*zk(ff).*qk_sum(:,ff)./gk(ff)-2*b_/F*zk(ff)+b_/F*qk_sum(:,ff)./gk(:,ff);
       C_k(:,ff)=(M*L+a_)*qk_sum(:,ff)./gk(ff)-b_/F-L*zk(ff);
        end
spre= max( (- sum(B_k,2)-sqrt(sum(B_k,2).^2-4*sum(A_k,2).*sum(C_k,2)))./(2*sum(A_k,2)), (- sum(B_k,2)+sqrt(sum(B_k,2).^2-4*sum(A_k,2).*sum(C_k,2)))./(2*sum(A_k,2))) ; %         for ff=1:F
        for ff=1:F
L_off(ff) =-(M*L+a_)*log( 1-qk_sum(:,ff)./gk(ff)./ (spre.^(-1)+zk(ff) ) ) -L*log(1+spre.*zk(ff))-b_/F*spre;
         end
        L_offgrid(ri,di)= real(sum(L_off,2));   
    end
    end
          [r_max(i),d_max(i)]=find((L_offgrid)==max(max(L_offgrid)));

    ll1=0.1:step:0.9;
    ll2=0.1:step:0.9;
    r_offgrid(i)=rr(Position_id(i,1))+0.05*ll1(r_max(i));
    d_offgrid(i)=dp(Position_id(i,2))+5*ll2(d_max(i));
    rroff_grid(Position_id(i,1))=r_offgrid(i);
    dpoff_grid(Position_id(i,2))= d_offgrid(i);
      end

end
end
    
end