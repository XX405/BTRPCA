% function [X]=truncatedL1L2(M, L, Z,par)
% t_now = 20;
% %% parameter
% % time_DCA = tic;
% 
% [m,n,h] = size(M);
% 
% 
% 
% [m,n,h]=size(L);
% L=[m*h,n];
% [m,n,h]=size(Z);
% Z=[m*h,n];
% lambda = 0.003;
% 
% if ~isfield(par,'iter_out_min') 
%     par.iter_out_min=5;
% end
% iter_out_min=par.iter_out_min;
% if ~isfield(par,'beta');par.beta=3e-7;end;beta1=par.beta;
% if ~isfield(par,'eps_min');par.eps_min=1e-10;end
% if ~isfield(par,'tol');par.tol=[1e-5,1e-5];end;tol_out=par.tol(1);tol_in=par.tol(2);
% if ~isfield(par,'c');par.c=0;end;c=par.c;if(par.c<0);error('wrong input: par.c');end
% if ~isfield(par,'show_fval');par.show_fval=1;end;show_fval=par.show_fval;
% if ~isfield(par,'detail');par.detail=1;end
% if ~isfield(par,'test_time');par.test_time=0;end
% if ~isfield(par,'maxiter');par.maxiter=[100,1000];end
% maxiter_out=100;maxiter_in=1000;
% iter_in_min=5; tol_out=1e-5;iter_out_min=5;
% %% initialization
% % cr_in=zeros(1,maxiter_out);
% % cr_x=zeros(1,maxiter_out);
% % cr_f=zeros(1,maxiter_out);
% % iter_in=zeros(1,maxiter_out);
% % con_Y=zeros(1,maxiter_out);
% % err=zeros(1,maxiter_out);
% % fval=zeros(1,maxiter_out);
% % t_all=zeros(1,maxiter_out);
% % Lre=zeros(1,maxiter_out);
% 
% X = zeros( m*n,h);
% S = zeros( m*n,h);
% % for t=1:h
% %     Xt=X(:,:,t);
% %     St=S(:,:,t);
% % end
% %% DCA begin
% 
% for i=1:maxiter_out
% 
%     %% compute P
%      % for t=1:h
%      %    S_h=S_temp(:,t);
%      % end
%     [U_mc,S_temp,V_mc]=svd(X,'econ');
% 
%     if all(size(S_temp)-1)
% 
%          St=S_temp(:,h);
% 
%         S_temp=diag(St);
%      end
%     S_temp=S_temp(:);
%     eps_min=1e-10;
%     T=S_temp>eps_min;
%     r_temp=sum(T);%rank of X
% 
%     if r_temp<=t_now
%         P_temp=ones(r_temp,1);
%     else
%         P_temp=[ones(t_now,1);...
%             alpha*S_temp(t_now+1:r_temp)/norm(S_temp(t_now+1:r_temp),'fro')];
%     end
%     P=U_mc(:,1:r_temp)*diag(P_temp)*V_mc(:,1:r_temp)';
%     c=1;
%     P=P+c/lambda*X;
% 
%     %% ADMM
%     for j=1:maxiter_in
%         mu = 1e-4; 
% 
%          beta=3e-7;
%      lambda = 0.003;  
%          Y = (mu*(M-L+Z/mu)+lambda*P+beta*(X+S))/(mu+c+beta);
%         % Y = (par.mu*(D-E+L/par.mu)+lambda*P+beta1*(X+S))/(par.mu+c+beta1);
%         %% compute X
% 
%         X=SVT(Y-S,lambda/beta,'L1');
% 
%         %% compute S
%         S=S+X-Y;
%         %% stop inner loop
%     %     tol_in=1e-5;iter_out_min=5;
%     %     cr_in1=norm(X-Xpin,'fro')/max(norm(Xpin,'fro'),eps_min);
%     %     if cr_in1<tol_in
%     %         stop_flag= (i<iter_out_min+1 || j>=iter_in_min);
%     %         if i==1 && j<=50
%     %             stop_flag=0;
%     %         end
%     %         if stop_flag
%     %             break;
%     %         end
%     end
% 
%  end
%     % stop outer loop
%     % cr_in(i)=cr_in1;
%     % cr_x1=norm(X-Xpout,'fro')/max(norm(Xpout,'fro'),eps_min);
%     % cr_x(i)=cr_x1;
%     % stop_flag=cr_x1<tol_out && i>=iter_out_min+1;
%     show_fval=1;
%     if show_fval
%         fval(i)=fun_fval(M,X,L,Z,lambda,t_now);
%         if i>=2
%             cr_fval=abs(fval(i)-fval(i-1))/max(fval(i-1),eps_min);
%         else
%             cr_fval=1;
%         end
%         cr_f(i)=cr_fval;
%         % stop_flag=(cr_fval<tol_out || cr_x1<tol_out) && i>=iter_out_min+1;
%     end
% %     iter_in(i)=j;
% %     con_Y(i)=norm(X-Y,'fro');
% %     err(i)=norm(D-X-L,'fro')/norm(D,'fro');
% %     t_all(i)=t_now;
% %     detail=1;
% %     if detail
% %         if i==1
% %             result_title='OutIt| InnIt|    ReErr|   OutCrX|   OutCrF|   InnCrX| Trun#';
% %             disp(result_title)
% %         end
% %         result_sprintf=sprintf('%5d| %5d| %8.2e| %8.2e| %8.2e| %8.2e| %5d',...
% %             i,iter_in(i),err(i),cr_x(i),cr_f(i),cr_in(i),t_all(i));
% %         disp(result_sprintf)
% %     end
% %     test_time=0;
% %     if test_time%not record time
% %         if stop_flag
% %             break;
% %         end
% %     else%record time
% %         if norm(X-Xpout,'fro')/norm(Xpout,'fro')<test_time
% %             %         if norm(X-par.x_t,'fro')/norm(par.x_t,'fro')<par.test_thre
% %             break;
% %         end
%     end
% 
% 
% 
% %% SVT
% function X=SVT(Z,tau,type)
% %for computing SVT of a matrix
% if nargin<3
%     type='L1';
% end
% 
% 
% [m,n,h]=size(Z);
% for ii=h
% 
% Z=Z(:,:,ii);
% end
% if m <= n
%     AAT=Z*Z';
%     [S,V]=eig(AAT);
%     V=max(0,diag(V));
%     V=sqrt(V);
%     tol=n*eps(max(V));
%     mid=thresholdings(V,tau(end:-1:1),type);
%     ind=mid>tol;
%     if any(ind)
%         mid=mid(ind)./V(ind);
%         X=S(:,ind)*diag(mid)*S(:,ind)'*Z;
%     else
%         X=zeros(m,n);
%     end
% else
%     X=SVT(Z',tau,type);
%     X=X';
% end
% end
% 
% function x=thresholdings(y,tau,type)
% %solve: min_x: 1/2(x-y)^2+tau*f(x)
% %type:
% % 'L1': f(x)=|x|
% % 'L1/2': f(x)=sqrt(x)
% if ~exist('type','var')
%     type='L1';
% end
% switch type
%     case 'L1'
%         x=sign(y).*max(0,abs(y)-tau);
%     case 'L1/2'
%         ind=abs(y)>3/2*tau.^(2/3);
%         x=zeros(size(y));
%         y=y(ind);
%         x(ind)=2/3*y.*(1+cos(2/3*(pi-acos(tau/4.*(abs(y)/3).^(-3/2)))));
% end
% end
% 
% function v=fun_fval(D, X, E, L, lambda,  t_now)
% alpha=1e-4;mu=1e-4;
% X_temp=svd(X,'econ');
% v1 = lambda*(TL1L2norm(X_temp,t_now,alpha));
% v2 = mu/2*norm(D-X-E+L/mu,'fro')^2;
% v = v1+v2;
% end
% 
% function val=TL1L2norm(Z,t,alpha)
% %compute truncated l_{1-2} metric
% Z=sort(abs(Z),'descend');
% val=sum(Z(t+1:end))-alpha*sqrt(sum(Z(t+1:end).^2));
% end
function [X]=My_truncatedL1L2(M, L, Z, par)
t_now = 20;
%% parameter
% time_DCA = tic;
[m,n,h] = size(M);
lambda = 0.003;
if ~isfield(par,'iter_out_min');par.iter_out_min=5;end
iter_out_min=par.iter_out_min;
if ~isfield(par,'beta');par.beta=3e-7;end;beta1=par.beta;
if ~isfield(par,'eps_min');par.eps_min=1e-10;end
if ~isfield(par,'tol');par.tol=[1e-5,1e-5];end;tol_out=par.tol(1);tol_in=par.tol(2);
if ~isfield(par,'c');par.c=0;end;c=par.c;if(par.c<0);error('wrong input: par.c');end
if ~isfield(par,'show_fval');par.show_fval=1;end;show_fval=par.show_fval;
if ~isfield(par,'detail');par.detail=1;end
if ~isfield(par,'test_time');par.test_time=0;end
if ~isfield(par,'maxiter');par.maxiter=[100,1000];end
maxiter_out=100;maxiter_in=1000;
iter_in_min=5;
%% initialization
cr_in=zeros(1,maxiter_out);
cr_x=zeros(1,maxiter_out);
cr_f=zeros(1,maxiter_out);
iter_in=zeros(1,maxiter_out);
con_Y=zeros(1,maxiter_out);
err=zeros(1,maxiter_out);
fval=zeros(1,maxiter_out);
t_all=zeros(1,maxiter_out);
Lre=zeros(1,maxiter_out);

X = zeros( m,n, h);
S = zeros( m,n, h);
%% DCA begin

for i=1:maxiter_out
    Xpout=X;
    X = fft(X,[],3);
    %% compute P
    t=round(h/2);
    [U_mc,S_temp,V_mc]=svd(X(:,:,t),'econ');
    if all(size(S_temp)-1)
      % [m,n,h]  =size(S_temp);
        % 
        S_temp=diag(S_temp);
    end
    S_temp=S_temp(:);
    r_temp=sum(S_temp>par.eps_min);%rank of X
    
    if r_temp<=t_now
        P_temp=ones(r_temp,1);
    else
        P_temp=[ones(t_now,1);...
            par.alpha*S_temp(t_now+1:r_temp)/norm(S_temp(t_now+1:r_temp),'fro')];
    end
    P=U_mc(:,1:r_temp)*diag(P_temp)*V_mc(:,1:r_temp)';
    P=P+c/lambda*X;
    %% ADMM
    for j=1:maxiter_in
        Xpin=X;
        
        Y = (par.mu*(M-L+Z/par.mu)+lambda*P+beta1*(X+S))/(par.mu+c+beta1);
        %% compute X
         X=SVT(Y-S,lambda/beta1,'L1');
          X = ifft(X,[],h);
        %% compute S
        X=reshape(X,[m,n,h]);
        S=S+X-Y;
        %% stop inner loop
        cr_in1=norm(Xpin-Xpin,'fro')/max(norm(Xpin,'fro'),par.eps_min);
        if cr_in1<tol_in
            stop_flag= (i<iter_out_min+1 || j>=iter_in_min);
            if i==1 && j<=50
                stop_flag=0;
            end
            if stop_flag
                break;
            end
        end
    end
    %% stop outer loop
    cr_in(i)=cr_in1;
    cr_x1=norm(X-Xpout,'fro')/max(norm(Xpout,'fro'),par.eps_min);
    cr_x(i)=cr_x1;
    stop_flag=cr_x1<tol_out && i>=iter_out_min+1;
    if show_fval
        fval(i)=fun_fval(M,X,L,Z,lambda,par,t_now);
        if i>=2
            cr_fval=abs(fval(i)-fval(i-1))/max(fval(i-1),par.eps_min);
        else
            cr_fval=1;
        end
        cr_f(i)=cr_fval;
        stop_flag=(cr_fval<tol_out || cr_x1<tol_out) && i>=iter_out_min+1;
    end
    iter_in(i)=j;
    con_Y(i)=norm(X-Y,'fro');
    err(i)=norm(M-X-L,'fro')/norm(M,'fro');
    t_all(i)=t_now;
    if par.detail
        if i==1
            % result_title='OutIt| InnIt|    ReErr|   OutCrX|   OutCrF|   InnCrX| Trun#';
            % disp(result_title)
        end
        % result_sprintf=sprintf('%5d| %5d| %8.2e| %8.2e| %8.2e| %8.2e| %5d',...
        %     i,iter_in(i),err(i),cr_x(i),cr_f(i),cr_in(i),t_all(i));
        % disp(result_sprintf)
    end
    if ~par.test_time%not record time
        if stop_flag
            break;
        end
    else%record time
        if norm(X-Xpout,'fro')/norm(Xpout,'fro')<par.test_thre
            %         if norm(X-par.x_t,'fro')/norm(par.x_t,'fro')<par.test_thre
            break;
        end
    end
end
end


%% SVT
function X=SVT(Z,tau,type)
%for computing SVT of a matrix
if nargin<3
    type='L1';
end
[m,n,h]=size(Z);
% sa=size(Z);
% la=length(sa);
% faces = prod(sa(3:la));
if m <= n
   
    % [m,n,h]=size(Z);
   z= reshape(Z,[],n*h/2);
   AAT =z*z';
    % end
    
    [S,V]=eig(AAT);
    V=max(0,diag(V));
    V=sqrt(V);
    tol=n*eps(max(V));
    mid=thresholdings(V,tau(end:-1:1),type);
    % mid=SoftThresh(V,tau(end:-1:1));
    ind=mid>tol;
    if any(ind)
        mid=mid(ind)./V(ind);
        
        X=S(:,ind)*diag(mid)*S(:,ind)'*z;
       
    else
        X=zeros(m,n,h);
    end
else
    X=SVT(Z,tau,type);
    X=X';
end
end

function x=thresholdings(y,tau,type)
%solve: min_x: 1/2(x-y)^2+tau*f(x)
%type:
% 'L1': f(x)=|x|
% 'L1/2': f(x)=sqrt(x)
if ~exist('type','var')
    type='L1';
end
switch type
    case 'L1'
        x=sign(y).*max(0,abs(y)-tau);
    case 'L1/2'
        ind=abs(y)>3/2*tau.^(2/3);
        x=zeros(size(y));
        y=y(ind);
        x(ind)=2/3*y.*(1+cos(2/3*(pi-acos(tau/4.*(abs(y)/3).^(-3/2)))));
end
end

function v=fun_fval(D, X, E, L, lambda, par, t_now)
X_temp=tsvd(X,'econ');
v1 = lambda*(TL1L2norm(X_temp,t_now,par.alpha));
v2 = par.mu/2*norm(D-X-E+L/par.mu,'fro')^2;
v = v1+v2;
end

function val=TL1L2norm(Z,t,alpha)
%compute truncated l_{1-2} metric
Z=sort(abs(Z),'descend');
val=sum(Z(t+1:end))-alpha*sqrt(sum(Z(t+1:end).^2));
end