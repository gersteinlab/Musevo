function A_trainModel_4mer

% prune duplicates

global lls_tr lls_te X_tr X_te S Nh Nw;
rng(100);

fin = 'data/dataset_4_100';
fout = 'models/model_4_100';

Dx = nan; % set from dataset
Dz = 5;
N = nan; % set from dataset
Ns = 2^Dz;
Nh = 1; Nw = 5; % NN height & width
ew = 4; % evolutionary window size
ec = 1; % evolutionary coupling

trStep = 5; % steps betw tr samps

% load / prune / prepare dataset
load(fin,'X','feat_key');
% preproc
nFile = size(X,1);
teIdx = trStep:trStep:nFile;
trIdx = setdiff(1:nFile,teIdx);
X(:) = X(:) > 0;
X_tr = X(trIdx,:);
X_te = X(teIdx,:);
nFeat = size(feat_key,2);

%%%%%%%%%%%% train evolutionary model

Dx = nFeat;
N = size(X_tr,1);
N_te = size(X_te,1);
S = Z_de2bi(0:Ns-1);

if Nh==1
    Ntheta = Dz*Dx + Dx;
else
    Ntheta = Dz*Nw + Nw + (Nh-2) * (Nw*Nw + Nw) + Nw*Dx + Dx;
end
theta0 = randn(Ntheta,1);
Q00 = ones(N,Ns)*(1/Ns);
Q = Q00;
Ls = [];
corrs = [];
J = zeros(N,N);
for n = 1:N
    n2 = [n-ew:n+ew];
    n2 = n2(n2>0);
    n2 = n2(n2<=N);
    J(n,n2) = 1;
end
V = zeros(Ns,Ns);
for i = 1:Ns
    for j = 1:Ns
        h = sum(S(i,:)~=S(j,:));
        V(i,j) = exp(-h*ec);
    end
end
thetas = [];
Qs = [];
for it = 1:10

    % optimize thetas
    fn = @(theta) -optFn(theta,X_tr,Q);
    if it==1
        theta_in = theta0;
    else
        theta_in = theta1;
    end
    opts = optimset('Display','iter','MaxIter',5,...
        'MaxFunEvals',10000);
    [theta1] = fminunc(fn,theta_in,opts);
    [val S_idx_tr U] = optFn(theta1,X_tr,Q);
    
    % optimize Q
    for it2=1:(N)*100
        n = ceil(rand*N);              
        vec = U(n,:);
        idx = find(J(n,:)==1);
        for i = 1:length(idx)
            vec = vec + Q(idx(i),:) * V;
        end
        vec = exp(vec);
        vec = vec ./ sum(vec);           
        Q(n,:) = vec;        
    end
    
    % calc ELBO bound
    L = optFn(theta1,X_tr,Q);
    for n = 1:N
        idx = find(J(n,:)==1);
        idx = idx(idx>n);
        for i = 1:length(idx)
            L = L + (Q(idx(i),:) * V) * Q(n,:)';
        end    
        L = L - sum(Q(n,:).*log(Q(n,:)));
    end
    
    Ls = [Ls L];
    display('ELBOs:');
    Ls    
    
    thetas = [thetas theta1];
    Qs = cat(3,Qs,Q);
    
    S_idx_tr = [];
    for n=1:N
        S_idx_tr = [S_idx_tr find(Q(n,:)==max(Q(n,:)),1)-1];
    end
    Q1 = ones(N,1) * mean(Q);
    [val S_idx_te] = optFn(theta1,X_te,Q1);
    S_tr = Z_de2bi(S_idx_tr,Dz);
    S_te = Z_de2bi(S_idx_te,Dz);
    y_est = zeros(N_te,1);    
    beta = 0.1;
    for i = 1:N_te
        h = (ones(N,1) * S_te(i,:)) ~= S_tr;
        h = sum(h,2);
        w = exp(-h*beta);
        w = w ./ sum(w);
        y_est(i) = sum(w.*(1:N)');
    end    
    corr = corrcoef(y_est,1:N_te);
    corr = corr(1,2);
    corrs = [corrs corr];
    display('corrs:');
    corrs
    pause(0.5);
 
end

save(fout,'theta0','Q00','thetas','Qs','corrs','Ls','X_tr','X_te','feat_key');


function [val S_idx U] = optFn(theta,X,Q)

global S sig Nh Nw;

[N Dx] = size(X);
[Ns Dz] = size(S);

a = S;

for h = 1:Nh
    if h==1 && Nh>1
        W = reshape(theta(1:Dz*Nw),[Dz,Nw]);
        theta = theta(Dz*Nw+1:end);
        b = theta(1:Nw)';
        theta = theta(Nw+1:end);
    elseif h==1 && Nh==1
        W = reshape(theta(1:Dz*Dx),[Dz,Dx]);
        theta = theta(Dz*Dx+1:end);
        b = theta(1:Dx)';
        theta = theta(Dx+1:end);        
    elseif h<Nh
        W = reshape(theta(1:Nw*Nw),[Nw,Nw]);
        theta = theta(Nw*Nw+1:end);
        b = theta(1:Nw)';
        theta = theta(Nw+1:end); 
    else
        W = reshape(theta(1:Nw*Dx),[Nw,Dx]);
        theta = theta(Nw*Dx+1:end);
        b = theta(1:Dx)';
        theta = theta(Dx+1:end);         
    end
    a = a*W + ones(Ns,1)*b;
    a = 1 ./ (1+exp(-a));

end

ll = 0;

S_idx = [];
U = zeros(N,Ns);
for n = 1:N
%     cll = lognpdf(ones(Ns,1)*X(n,:),a,sig*ones(Ns,Dx));
    cll = log(binopdf(ones(Ns,1)*X(n,:),ones(Ns,Dx),a));
    cll = sum(cll,2);
    U(n,:) = cll';
    ll = ll + sum(cll.*Q(n,:)');
    S_idx = [S_idx find((cll.*Q(n,:))==max(cll.*Q(n,:)),1)-1];
end

val = ll;

