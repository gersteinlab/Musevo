function A_trainModel_Structs_Zscore

% Train evolutionary model

% OUTPUT:

% ELBOs:
% 
% Ls =
% 
%    1.0e+03 *
% 
%   Columns 1 through 9
% 
%     2.3013    2.3744    2.3761    2.3761    2.3761    2.3761    2.3761    2.3761    2.3761
% 
%   Column 10
% 
%     2.3761
% 
% corrs:
% 
% corrs =
% 
%      0     0     0     0     0     0     0     0     0     0

tic
global lls_tr lls_te X_tr_struct X_te_struct S Nh Nw;

% choose model
load('data/structs_kms_A');
load('data/structs_dataset_A','X_struct','feat_key_struct','ids');
fout = ['models/k1_dedup_structs_zscore_A'];
fout_best = ['models/model_best_structs_zscore_A'];

% load('data/structs_kms_B');
% load('data/structs_dataset_B','X_struct','feat_key_struct','ids');
% fout = ['models/k1_dedup_structs_zscore_B'];
% fout = ['models/k1_dedup_structs_zscore_B'];
% fout_best = ['models/model_best_structs_zscore_B'];
% 
% load('data/structs_kms_None');
% load('data/structs_dataset_None','X_struct','feat_key_struct','ids');
% fout = ['models/k1_dedup_structs_zscore_None'];
% fout = ['models/k1_dedup_structs_zscore_B'];
% fout_best = ['models/model_best_structs_zscore_B'];

rng(100);

Dx = nan; % set from dataset
Dz = 5;
N = nan; % set from dataset
Ns = 2^Dz;
Nh = 1; Nw = 5; % NN height & width
ew = 4; % evolutionary window size
ec = 1; % evolutionary coupling

trStep = 5; % steps betw tr samps
featThres = 50; % feat must occur in # songs

% prune / prepare dataset
ls = [];
for i = 1:length(ids)
    ls = [ls str2num(ids{i})];
end
ids0 = ls;
% remove dups
load('data/ids_dedup','ids');
ls = [];
for i = 1:length(ids)
    ls = [ls str2num(ids{i})];
end
idx = ismember(ids0,ls);
X_struct = X_struct(idx,:);

% preproc
nFile = size(X_struct,1);
teIdx = trStep:trStep:nFile;
trIdx = setdiff(1:nFile,teIdx);

% NORMALIZATIONS
X_struct_bin = X_struct;
X_struct_bin(:) = X_struct_bin(:) > 0;
X_struct = X_struct - ones(nFile,1)*mean(X_struct); % Z-score features
X_struct = X_struct ./ ones(nFile,1).*std(X_struct);
 
X_tr_struct = X_struct(trIdx,:);
X_tr_struct_bin = X_struct_bin(trIdx,:);
X_te_struct = X_struct(teIdx,:);
ct = sum(X_tr_struct_bin);

selFeats = ct>=featThres;
X_struct = X_struct(:,selFeats);
X_tr_struct = X_tr_struct(:,selFeats);
X_te_struct = X_te_struct(:,selFeats);
feat_key_struct = feat_key_struct(selFeats,:);
nFeat = size(feat_key_struct,1);
% feat props
props = [sum(feat_key_struct(:,2:4)==-1) nFeat];
props(2:4) = props(2:4) - props(1:3);
props
% data points
sum(sum(X_struct,2)==0)
nFeat

% nFeat = props(1);
% X_tr_struct = X_tr_struct(:,1:nFeat);
% X_te_struct = X_te_struct(:,1:nFeat);
% feat_key_struct = feat_key_struct(1:nFeat,:);
save('datMats','X_tr_struct','X_te_struct');
%%%%%%%%%%%% train evolutionary model

Dx = nFeat;
N = size(X_tr_struct,1);
N_te = size(X_te_struct,1);
S = Z_de2bi(0:Ns-1);

% train DVAE
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
    fn = @(theta) -optFn(theta,X_tr_struct,Q);
    if it==1
        theta_in = theta0;
    else
        theta_in = theta1;
    end
    opts = optimset('Display','iter','MaxIter',5,...
        'MaxFunEvals',10000);
    [theta1] = fminunc(fn,theta_in,opts);
    [val S_idx_tr U] = optFn(theta1,X_tr_struct,Q);
    
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
    L = optFn(theta1,X_tr_struct,Q);
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
    [val S_idx_te] = optFn(theta1,X_te_struct,Q1);
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
    
    if corrs(end)==max(corrs)
        save(fout_best,'theta0','thetas','Qs','corrs','Ls','X_tr_struct','X_te_struct','feat_key_struct');
    end
    
    toc
 
end
fprintf(['\nProjection Matrix: ' projMat '\n']);

save(fout,'theta0','Q00','thetas','Qs','corrs','Ls','X_tr_struct','X_te_struct','feat_key_struct');


function [val S_idx U] = optFn(theta,X_struct,Q)

global S sig Nh Nw;

[N Dx] = size(X_struct);
[Ns Dz] = size(S);
sig = 0.1; % sig = 10

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
%     if h<Nh
        a = 1 ./ (1+exp(-a));
%     end
end

ll = 0;

S_idx = [];
U = zeros(N,Ns);
for n = 1:N
    cll = lognpdf(ones(Ns,1)*X_struct(n,:),a,sig*ones(Ns,Dx));
%     cll = log(binopdf(ones(Ns,1)*X_struct(n,:),ones(Ns,Dx),a));
    cll = sum(cll,2);
    U(n,:) = cll';
    ll = ll + sum(cll.*Q(n,:)');
    S_idx = [S_idx find((cll.*Q(n,:))==max(cll.*Q(n,:)),1)-1];
end

val = ll;

