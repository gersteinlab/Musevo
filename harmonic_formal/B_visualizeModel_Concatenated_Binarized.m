function B_visualizeModel_Concatenated_Binarized

global lls_tr lls_te X_tr_struct X_te_struct S Nh Nw;
rng(100);
close all;

Dx = nan; % set from dataset
Dz = 5;
N = nan; % set from dataset
Ns = 2^Dz;
Nh = 1; Nw = 5; % NN height & width
ew = 4; % evolutionary window size
ec = 1; % evolutionary coupling

% choose model
load('data/structs_kms_A');
load('data/structs_dataset_A','X_struct','feat_key_struct','ids');
fout = ['models/k1_dedup_concatenated_binarized_A'];

% load('data/structs_kms_B');
% load('data/structs_dataset_B','X_struct','feat_key_struct','ids');
% fout = ['models/k1_dedup_concatenated_binarized_B'];
% 
% load('data/structs_kms_None');
% load('data/structs_dataset_None','X_struct','feat_key_struct','ids');
% fout = ['models/k1_dedup_concatenated_binarized_None'];

trStep = 5; % steps betw tr samps
featThres = 40; % feat must occur in # songs

% load / prune / prepare dataset
load('models/k1_dedup','X_tr','X_te','nFeat','feat_key');
cats = [cat '1' '2' '3' '4' '5' '6' '7' '8' '9' ...
    '10' '11' '12' '13' '14' '15' '16' '17' '18' ...
    '19' '20' '21' '22' '23' '24'];

nFeat_1 = nFeat;

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

% binarize X_struct
X_struct(:) = X_struct(:) > 0;
X_struct_bin = X_struct;

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

nFeat = nFeat + nFeat_1;
nFeat

% read in X_tr, X_te
X_tr_all = [X_tr X_tr_struct];
X_te_all = [X_te X_te_struct];


%%%%%%%%%%%% train evolutionary model

Dx = nFeat;
N = size(X_tr_all,1);
N_te = size(X_te_all,1);
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

load(fout,'theta0','Q00','thetas','Qs','corrs','Ls','X_tr_all','X_te_all','feat_key_struct');

figure(1000);
subplot(1,2,1);
plot(Ls,'b-','LineWidth',2);
subplot(1,2,2);
plot(corrs,'b-','LineWidth',2);

% test representations
idx = find(corrs==max(corrs),1);
Q = Qs(:,:,idx);
theta1 = thetas(:,idx);
for it = 1:2
    
    if it==1
        [val S_idx_tr] = optFn(theta0,X_tr_all,Q00);
        [val S_idx_te] = optFn(theta0,X_te_all,Q00(1:N_te,:));
    else
        Q1 = ones(N,1) * mean(Q);        
        [val S_idx_te] = optFn(theta1,X_te_all,Q1); 
        S_idx_tr = [];
        for n=1:N
            S_idx_tr = [S_idx_tr find(Q(n,:)==max(Q(n,:)),1)-1];
        end        
    end
    
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
%     display(['corrcoef: ' num2str(corr)]);

    S_trte = S_tr; 
    figure(it);
    for i = 1:4
        subplot(1,4,i);
        mat = S_trte((i-1)*200+1:min(i*200,N),:);
        if size(mat,1)<200
            mat = [mat ; zeros(200-size(mat,1),size(mat,2))];
        end
        imshow(kron(mat,ones(100,100)));
    end    
end

% track feats throughout time
tw = 40; % time window
avs = zeros(N,Dz);
for n = 1:N
    n2 = [n-tw:n+tw];
    n2 = n2(n2>0);
    n2 = n2(n2<=N);
    avs(n,:) = mean(S_tr(n2,:));
end
figure(100);
cols = {'r' 'g' 'b' 'm' 'c'};
for i = 1:Dz
    plot(avs(:,i),[cols{i} '-'],'LineWidth',2);
    hold on;
end
figure(101);
cols = {'r' 'g' 'b' 'm' 'c'};
avs = avs ./ (mean(avs,2)*ones(1,Dz));
avs = avs ./ (ones(N,1)*mean(avs,1));
for i = 1:Dz
    plot(avs(:,i),[cols{i} '-'],'LineWidth',2);
    hold on;
end

% convert numeric labels to words
feat_key_trunc = feat_key;
feat_key_trunc(:,5) = [];
feats = [feat_key_struct' feat_key_trunc'];
labels = cell(4,Dx);
for i = 1:((Dx-nFeat_1)*4)+1
    if feats(i) ~= -1
        labels(i) = cats(feats(i));
    else
        labels(i) = {' '};
    end
end

j = 1;
for i = 1+((Dx-nFeat_1)*4):(Dx*4)
    if feats(i) ~= -1
        labels(i) = {feats(i)-1};
        if feat_key(j,5) == 1
            labels(i) = {string(labels(i)) + ' maj'};
        else
            labels(i) = {string(labels(i)) + ' min'};
        end
        
        if j < length(feat_key)
            j = j+1;
        end
    else
        labels(i) = {'-'};
    end
end

for i = 1:4:(Dx*4)
    labels(i) = {string(labels(i)) + ' ' + string(labels(i+1)) + ...
       ' ' + string(labels(i+2)) + ' ' + string(labels(i+3))};
end

labels = labels(1,:);

% plot combined struct + harmonic features
theta = theta1;
W = reshape(theta(1:Dz*Dx),[Dz,Dx]);
theta = theta(Dz*Dx+1:end);
b = theta(1:Dx)';
theta = theta(Dx+1:end);

% Z-score data normalization
W = W - (ones(Dz,1)*mean(W)); 
W = W ./ (ones(Dz,1)*std(W));

% JW's more efficient code
[mat idx] = sort(W);
for i = 1:size(W,2)
    W(idx(:,i),i) = 1:size(W,1);
end

figure(199);
for i = 1:Dz
    vec = W(i,:); 
    subplot(Dz,1,i);
    plot(vec(:),'r-','linewidth',2);
    xticklabels(labels); 
    xticks(1:Dx); 
    xtickangle(30);
end


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