function B_visualizeModel_Structs_Binarized
% visualize learned model latent factors

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
fout = ['models/k1_dedup_structs_binarized_A'];

% load('data/structs_kms_B');
% load('data/structs_dataset_B','X_struct','feat_key_struct','ids');
% fout = ['models/k1_dedup_structs_binarized_B'];
% 
% load('data/structs_kms_None');
% load('data/structs_dataset_None','X_struct','feat_key_struct','ids');
% fout = ['models/k1_dedup_structs_binarized_None'];

trStep = 5; % steps betw tr samps
featThres = 50; % feat must occur in # songs

% load / prune / prepare dataset
load('models/k1_dedup','X_tr','X_te','nFeat','feat_key');
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
sum(sum(X_struct,2)==0);
nFeat

% nFeat = props(1);
% X_tr = X_tr(:,1:nFeat);
% X_te = X_te(:,1:nFeat);
% feat_key = feat_key(1:nFeat,:);

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

load(fout,'theta0','Q00','thetas','Qs','corrs','Ls','X_tr_struct','X_te_struct','feat_key_struct');

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
        [val S_idx_tr] = optFn(theta0,X_tr_struct,Q00);
        [val S_idx_te] = optFn(theta0,X_te_struct,Q00(1:N_te,:));
    else
        Q1 = ones(N,1) * mean(Q);        
        [val S_idx_te] = optFn(theta1,X_te_struct,Q1); 
        S_idx_tr = [];
        for n=1:N
            S_idx_tr = [S_idx_tr find(Q(n,:)==max(Q(n,:)),1)-1];
        end        
    end
    
    S_tr = Z_de2bi(S_idx_tr,Dz);
    S_te = Z_de2bi(S_idx_te,Dz);
%     if it==1
%         S_tr = Z_tr;
%         S_te = Z_te;
%     end
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

    S_trte = S_tr; %zeros(N,Dz);
%     S_trte(1:2:N,:) = S_tr;
%     S_trte(2:2:N,:) = S_te;
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
feats = feat_key_struct';
labels = cell(4,Dx);
for i = 1:(Dx*4)
    if feats(i) ~= -1
        labels(i) = cat(feats(i));
    else
        labels(i) = {' '};
    end
end

for i = 1:4:(Dx*4)
    labels(i) = {string(labels(i)) + ' ' + string(labels(i+1)) + ...
       ' ' + string(labels(i+2)) + ' ' + string(labels(i+3))};
end
labels = labels(1,:);

% plot combined features
theta = theta1;
W = reshape(theta(1:Dz*Dx),[Dz,Dx]);
theta = theta(Dz*Dx+1:end);
b = theta(1:Dx)';
theta = theta(Dx+1:end);
W = W + ones(Dz,1)*b;

figure(199);
for i = 1:Dz
    vec = W(i,:); 
    subplot(Dz,1,i);
    plot(vec(:),'r-','linewidth',2);
    xticklabels(labels); 
    xticks(1:34); 
    xtickangle(30);
end

% plot feature profiles
% theta = theta1;
% W = reshape(theta(1:Dz*Dx),[Dz,Dx]);
% theta = theta(Dz*Dx+1:end);
% b = theta(1:Dx)';
% theta = theta(Dx+1:end);
% W = W + ones(Dz,1)*b;
% figure(200);
% for i = 1:Dz
%     vec = W(i,:);
% %     vec2 = ones(nCh,1) * (min(vec)-1);
%     vec2 = ones(nCh,1) * NaN;
%     for j = 1:nFeat
%         if feat_key_struct(j,2) == -1
%             vec2(feat_key_struct(j,1)) = vec(j);
%         end
%     end
%     subplot(Dz,1,i);
%     sgtitle('Km=1')
%     plot(vec2(isfinite(vec2)),'r-','linewidth',2);
%     
%     xt = [];
%     mxs = islocalmax(vec2);
%     for k = 1:size(mxs)
%         if mxs(k) == 1
%             xt = [xt k];
%         end
%     end
%     xticks(xt);   
% end
% 
% figure(201);
% for i = 1:Dz
%     vec = W(i,:);
% %     vec2 = ones(nCh,nCh) * (min(vec)-1);
%     vec2 = ones(nCh,nCh) * NaN;
%     for j = 1:nFeat
%         if feat_key_struct(j,2) ~= -1 && feat_key_struct(j,3) == -1
%                 vec2(feat_key_struct(j,1),feat_key_struct(j,2)) = vec(j);
%         end
%     end
%     subplot(Dz,1,i);
%     sgtitle('Km=2');
%     plot(vec2(isfinite(vec2)),'r-','linewidth',2);
%     
% %     xt = [];
% %     mxs = islocalmax(vec2(:));
% %     for k = 1:size(mxs)
% %         if mxs(k) == 1
% %             xt = [xt k];
% %         end
% %     end
% %     xticks(xt); 
% %     
% %     vec3 = [];
% %     xl = xticklabels;
% %     for k = 1:size(xl)
% %         [a b] = ind2sub([nCh nCh],str2double(xl(k)));
% %         vec3(k) = '(' + string(a) + ', ' + string(b) + ')';
% %     end
% %     xticklabels(vec3)
% %     xtickangle(90)
% end
% 
% figure(202);
% for i = 1:Dz
%     vec = W(i,:);
%     vec2 = ones(nCh,nCh,nCh) * (min(vec)-1);
%     for j = 1:nFeat
%         if feat_key_struct(j,3) ~= -1 && feat_key_struct(j,4) == -1
%                 vec2(feat_key_struct(j,1),feat_key_struct(j,2), ...
%                     feat_key_struct(j,3)) = vec(j);
%         end
%     end
%     subplot(Dz,1,i);
%     sgtitle('Km=3')
%     plot(vec2(:),'r-','linewidth',2);
%     
%     xt = [];
%     mxs = islocalmax(vec2(:));
%     for k = 1:size(mxs)
%         if mxs(k) == 1
%             xt = [xt k];
%         end
%     end
%     xticks(xt); 
%     
%     xl = xticklabels;
%     for k = 1:size(xl)
%         [a b c] = ind2sub([nCh nCh nCh],str2double(xl(k)));
%         vec3(k) = '(' + string(a) + ', ' + string(b) + ', ' + string(c) + ')';
%     end
%     xticklabels(vec3)
%     xtickangle(90)
% end
% 
% figure(203);
% for i = 1:Dz
%     vec = W(i,:);
%     vec2 = ones(nCh,nCh,nCh,nCh) * (min(vec)-1);
%     for j = 1:nFeat
%         if feat_key_struct(j,4) ~= -1
%                 vec2(feat_key_struct(j,1),feat_key_struct(j,2), ...
%                     feat_key_struct(j,3),feat_key_struct(j,4)) = vec(j);
%         end
%     end
%     subplot(Dz,1,i);
%     sgtitle('Km=4')
%     plot(vec2(:),'r-','linewidth',2);
% 
%     xt = [];
%     mxs = islocalmax(vec2(:));
%     for k = 1:size(mxs)
%         if mxs(k) == 1
%             xt = [xt k];
%         end
%     end
%     xticks(xt); 
%     
%     xl = xticklabels;
%     for k = 1:size(xl)
%         [a b c d] = ind2sub([nCh nCh nCh nCh],str2double(xl(k)));
%         vec3(k) = '(' + string(a) + ', ' + string(b) + ', ' ...
%             + string(c) + ', ' + string(d) + ')';
%     end
%     xticklabels(vec3)
%     xtickangle(90)
% end


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