function structs_kmer_analysis

inDir = '../McGill-Billboard/';
inDirs = dir([inDir '*']);
inDirs = inDirs(4:end);
nFile = length(inDirs);

load('data/structs','structs');
load('data/projections','cat','projMat_A','cat_A','projMat_B','cat_B');


% get km histograms
nCh = 46;
km1 = zeros(nCh,1);
km2 = zeros(nCh,nCh,1);
km3 = zeros(nCh,nCh,nCh,1);
km4 = zeros(nCh,nCh,nCh,nCh,1);

for i = 1:nFile
    [i nFile];
    struct = structs{4,i};
    if isempty(struct)
        continue;
    end
        
    for j = 1:length(struct)
        idx1 = struct(j);    
        km1(idx1) = km1(idx1) + 1;
    end
    
    for j = 1:length(struct)-1
        idx1 = struct(j);  
        idx2 = struct(j+1);
        km2(idx1,idx2) = km2(idx1,idx2) + 1;
    end
    
    for j = 1:length(struct)-2
        idx1 = struct(j);  
        idx2 = struct(j+1);
        idx3 = struct(j+2);
        km3(idx1,idx2,idx3) = km3(idx1,idx2,idx3) + 1;
    end
        
    for j = 1:length(struct)-3
        idx1 = struct(j);  
        idx2 = struct(j+1);
        idx3 = struct(j+2);
        idx4 = struct(j+3);
        km4(idx1,idx2,idx3,idx4) = km4(idx1,idx2,idx3,idx4) + 1;
    end
          
end

projMat = 'None';

% % Apply Projection Matrix A
% 
%     cat = cat_A;
%     nCh_A = 18;
%     km1 = projMat_A * km1;
% 
%     for i = 1:nCh  
%     km2_redux1(i,:) = projMat_A * km2(:,i);
%     end 
% 
%     for i = 1:nCh_A
%     km2_redux2(i,:) = projMat_A * km2_redux1(:,i);
%     end  
%     km2 = km2_redux2;
% 
%     for i = 1:nCh
%         for j = 1:nCh
%             km3_redux1(i,j,:) = projMat_A * km3(:,i,j);
%         end
%     end
%     for i = 1:nCh
%         for j = 1:nCh_A
%             km3_redux2(i,j,:) = projMat_A * km3_redux1(:,i,j);
%         end
%     end
%     for i = 1:nCh_A
%         for j = 1:nCh_A
%             km3_redux3(i,j,:) = projMat_A * km3_redux2(:,i,j);
%         end
%     end
%     km3 = km3_redux3;
% 
%     for i = 1:nCh
%         for j = 1:nCh
%             for k = 1:nCh
%                 km4_redux1(i,j,k,:) = projMat_A * km4(:,i,j,k);
%             end
%         end
%     end 
%     for i = 1:nCh
%         for j = 1:nCh
%             for k = 1:nCh_A
%                 km4_redux2(i,j,k,:) = projMat_A * km4_redux1(:,i,j,k);
%             end
%         end
%     end 
%     for i = 1:nCh
%         for j = 1:nCh_A
%             for k = 1:nCh_A
%                 km4_redux3(i,j,k,:) = projMat_A * km4_redux2(:,i,j,k);
%             end
%         end
%     end
%     for i = 1:nCh_A
%         for j = 1:nCh_A
%             for k = 1:nCh_A
%                 km4_redux4(i,j,k,:) = projMat_A * km4_redux3(:,i,j,k);
%             end
%         end
%     end
%     km4 = km4_redux4;
% 
%     nCh = nCh_A;
%     projMat = 'A';
    
% % Apply Projection Matrix B
% 
%     cat = cat_B;
%     nCh_B = 17;
%     km1 = projMat_B * km1;
% 
%     for i = 1:nCh  
%     km2_redux1(i,:) = projMat_B * km2(:,i);
%     end   
% 
%     for i = 1:nCh_B
%     km2_redux2(i,:) = projMat_B * km2_redux1(:,i);
%     end  
%     km2 = km2_redux2;
% 
%     for i = 1:nCh
%         for j = 1:nCh
%             km3_redux1(i,j,:) = projMat_B * km3(:,i,j);
%         end
%     end
%     for i = 1:nCh
%         for j = 1:nCh_B
%             km3_redux2(i,j,:) = projMat_B * km3_redux1(:,i,j);
%         end
%     end
%     for i = 1:nCh_B
%         for j = 1:nCh_B
%             km3_redux3(i,j,:) = projMat_B * km3_redux2(:,i,j);
%         end
%     end
%     km3 = km3_redux3;
% 
%     for i = 1:nCh
%         for j = 1:nCh
%             for k = 1:nCh
%                 km4_redux1(i,j,k,:) = projMat_B * km4(:,i,j,k);
%             end
%         end
%     end 
%     for i = 1:nCh
%         for j = 1:nCh
%             for k = 1:nCh_B
%                 km4_redux2(i,j,k,:) = projMat_B * km4_redux1(:,i,j,k);
%             end
%         end
%     end 
%     for i = 1:nCh
%         for j = 1:nCh_B
%             for k = 1:nCh_B
%                 km4_redux3(i,j,k,:) = projMat_B * km4_redux2(:,i,j,k);
%             end
%         end
%     end
%     for i = 1:nCh_B
%         for j = 1:nCh_B
%             for k = 1:nCh_B
%                 km4_redux4(i,j,k,:) = projMat_B * km4_redux3(:,i,j,k);
%             end
%         end
%     end
%     km4 = km4_redux4;
% 
%     nCh = nCh_B;
%     projMat = 'B';  

save(['data/structs_kms_' projMat],'km1','km2','km3','km4','nCh','cat','projMat');

close all;

figure(1);
km1a = km1./sum(km1(:));
km1a = km1a(:)';
x = categorical(cat);
bar(x, km1a);

figure(2);
nList = 10;
x = categorical;
y = [];
[props idxs] = sort(km1a(:),'descend');
for i = 1:nList
    idx = idxs(i);
    x(i) = cat(idx);
    y(i) = props(i);
end  
bar(x,y);

figure(3);
km2a = km2./sum(km2(:));
imagesc(km2a);
colormap(hot);
xticks(1:46);
xticklabels(cat);
xtickangle(90);
yticks(1:46);
yticklabels(cat);

% Printouts

% Matrix tag
fprintf(['\nProjection Matrix: ' projMat '\n']);

% km1
fprintf(['\nkm=1:\n']);
nList = 10;
km1a = km1 ./ sum(km1(:));
[props idxs] = sort(km1a(:),'descend');
for i = 1:nList
    idx = idxs(i);
    symb1 = cat(idx);
    symb1 = string(symb1);
    
    fprintf([num2str(i) '\t' num2str(symb1) '\t' num2str(props(i)) '\n']);
end

% km2
fprintf(['\nkm=2:\n']);
nList = 16;
km2a = km2 ./ sum(km2(:));
[props idxs] = sort(km2a(:),'descend');
for i = 1:nList
    idx = idxs(i);
    [i1 i2] = ind2sub([nCh*ones(1,2)],idx);
    
    symb1 = cat(i1);
    symb2 = cat(i2);
    symb1 = string(symb1);
    symb2 = string(symb2);
    
    fprintf([num2str(i) '\t' num2str(symb1) ' : ' num2str(symb2) ...
        '\t' num2str(props(i)) '\n']);
end

% km3
fprintf(['\nkm=3:\n']);
nList = 4;
km3a = km3 ./ sum(km3(:));
[props idxs] = sort(km3a(:),'descend');
for i = 1:nList
    idx = idxs(i)';
    [i1 i2 i3] = ind2sub([nCh*ones(1,3)],idx);

    symb1 = cat(i1);
    symb2 = cat(i2);
    symb3 = cat(i3);
    symb1 = string(symb1);
    symb2 = string(symb2);
    symb3 = string(symb3);
    
    fprintf([num2str(i) '\t' num2str(symb1) ' : ' num2str(symb2) ...
        ' : ' num2str(symb3) '\t' num2str(props(i)) '\n']);
end

% km4
fprintf(['\nkm=4:\n']);
nList = 4;
km4 = km4 ./ sum(km4(:));
[props idxs] = sort(km4(:),'descend');
for i = 1:nList
    idx = idxs(i);
    [i1 i2 i3 i4] = ind2sub([nCh*ones(1,4)],idx);
    
    symb1 = cat(i1);
    symb2 = cat(i2);
    symb3 = cat(i3);
    symb4 = cat(i4);
    symb1 = string(symb1);
    symb2 = string(symb2);
    symb3 = string(symb3);
    symb4 = string(symb4);
    
    fprintf([num2str(i) '\t' num2str(symb1) ' : ' num2str(symb2) ...
        ' : ' num2str(symb3) ' : ' num2str(symb4) '\t' ...
        num2str(props(i)) '\n']);
end
