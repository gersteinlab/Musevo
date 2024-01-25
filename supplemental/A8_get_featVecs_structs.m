function A8_get_featVecs_structs

inDir = '../McGill-Billboard/';
inDirs = dir([inDir '*']);
inDirs = inDirs(4:end);
nFile = length(inDirs);

projMat

load('data/seqs','ids');
load('data/structs', 'structs');
load('data/structs_kms_' + projMat,'km1','km2','km3','km4','nCh','cat');

tb = readtable('billboard-2.0-index.csv');
dt = table2cell(tb(:,2));
dt2 = zeros(1,length(dt));
for i = 1:size(dt,1)
    dt2(i) = datenum(dt{i});
end
dt3 = zeros(1,nFile);
dt4 = cell(1,nFile);
for i = 1:nFile
    dt3(i) = dt2(str2num(ids{i}));
    dt4(i) = dt(str2num(ids{i}));
end
[dum idxs] = sort(dt3);
ids = ids(idxs);
structs = structs(idxs);
dates = dt4(idxs);


% get km histograms
km1 = zeros(nCh,nFile);
km2 = zeros(nCh,nCh,nFile);
km3 = zeros(nCh,nCh,nCh,nFile);
km4 = zeros(nCh,nCh,nCh,nCh,nFile,'uint8');
for i = 1:nFile
    [i nFile]
    struct = structs{i};
    if isempty(struct)
        continue;
    end
      
    for j = 1:length(struct)
        vec = strcmp(cat,struct(j));   
        for k = 1:nCh
            if vec(k) == 1
                idx1 = k;  
            end       
        end      
        km1(idx1,i) = km1(idx1,k) + 1;
    end
    
    for j = 1:length(struct)-1
        
        vec1 = strcmp(cat,struct(j)); 
        for k = 1:nCh
            if vec1(k) == 1
                idx1 = k;  
            end       
        end       
        vec2 = strcmp(cat,struct(j+1)); 
        for k = 1:nCh
            if vec2(k) == 1
                idx2 = k;  
            end       
        end
        km2(idx1,idx2,i) = km2(idx1,idx2,i) + 1;
    end 
    
    for j = 1:length(struct)-2
        
        vec1 = strcmp(cat,struct(j)); 
        for k = 1:nCh
            if vec1(k) == 1
                idx1 = k;  
            end       
        end      
        vec2 = strcmp(cat,struct(j+1)); 
        for k = 1:nCh
            if vec2(k) == 1
                idx2 = k;  
            end       
        end    
        vec3 = strcmp(cat,struct(j+2)); 
        for k = 1:nCh
            if vec3(k) == 1
                idx3 = k;  
            end       
        end
        km3(idx1,idx2,idx3,i) = km3(idx1,idx2,idx3,i) + 1;
    end 
    
    for j = 1:length(struct)-3
        
        vec1 = strcmp(cat,struct(j)); 
        for k = 1:nCh
            if vec1(k) == 1
                idx1 = k;  
            end       
        end      
        vec2 = strcmp(cat,struct(j+1)); 
        for k = 1:nCh
            if vec2(k) == 1
                idx2 = k;  
            end       
        end    
        vec3 = strcmp(cat,struct(j+2)); 
        for k = 1:nCh
            if vec3(k) == 1
                idx3 = k;  
            end       
        end
        vec4 = strcmp(cat,struct(j+3)); 
        for k = 1:nCh
            if vec4(k) == 1
                idx4 = k;  
            end       
        end
        km4(idx1,idx2,idx3,idx4,i) = km4(idx1,idx2,idx3,idx4,i) + 1;
    end
          
end

feat_key_struct = zeros(0,4);
X_struct = zeros(nFile,0);
cFeat = 1;

fprintf(['\nProjection Matrix: ' projMat '\n']);

display('1-km');

    for j1 = 1:nCh
        if sum(km1(j1,:))>0
            feat_key_struct(cFeat,:) = [j1 -1 -1 -1];
            X_struct(:,cFeat) = squeeze(km1(j1,:));
            cFeat = cFeat + 1;
        end
    end

display('2-km');

    for j1 = 1:nCh
        for j2 = 1:nCh
            if sum(km2(j1,j2,:))>0
                feat_key_struct(cFeat,:) = [j1 j2 -1 -1];
                X_struct(:,cFeat) = squeeze(km2(j1,j2,:));
                cFeat = cFeat + 1;
            end
        end
    end

display('3-km');

    for j1 = 1:nCh
        for j2 = 1:nCh
            for j3 = 1:nCh
                if sum(km3(j1,j2,j3,:))>0
                    feat_key_struct(cFeat,:) = [j1 j2 j3 -1];
                    X_struct(:,cFeat) = squeeze(km3(j1,j2,j3,:));
                    cFeat = cFeat + 1;
                end
            end
        end
    end

display('4-km');

    for j1 = 1:nCh
        for j2 = 1:nCh
            for j3 = 1:nCh
                for j4 = 1:nCh
                    if sum(km4(j1,j2,j3,j4,:))>0
                        feat_key_struct(cFeat,:) = [j1 j2 j3 j4];
                        X_struct(:,cFeat) = squeeze(km4(j1,j2,j3,j4,:));
                        cFeat = cFeat + 1;
                    end
                end
            end
        end
    end


save('data/structs_dataset_' + projMat,'X_struct','feat_key_struct','dates','ids');
