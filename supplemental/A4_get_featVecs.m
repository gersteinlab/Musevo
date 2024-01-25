function A4_get_featVecs

% Using filters from Eq. 2 (tau-normalized)

inDir = '../McGill-Billboard/';
inDirs = dir([inDir '*']);
inDirs = inDirs(4:end);
nFile = length(inDirs);

key = {'C' 'C#' 'Db' 'D' 'D#' 'Eb' 'E' 'F' 'F#' 'Gb' 'G' 'G#' 'Ab' ...
    'A' 'A#' 'Bb' 'B' 'Cb'};
val = [0 1 1 2 3 3 4 5 6 6 7 8 8 9 10 10 11 11];

revKey = {'C' 'C#' 'D' 'Eb' 'E' 'F' 'F#' 'G' 'Ab' 'A' 'Bb' 'B'};

load('data/seqs','ids','seqs');
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
seqs = seqs(idxs);
dates = dt4(idxs);

% get keys
music_keys = nan * ones(1,nFile);
for i = 1:nFile
    [i nFile]
    seq = seqs{i};
    if isempty(seq)
        continue;
    end
    reps = (seq(1:end-1)==seq(2:end));
    seq = seq([reps==0 true]); 
    fid = fopen(['../McGill-Billboard/' ids{i} '/salami_chords.txt']);
    for j = 1:10
        ln = fgetl(fid); 
        if strcmp(ln(1:7),'# tonic')
            break;
        end    
    end
    fclose(fid);
    if length(ln)==10
        k = ln(end);
    else
        k = ln(end-1:end);
    end
    v = val(strcmp(key,k));
    h1 = sum(seq==v);
    h2 = sum(seq==(v+12));
    if h2>h1
        v = v + 12;
    end
    music_keys(i) = v;
end

save('data/music_keys','music_keys','ids');

% get km histograms
nCh = 24;
km1 = zeros(nCh,2,nFile);
km2 = zeros(nCh,nCh,2,nFile);
km3 = zeros(nCh,nCh,nCh,2,nFile);
km4 = zeros(nCh,nCh,nCh,nCh,2,nFile);
for i = 1:nFile
    [i nFile]
    seq = seqs{i};
    if isempty(seq)
        continue;
    end
    reps = (seq(1:end-1)==seq(2:end));
    seq = seq([reps==0 true]);        
    for j = 1:length(seq)-1
        idx1 = music_keys(i);   
        idx2 = seq(j);         
        if idx1>=12
            typ1 = 2;
            idx1 = idx1 - 12;
        else
            typ1 = 1;
        end
        if idx2>=12
            thres = 12;
        else
            thres = 0;
        end        
        idx2 = idx2 - idx1;
        if idx2<thres
            idx2 = idx2 + 12;
        end
        km1(idx2+1,typ1,i) = km1(idx2+1,typ1,i) + 1;
    end
    for j = 1:length(seq)-2
        idx1 = music_keys(i);
        idx2 = seq(j);
        idx3 = seq(j+1);
        vec = [idx2 idx3];
        if idx1>=12
            typ1 = 2;
            idx1 = idx1 - 12;
        else
            typ1 = 1;
        end
        for k = 1:length(vec)        
            if vec(k)>=12
                thres = 12;
            else
                thres = 0;
            end
            vec(k) = vec(k) - idx1;
            if vec(k)<thres
                vec(k) = vec(k) + 12;
            end
        end
        vec = vec + 1;
        km2(vec(1),vec(2),typ1,i) = km2(vec(1),vec(2),typ1,i) + 1;
    end    
    for j = 1:length(seq)-3
        idx1 = music_keys(i);
        idx2 = seq(j);
        idx3 = seq(j+1);
        idx4 = seq(j+2);        
        vec = [idx2 idx3 idx4];
        if idx1>=12
            typ1 = 2;
            idx1 = idx1 - 12;
        else
            typ1 = 1;
        end
        for k = 1:length(vec)        
            if vec(k)>=12
                thres = 12;
            else
                thres = 0;
            end
            vec(k) = vec(k) - idx1;
            if vec(k)<thres
                vec(k) = vec(k) + 12;
            end
        end
        vec = vec + 1;
        km3(vec(1),vec(2),vec(3),typ1,i) = km3(vec(1),vec(2),vec(3),typ1,i) + 1;
    end       
    for j = 1:length(seq)-4
        idx1 = music_keys(i);
        idx2 = seq(j);
        idx3 = seq(j+1);
        idx4 = seq(j+2);
        idx5 = seq(j+3);
        vec = [idx2 idx3 idx4 idx5];
        if idx1>=12
            typ1 = 2;
            idx1 = idx1 - 12;
        else
            typ1 = 1;
        end
        for k = 1:length(vec)        
            if vec(k)>=12
                thres = 12;
            else
                thres = 0;
            end
            vec(k) = vec(k) - idx1;
            if vec(k)<thres
                vec(k) = vec(k) + 12;
            end
        end
        vec = vec + 1;
        km4(vec(1),vec(2),vec(3),vec(4),typ1,i) = km4(vec(1),vec(2),vec(3),vec(4),typ1,i) + 1;
    end       
end

feat_key = zeros(0,5);
X = zeros(nFile,0);
cFeat = 1;
display('1-km');
for i = 1:2
    for j1 = 1:nCh
        if sum(km1(j1,i,:))>0
            feat_key(cFeat,:) = [j1 -1 -1 -1 i];
            X(:,cFeat) = squeeze(km1(j1,i,:));
            cFeat = cFeat + 1;
        end
    end
end
display('2-km');
for i = 1:2
    for j1 = 1:nCh
        for j2 = 1:nCh
            if sum(km2(j1,j2,i,:))>0
                feat_key(cFeat,:) = [j1 j2 -1 -1 i];
                X(:,cFeat) = squeeze(km2(j1,j2,i,:));
                cFeat = cFeat + 1;
            end
        end
    end
end
display('3-km');
for i = 1:2
    for j1 = 1:nCh
        for j2 = 1:nCh
            for j3 = 1:nCh
                if sum(km3(j1,j2,j3,i,:))>0
                    feat_key(cFeat,:) = [j1 j2 j3 -1 i];
                    X(:,cFeat) = squeeze(km3(j1,j2,j3,i,:));
                    cFeat = cFeat + 1;
                end
            end
        end
    end
end
display('4-km');
for i = 1:2
    for j1 = 1:nCh
        for j2 = 1:nCh
            for j3 = 1:nCh
                for j4 = 1:nCh
                    if sum(km4(j1,j2,j3,j4,i,:))>0
                        feat_key(cFeat,:) = [j1 j2 j3 j4 i];
                        X(:,cFeat) = squeeze(km4(j1,j2,j3,j4,i,:));
                        cFeat = cFeat + 1;
                    end
                end
            end
        end
    end
end

save('data/dataset','X','feat_key','dates','ids');
