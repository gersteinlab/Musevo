function kmer_analysis4

% Filters from Eq. 2 (tau-normalized)

inDir = '../McGill-Billboard/';
inDirs = dir([inDir '*']);
inDirs = inDirs(4:end);
nFile = length(inDirs);
% nFile  = 10;
nList = 20;

key = {'C' 'C#' 'Db' 'D' 'D#' 'Eb' 'E' 'F' 'F#' 'Gb' 'G' 'G#' 'Ab' ...
    'A' 'A#' 'Bb' 'B' 'Cb'};
val = [0 1 1 2 3 3 4 5 6 6 7 8 8 9 10 10 11 11];

revKey = {'C' 'C#' 'D' 'Eb' 'E' 'F' 'F#' 'G' 'Ab' 'A' 'Bb' 'B'};

load('data/seqs','ids','seqs');

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

% get km histograms
nCh = 24;
km1 = zeros(nCh,2);
km2 = zeros(nCh,nCh,2);
km3 = zeros(nCh,nCh,nCh,2);
km4 = zeros(nCh,nCh,nCh,nCh,2);
for i = 1:nFile
    [i nFile]
    seq = seqs{i};
    if isempty(seq)
        continue;
    end
%     seq = seq + 1;
    reps = (seq(1:end-1)==seq(2:end));
    seq = seq([reps==0 true]);        
    for j = 1:length(seq)-1
%         idx1 = seq(j);   
%         idx2 = seq(j+1); 
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
        km1(idx2+1,typ1) = km1(idx2+1,typ1) + 1;
    end
    for j = 1:length(seq)-2
%         idx1 = seq(j);
%         idx2 = seq(j+1);
%         idx3 = seq(j+2);
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
        km2(vec(1),vec(2),typ1) = km2(vec(1),vec(2),typ1) + 1;
    end    
    for j = 1:length(seq)-3
%         idx1 = seq(j);
%         idx2 = seq(j+1);
%         idx3 = seq(j+2);
%         idx4 = seq(j+3);
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
        km3(vec(1),vec(2),vec(3),typ1) = km3(vec(1),vec(2),vec(3),typ1) + 1;
    end       
    for j = 1:length(seq)-4
%         idx1 = seq(j);
%         idx2 = seq(j+1);
%         idx3 = seq(j+2);
%         idx4 = seq(j+3);
%         idx5 = seq(j+4);
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
        km4(vec(1),vec(2),vec(3),vec(4),typ1) = km4(vec(1),vec(2),vec(3),vec(4),typ1) + 1;
    end       
end

save('data/kms_norm_tonic','km1','km2','km3','km4');

close all;
figure(1);
km1 = km1./sum(km1(:));
subplot(1,2,1);
km1a = km1(:,1)';
plot(0:11,km1a(1:12),'r-','linewidth',2); hold on;
plot(0:11,km1a(13:24),'b-','linewidth',2);
subplot(1,2,2);
km1a = km1(:,2)';
plot(0:11,km1a(1:12),'r-','linewidth',2); hold on;
plot(0:11,km1a(13:24),'b-','linewidth',2);

figure(2);
km2 = km2./sum(km2(:));
subplot(1,2,1);
km2a = km2(:,:,1);
imagesc(km2a); hold on;
colormap(hot);
subplot(1,2,2);
km2a = km2(:,:,2);
imagesc(km2a); hold on;
colormap(hot);

figure(3)
km2 = km2./sum(km2(:));
km2a = km2(:,:,1);
km2b = km2(:,:,2);

TwelveEdges = [0 pi/6 pi/3 pi/2 2*pi/3 5*pi/6 pi ...
    7*pi/6 4*pi/3 3*pi/2 5*pi/3 11*pi/6 2*pi]; 

for i = 1:24
subplot(4,12,i)

    polarhistogram('FaceColor','red', ...
        'BinEdges', TwelveEdges, ...
        'BinCounts', km2a(1:12,i), ...
        'FaceAlpha', 0.2, ... 
        'LineWidth', 1); hold on;
    
    polarhistogram('FaceColor','blue', ...
        'BinEdges', TwelveEdges, ...
        'BinCounts', km2a(13:24,i), ...
        'FaceAlpha', 0.2, ... 
        'LineWidth', 1) ; hold on;
    
    thetaticks([15 45 75 105 135 165 195 225 255 285 315 345]);
    thetaticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'});
    rticks([0.0]);
    rticklabels({});    
    pax = gca;
    pax.FontSize = 5;
    
end

for i = 1:24
subplot(4,12,i+24)
    polarhistogram('FaceColor','red', ...
        'BinEdges', TwelveEdges, ...
        'BinCounts', km2b(1:12,i), ...
        'FaceAlpha', 0.2, ... 
        'LineWidth', 1); hold on;
    
    polarhistogram('FaceColor','blue', ...
        'BinEdges', TwelveEdges, ...
        'BinCounts', km2b(13:24,i), ...
        'FaceAlpha', 0.2, ... 
        'LineWidth', 1)
    
    thetaticks([15 45 75 105 135 165 195 225 255 285 315 345]);
    thetaticklabels({'0','1','2','3','4','5','6','7','8','9','10','11'});
    rticks([0.0]);
    rticklabels({});
    pax = gca;
    pax.FontSize = 5;
    
end
sgtitle('Km=2 normalized to tonic')

annotation('textbox', [0 0.83 1 0.1], ...
    'String', 'major key context', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

annotation('textbox', [0 0.4 1 0.1], ...
    'String', 'minor key context', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')


annotation('textbox', [0.05 0.76 1 0.1], ...
    'String', 'major first chord', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

annotation('textbox', [0.05 0.54 1 0.1], ...
    'String', 'minor first chord', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')


annotation('textbox', [0.05 0.32 1 0.1], ...
    'String', 'major first chord', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

annotation('textbox', [0.05 0.1 1 0.1], ...
    'String', 'minor first chord', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

annotation('textbox', [0.05 0 1 0.1], ...
    'String', 'Root of first chord:', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

annotation('textbox', [0.148 0 1 0.1], ...
    'String', '0', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

annotation('textbox', [0.215 0 1 0.1], ...
    'String', '1', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

annotation('textbox', [0.28 0 1 0.1], ...
    'String', '2', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

annotation('textbox', [0.345 0 1 0.1], ...
    'String', '3', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

annotation('textbox', [0.41 0 1 0.1], ...
    'String', '4', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

annotation('textbox', [0.475 0 1 0.1], ...
    'String', '5', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

annotation('textbox', [0.545 0 1 0.1], ...
    'String', '6', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

annotation('textbox', [0.61 0 1 0.1], ...
    'String', '7', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

annotation('textbox', [0.675 0 1 0.1], ...
    'String', '8', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

annotation('textbox', [0.74 0 1 0.1], ...
    'String', '9', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

annotation('textbox', [0.805 0 1 0.1], ...
    'String', '10', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

annotation('textbox', [0.87 0 1 0.1], ...
    'String', '11', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'left')

% km1
fprintf(['km=1:\n']);
% nList = 10;
km1 = km1 ./ sum(km1(:));
[props idxs] = sort(km1(:),'descend');
for i = 1:nList
    idx = idxs(i);
    [i1 i0] = ind2sub([nCh*ones(1,1) 2],idx);
    if i0==1
        nm = 'C:maj ';
        fl = 0;
    else
        nm = 'A:min ';
        fl = 1;
    end
    vec = [i1];
    for j = 1:length(vec)
        if vec(j)>=12
            typ = 'min';
            vec(j) = vec(j) - 12;
        else
            typ = 'maj';
        end
        if fl==1
            if vec(j)>=13
                thres = 13;
            else
                thres = 1;
            end
            vec(j) = vec(j) - 3;
            if vec(j)<thres
                vec(j) = vec(j) + 12;
            end
        end
        nm = [nm revKey{vec(j)} ':' typ ' '];
    end
    nm = nm(1:end-1);
    fprintf([num2str(i) '\t' nm '\t' num2str(props(i)) '\n']);
end

% km2
fprintf(['\nkm=2:\n']);
% nList = 10;
km2 = km2 ./ sum(km2(:));
[props idxs] = sort(km2(:),'descend');
for i = 1:nList
    idx = idxs(i);
    [i1 i2 i0] = ind2sub([nCh*ones(1,2) 2],idx);
    if i0==1
        nm = 'C:maj ';
        fl = 0;
    else
        nm = 'A:min ';
        fl = 1;
    end
    vec = [i1 i2];
    for j = 1:length(vec)
        if vec(j)>=12
            typ = 'min';
            vec(j) = vec(j) - 12;
        else
            typ = 'maj';
        end
        if fl==1
            if vec(j)>=13
                thres = 13;
            else
                thres = 1;
            end
            vec(j) = vec(j) - 3;
            if vec(j)<thres
                vec(j) = vec(j) + 12;
            end
        end        
        nm = [nm revKey{vec(j)} ':' typ ' '];
    end    
    nm = nm(1:end-1);
    fprintf([num2str(i) '\t' nm '\t' num2str(props(i)) '\n']);
end

% km3
fprintf(['\nkm=3:\n']);
% nList = 10;
km3 = km3 ./ sum(km3(:));
[props idxs] = sort(km3(:),'descend');
for i = 1:nList
    idx = idxs(i);
    [i1 i2 i3 i0] = ind2sub([nCh*ones(1,3) 2],idx);
    if i0==1
        nm = 'C:maj ';
        fl = 0;
    else
        nm = 'A:min ';
        fl = 1;
    end
    vec = [i1 i2 i3];
    for j = 1:length(vec)
        if vec(j)>=12
            typ = 'min';
            vec(j) = vec(j) - 12;
        else
            typ = 'maj';
        end
        if fl==1
            if vec(j)>=13
                thres = 13;
            else
                thres = 1;
            end
            vec(j) = vec(j) - 3;
            if vec(j)<thres
                vec(j) = vec(j) + 12;
            end
        end        
        nm = [nm revKey{vec(j)} ':' typ ' '];
    end
    nm = nm(1:end-1);
    fprintf([num2str(i) '\t' nm '\t' num2str(props(i)) '\n']);
end

% km4
fprintf(['\nkm=4:\n']);
% nList = 10;
km4 = km4 ./ sum(km4(:));
[props idxs] = sort(km4(:),'descend');
for i = 1:nList
    idx = idxs(i);
    [i1 i2 i3 i4 i0] = ind2sub([nCh*ones(1,4) 2],idx);
    if i0==1
        nm = 'C:maj ';
        fl = 0;
    else
        nm = 'A:min ';
        fl = 1;
    end    
    vec = [i1 i2 i3 i4];
    for j = 1:length(vec)
        if vec(j)>=12
            typ = 'min';
            vec(j) = vec(j) - 12;
        else
            typ = 'maj';
        end
        if fl==1
            if vec(j)>=13
                thres = 13;
            else
                thres = 1;
            end
            vec(j) = vec(j) - 3;
            if vec(j)<thres
                vec(j) = vec(j) + 12;
            end
        end        
        nm = [nm revKey{vec(j)} ':' typ ' '];
    end
    nm = nm(1:end-1);
    fprintf([num2str(i) '\t' nm '\t' num2str(props(i)) '\n']);
end
