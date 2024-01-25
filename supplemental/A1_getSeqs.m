function A1_getSeqs

% get chord seqs (C)

% inDir is the directory with our data %
inDir = '../McGill-Billboard1/';

% inDirs is the list of files and folders in our directory %
inDirs = dir([inDir '*']);
inDirs = inDirs(4:end);

% nFile is the number of files in our directory %
nFile = length(inDirs);

key = {'C' 'C#' 'Db' 'D' 'D#' 'Eb' 'E' 'F' 'F#' 'Gb' 'G' 'G#' 'Ab' ...
    'A' 'A#' 'Bb' 'B'};
val = [0 1 1 2 3 3 4 5 6 6 7 8 8 9 10 10 11];

revKey = {'C' 'C#' 'D' 'Eb' 'E' 'F' 'F#' 'G' 'Ab' 'A' 'Bb' 'B'};

% ids and seqs are both lists of length nFile %
ids = cell(1,nFile);
seqs = cell(1,nFile);


for i = 1:nFile
    [i nFile]
    
    % Concatenate elements of a file name to open %
    fin = [inDir inDirs(i).name '/majmin.lab'];
    
    % populate ith element of ids with the corresponding file name %
    ids{i} = inDirs(i).name;
    
    % fid is the text of the file %
    fid = fopen(fin,'r');
    
    while(1)
        % read a line from our file %
        ln = fgetl(fid);
        
        % exit the loop if the line is empty %
        if isempty(ln)
            break;
        end
        
        % sp is the indices of each horizontal tab %
        sp = strfind(ln,sprintf('\t'));
        
        % ch is a tabbed line after the tab %
        ch = ln(sp(2)+1:end);
        
        % cln is the indices of each : in a given ch %
        cln = strfind(ch,':');
        
        % if cln is empty, continue %
        if isempty(cln)
            continue;
        end
        
        % let is the chord name %
        let = ch(1:cln-1);
        
        % typ is the chord quality %
        typ = ch(cln+1:end);
        
        % idx is the pitch class of the chord? (revisit) %
        idx = find(strcmp(let,key));
        bin = -1;
        
        % set bin to 0 for major chords, 1 for minor chords %
        if strcmp(typ,'maj')
            bin = 0;
        else
            bin = 1;
        end
        
        % populate sequence with codified code value! %
        seqs{i} = [seqs{i} (12*bin)+val(idx)];
    end
    
    % close our file
    fclose(fid);    
end

save('data/seqs','ids','seqs');
