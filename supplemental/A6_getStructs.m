function A6_getStructs

inDir = '../McGill-Billboard/';
inDirs = dir([inDir '*']);
inDirs = inDirs(4:end);
nFile = length(inDirs);

structs = cell(4,nFile);

load('data/projections','cat');

for i = 1:nFile
    [i nFile]
    
    symbs = cell(0,0);
    timestamps = cell(0,0);
    letters = cell(0,0);
    
    nCh = 46;
    strseq = [];
    
    fin = [inDir inDirs(i).name '/salami_chords.txt'];
    ids{i} = inDirs(i).name;
    fid = fopen(fin,'r');
    
    for j = 1:5
        ln = fgetl(fid);
    end
    
    while(1)
        ln = fgetl(fid);
        
        if isempty(ln) | ln==-1
            break;
        end
        
        sp = strfind(ln,sprintf('\t'));
        
        if ln(sp+1)=='|'
            continue;
        end
        
        timestamp = ln(1:sp);
        ln = ln(sp+1:end);
        com = strfind(ln,',');
        
        if isempty(com)
            sym = ln;
            letter = '';
            
        elseif length(com)==1
            sym = ln(1:com(1)-1);
            letter = '';
            
        else
            sym = ln(com(1)+2:com(2)-1);
            letter = ln(com(1)-1);
            
        end
        
        idx = find(strcmp(sym, cat));
        strseq = [strseq idx];
        
        if isempty(sym)
            continue;
            
        elseif sym(1)=='|'
            continue;
            
        elseif length(sym)==1
            continue;
            
        elseif sym(1)=='Z'
            continue;
            
        else
            symbs{end+1} = sym;
            timestamps{end+1} = timestamp;
            letters{end+1} = letter;
            
        end
       
        if strcmp(sym,'chorus a')
            pause(0);
        end     
        
    end
    
    structs{1,i} = [symbs];
    structs{2,i} = [timestamps];
    structs{3,i} = [letters];
    structs{4,i} = strseq;
    fclose(fid);
    
end


save('data/structs','structs');
