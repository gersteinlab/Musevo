function A5_get_gnrs_data

% get genres data

fin = 'data/parsed_songs_log_0_3_0_1_100';

fid = fopen(fin);
ln = fgetl(fid);
ids = {};
dts = {};
sqs = {};
vals = [];
gnrs = {};

sp = find(ln==sprintf('\t'));
% sqs = strsplit(ln(sp(4)+1:end));
sqs = strsplit(ln(sp(5)+1:end));
nsq = length(sqs);

i = 1;
while(1)
    ln = fgetl(fid);
    if ln==-1
        break;
    end
    sp = find(ln==sprintf('\t'));
    ids{i} = ln(1:sp(1)-1);
    dts{i} = ln(sp(3)+1:sp(4)-1);
    gnrs{i} = ln(sp(4)+1:sp(5)-1);
    if dts{i}=='0'
        val = zeros(1,size(vals,2));
    else
        val = ln(sp(5)+1:end);
        val = str2num(val);
    end
    vals = [vals ; val];
    i = i + 1;
end

fclose(fid);

validx = find(~strcmp(dts,'0'));
ids = ids(validx);
dts = dts(validx);
vals = vals(validx,:);
% 
% save('data/ids_dedup','ids');
% save('data/gnrs','gnrs');

