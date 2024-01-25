function projMats

% cat is a hard-coded list of all possible structs
cat = {'(secondary) theme', ... % 1 
    'bridge', ... % 2
    'chorus', ... % 3
    'chorus a', ... % 4
    'chorus b', ... % 5
    'coda', ... % 6
    'end', ... % 7
    'ending', ... % 8
    'fade in', ... % 9
    'fadein', ... % 10
    'fadeout', ... % 11
    'flute)', ... % 12
    'instrumental', ... % 13
    'instrumental break', ... % 14
    'interlude', ... % 15
    'intro', ... % 16
    'intro a', ... % 17
    'intro b', ... % 18
    'intro-a', ... % 19
    'intro-b', ... % 20
    'key change', ... % 21
    'main theme', ... % 22
    'modulation', ... % 23
    'outro', ... % 24
    'pre chorus', ... % 25
    'pre-chorus', ... % 26
    'pre-chorus two', ... % 27
    'pre-intro', ... % 28
    'pre-verse', ... % 29
    'prechorus', ... % 30
    'refrain', ... % 31
    'secondary theme', ... % 32
    'silence', ... % 33
    'solo', ... % 34
    'spoken', ... % 35
    'spoken verse', ... % 36
    'theme', ... % 37
    'trans', ... % 38
    'transition', ... % 39
    'verse', ... % 40
    'verse five', ... % 41
    'verse four', ... % 42
    'verse one', ... % 43 
    'verse three', ... % 44
    'verse two', ... % 45
    'vocal'}; % 46

% proj_Mat_A retains only the structs that comprise 99% of our data; the
% rest are collapsed into "other"
projMat_A = zeros(46,18);
projMat_A(1,18) = 1;
projMat_A(2,1) = 1;
projMat_A(3,2) = 1;
projMat_A(4,18) = 1;
projMat_A(5,18) = 1;
projMat_A(6,3) = 1;
projMat_A(7,4) = 1;
projMat_A(8,18) = 1;
projMat_A(9,18) = 1;
projMat_A(10,18) = 1;
projMat_A(11,5) = 1;
projMat_A(12,18) = 1;
projMat_A(13,6) = 1;
projMat_A(14,18) = 1;
projMat_A(15,7) = 1;
projMat_A(16,8) = 1;
projMat_A(17,18) = 1;
projMat_A(18,18) = 1;
projMat_A(19,18) = 1;
projMat_A(20,18) = 1;
projMat_A(21,18) = 1;
projMat_A(22,18) = 1;
projMat_A(23,18) = 1;
projMat_A(24,9) = 1;
projMat_A(25,10) = 1;
projMat_A(26,11) = 1;
projMat_A(27,18) = 1;
projMat_A(28,18) = 1;
projMat_A(29,12) = 1;
projMat_A(30,18) = 1;
projMat_A(31,18) = 1;
projMat_A(32,18) = 1;
projMat_A(33,13) = 1;
projMat_A(34,14) = 1;
projMat_A(35,18) = 1;
projMat_A(36,18) = 1;
projMat_A(37,18) = 1;
projMat_A(38,15) = 1;
projMat_A(39,16) = 1;
projMat_A(40,17) = 1;
projMat_A(41,18) = 1;
projMat_A(42,18) = 1;
projMat_A(43,18) = 1;
projMat_A(44,18) = 1;
projMat_A(45,18) = 1;
projMat_A(46,18) = 1;
projMat_A = projMat_A';

for i = 1:46
    vec = zeros(46,1);
    vec(i) = 1;
    idx = projMat_A * vec;
    
    for j = 1:17
        if idx(j) == 1
            cat_A(j) = cat(i);
        end
    end      
end

cat_A(18) = {'other'};


% proj_Mat_B retains only the structs that are not collapsible into other
% structs
projMat_B = zeros(46,17);
projMat_B(1,14) = 1;
projMat_B(2,1) = 1;
projMat_B(3,2) = 1;
projMat_B(4,2) = 1;
projMat_B(5,2) = 1;
projMat_B(6,8) = 1;
projMat_B(7,3) = 1;
projMat_B(8,8) = 1;
projMat_B(9,6) = 1;
projMat_B(10,6) = 1;
projMat_B(11,8) = 1;
projMat_B(12,4) = 1;
projMat_B(13,4) = 1;
projMat_B(14,4) = 1;
projMat_B(15,5) = 1;
projMat_B(16,6) = 1;
projMat_B(17,6) = 1;
projMat_B(18,6) = 1;
projMat_B(19,6) = 1;
projMat_B(20,6) = 1;
projMat_B(21,7) = 1;
projMat_B(22,14) = 1;
projMat_B(23,7) = 1;
projMat_B(24,8) = 1;
projMat_B(25,9) = 1;
projMat_B(26,9) = 1;
projMat_B(27,9) = 1;
projMat_B(28,10) = 1;
projMat_B(29,11) = 1;
projMat_B(30,9) = 1;
projMat_B(31,2) = 1;
projMat_B(32,14) = 1;
projMat_B(33,12) = 1;
projMat_B(34,4) = 1;
projMat_B(35,13) = 1;
projMat_B(36,13) = 1;
projMat_B(37,14) = 1;
projMat_B(38,15) = 1;
projMat_B(39,15) = 1;
projMat_B(40,16) = 1;
projMat_B(41,16) = 1;
projMat_B(42,16) = 1;
projMat_B(43,16) = 1;
projMat_B(44,16) = 1;
projMat_B(45,16) = 1;
projMat_B(46,17) = 1;

projMat_B = projMat_B';

cat_B = {'bridge', ...
    'chorus', ...
    'end', ...
    'instrumental', ...
    'interlude', ...
    'intro', ...
    'key change', ...
    'outro', ...
    'pre-chorus', ...
    'pre-intro', ...
    'pre-verse', ...
    'silence', ...
    'spoken', ...
    'theme', ...
    'transition', ...
    'verse', ...
    'vocal'};

projMat_C = zeros(46,4);
projMat_C(1,4) = 1;
projMat_C(2,4) = 1;
projMat_C(3,1) = 1;
projMat_C(4,4) = 1;
projMat_C(5,4) = 1;
projMat_C(6,4) = 1;
projMat_C(7,4) = 1;
projMat_C(8,4) = 1;
projMat_C(9,4) = 1;
projMat_C(10,4) = 1;
projMat_C(11,4) = 1;
projMat_C(12,4) = 1;
projMat_C(13,4) = 1;
projMat_C(14,4) = 1;
projMat_C(15,4) = 1;
projMat_C(16,4) = 1;
projMat_C(17,4) = 1;
projMat_C(18,4) = 1;
projMat_C(19,4) = 1;
projMat_C(20,4) = 1;
projMat_C(21,4) = 1;
projMat_C(22,4) = 1;
projMat_C(23,4) = 1;
projMat_C(24,4) = 1;
projMat_C(25,4) = 1;
projMat_C(26,4) = 1;
projMat_C(27,4) = 1;
projMat_C(28,4) = 1;
projMat_C(29,4) = 1;
projMat_C(30,4) = 1;
projMat_C(31,4) = 1;
projMat_C(32,4) = 1;
projMat_C(33,2) = 1;
projMat_C(34,4) = 1;
projMat_C(35,4) = 1;
projMat_C(36,4) = 1;
projMat_C(37,4) = 1;
projMat_C(38,4) = 1;
projMat_C(39,4) = 1;
projMat_C(40,3) = 1;
projMat_C(41,4) = 1;
projMat_C(42,4) = 1;
projMat_C(43,4) = 1;
projMat_C(44,4) = 1;
projMat_C(45,4) = 1;
projMat_C(46,4) = 1;

projMat_C = projMat_C';

cat_C = {'chorus'
    'silence'
    'verse'
    'other'};

cat_C = cat_C';


projMat_D = zeros(46,5);
projMat_D(1,5) = 1;
projMat_D(2,5) = 1;
projMat_D(3,1) = 1;
projMat_D(4,5) = 1;
projMat_D(5,5) = 1;
projMat_D(6,5) = 1;
projMat_D(7,2) = 1;
projMat_D(8,5) = 1;
projMat_D(9,5) = 1;
projMat_D(10,5) = 1;
projMat_D(11,5) = 1;
projMat_D(12,5) = 1;
projMat_D(13,5) = 1;
projMat_D(14,5) = 1;
projMat_D(15,5) = 1;
projMat_D(16,5) = 1;
projMat_D(17,5) = 1;
projMat_D(18,5) = 1;
projMat_D(19,5) = 1;
projMat_D(20,5) = 1;
projMat_D(21,5) = 1;
projMat_D(22,5) = 1;
projMat_D(23,5) = 1;
projMat_D(24,5) = 1;
projMat_D(25,5) = 1;
projMat_D(26,5) = 1;
projMat_D(27,5) = 1;
projMat_D(28,5) = 1;
projMat_D(29,5) = 1;
projMat_D(30,5) = 1;
projMat_D(31,5) = 1;
projMat_D(32,5) = 1;
projMat_D(33,3) = 1;
projMat_D(34,5) = 1;
projMat_D(35,5) = 1;
projMat_D(36,5) = 1;
projMat_D(37,5) = 1;
projMat_D(38,5) = 1;
projMat_D(39,5) = 1;
projMat_D(40,4) = 1;
projMat_D(41,5) = 1;
projMat_D(42,5) = 1;
projMat_D(43,5) = 1;
projMat_D(44,5) = 1;
projMat_D(45,5) = 1;
projMat_D(46,5) = 1;

projMat_D = projMat_D';

cat_D = {'chorus'
    'silence'
    'verse'
    'other'};

cat_D = cat_D';


save('data/projections','cat','projMat_A','cat_A','projMat_B', ...
    'cat_B','projMat_C','cat_C','projMat_D','cat_D');