file1 = 'Test.fas'; %name of the .fa file that contains sequences from 1 timespot
file2 = 'seq2.fas'; %name of the .fa file that contains sequences from 2 timespot
file3 = 'seq3.fas'; %name of the .fa file that contains sequences from 3 timespot
b1 = 1; %number of the first sequence to be read
b2 = 21; %number of the last sequnce to be read
am = (b2-b1+1); %amount of the sequences
tr1 = 0.05; % deletion threshold
fcut = 0.05; %monomorphous threshhold
gen = {}; %container for raw info
data = {}; %container for result info
%take data from files
gen{1} = fastaread(file1,BlockRead=[b1 b2],IgnoreGaps=false);
%gen{2} = fastaread(file2,BlockRead=[b1 b2],IgnoreGaps=false);
%gen{3} = fastaread(file3,BlockRead=[b1 b2],IgnoreGaps=false); BinData.mat
 l = length(gen{1}(1).Sequence);
for vsch = 2:am
 if l>length(gen{1}(vsch).Sequence)
l = length(gen{1}(vsch).Sequence);
 end
end
%l - length of the shortest sequence


for k = 1:1
    %data cleaning
    banlistSeq = [];
banlistSit = [];
erram = zeros(am,1);
potential = [];

for i0 = 1:l
    del = 0;
    sch = 1;
for j0 = 1:am
if isempty(banlistSeq) || j0 ~= banlistSeq
seq = gen{k}(j0).Sequence;

if seq(i0) =='-' ||  seq(i0) =='N'
potential = [potential j0];
    del = del+1;
end
else
    sch=sch+1;
end
end

if del ~=0
    
if del/(am-sum(erram)) >=tr1
  
banlistSit = [banlistSit i0];
potential = [];
else
    
erram(potential,1) =  erram(potential,1)+1;
banlistSeq = [banlistSeq potential];
banlistSeq = sort(banlistSeq);
potential = [];
end
end

end

if ~isempty(banlistSeq)

gen{k}(banlistSeq) = [];
end
[am,~] = size(gen{k});

if ~isempty(banlistSit)
    for sch3 = 1:am
  gen{k}(sch3).Sequence(banlistSit) = [];  
    end
end
l = length(gen{k}(1).Sequence);

data{k} = zeros(am,l);

%finding the common consensus 
cons = '';
for i = 1:l
    qA=0;
    qC=0;
    qT=0;
    qG=0;
    for j = 1:am
        seq = gen{k}(j).Sequence;
        switch seq(i)
            case 'A'
qA = qA+1;
            case 'C'
qC = qC+1;
            case 'T'
qT = qT+1;
            case 'G'
qG = qG+1;    
        end
    

    end

    [~,I] = max([qA qC qT qG]);
switch I
    case 1
cons = [cons 'A'];
    case 2
cons = [cons 'C'];        
    case 3
cons = [cons 'T'];        
    case 4
cons = [cons 'G']; 
end
end

%binarization
for i2 = 1:l
    for j2 = 1:am
        seq = gen{k}(j2).Sequence;         
        if seq(i2) == cons(i2)
            data{k}(j2,i2)=0;
        else
            data{k}(j2,i2)=1;
        end
    end
end

%Removing monomorph
mono = [];
for msch = 1:l
    if mean(data{k}(:,msch))>=1-fcut || mean(data{k}(:,msch))<=fcut
        mono = [mono msch];
    end
end
data{k}(:,mono) = [];
end
save(strcat('BinData','.mat'),'data');
disp(data);
