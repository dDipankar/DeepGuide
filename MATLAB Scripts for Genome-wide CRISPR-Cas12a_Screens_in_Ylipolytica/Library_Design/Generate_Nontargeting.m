% This script generates nontargeting sgRNA that function as negative
% controls in a genome-wide CRISPR screen. These sgRNAs ideally will not
% target anywhere within the genome.

% This is done by first breaking the genome (both top and bottom strands)
% into 10 nt chunks (base pairs 1-10, 2-11,3-12...etc). Then a random 25 nt
% 'sgRNA' is designed, and then first 10 nt of this sgRNA is compared to
% all the 10 nt chunks within the genome. If there is no match, then that
% sgRNA is considered a nontargeting guide. This is done for a total of
% 1,000,000 such randomly sequences, and any that pass the criteria for a
% nontargeting sgRNA are kept. The number that pass the check for a
% nontargeting sgRNA will be much smaller than 1,000,000. (Note: Should
% more nontargeting sgRNA be required, the number of random sequences
% generated can be increased from 1,000,000 to any desired number.)

% Inputs: Full Chromosome files of the genome (ending with _full.fasta)
% Important Outputs: List of nontargeting sgRNA (variable: Nontargeting)

% Note: This script may take a while to run and thus the parfor loops have
% been utilized to speed up runtime. This requires that the MATLAB version
% also has the Parallel Computing Toolbox installed. This script may also be
% run on the cluster to speed it up.


%Author: Adithya Ramesh
%PhD Candidate, Wheeldon Lab
%UC Riverside, 900 University Ave
%Riverside, CA-92507, USA
%Email: arame003@ucr.edu
%       wheeldon@ucr.edu
%% Importing each Chromosome from its fasta file into a string array variable

d=65;
for c=1:6
    [Chromofilename, Chromopathname] = uigetfile('.fasta', 'pick files ending with "_full.fasta"');
    Chromopathfile = [Chromopathname, Chromofilename];
    [~, Csequence] = fastaread(Chromopathfile);
    
    %Read the Chromosome sequences
    if ischar(Csequence)==1
        Cseq = Csequence;
        Csequence = cell(1);
        Csequence{1} = Cseq;
    end
    
    Cheader = char(d);
    d = d+1;
    
    Chromosome{c,1} = Csequence;
    Chromosome{c,2} = Cheader;
end
%% Generate list of all 10 bp substrings of the entire genome (Top and Bottom strands)
for i=1:length(Chromosome) 
    Seq=cell2mat(Chromosome{i,1});
    SeqComp=seqcomplement(Seq);
    for j=1:length(Seq)-10
        Allseqs{i,j}=Seq(j:j+9);        
    end
    c=j;
    for j=1:length(SeqComp)-10
        Allseqs{i,c}=SeqComp(j:j+9);
        c=c+1;
    end
end
%% Linearize all 10 bp substrings intoa variable with one column

idx=sum(~cellfun(@isempty,Allseqs),2);
c=1;

for i=1:6
    Allseqs_Linear(c:c+idx(i)-1,1)=Allseqs(i,1:idx(i))';
    c=c+idx(i);
end

S=unique(Allseqs_Linear); % Only the unique 10 bp substrings are kept
%% Random DNA sequence generator
%Varabile 'tot' decides total number of random sequences generated. This
%number can be increased to increase the number of nontargeting sgRNA that
%are finally generated.

tot=100000;
Rand=cell(tot,1);

for i=1:tot
    Rand{i,1}=randseq(25);
end
%% Check for uniqueness of random string

check=ones(tot,1);

parfor i=1:length(Rand)
    a=cell2mat(Rand(i,1));
    z=strfind(S,a);
    w=find(~cellfun(@isempty,z));
    tot(i,1)=numel(w);
end

idx=find(tot==0);
Nontargeting=Rand(idx);