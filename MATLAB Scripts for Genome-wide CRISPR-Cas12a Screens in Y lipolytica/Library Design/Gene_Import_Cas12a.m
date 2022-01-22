% This script imports all the names and sequences of the protein coding 
% CDS in the Y. lipolytica CLIB89 strain into a MATLAB variable. The full
% sequence of all 6 chromosomes in CLIB89 are also imported into
% another variable.

% Inputs: CDS files from each of the 6 chromosomes in CLIB89, as well as the 
% 6 full length chromosome sequences, that are provided in the same folder 
% as this script file.

% Ouputs: The variable 'CDS' gives the sequence of every CDS in the first
% column and its associated name/header in the corresponding cell of the 
% second column. The variable 'Chromosome' gives the sequence of each
% chromosome in the first column and its name in the second.


%Author: Adithya Ramesh
%PhD Candidate, Wheeldon Lab
%UC Riverside, 900 University Ave
%Riverside, CA-92507, USA
%Email: arame003@ucr.edu
%       wheeldon@ucr.edu
%% Importing each CDS from its fasta file into a string array variable

x=input('Input number of chromosomes for CDS to be imported from (CLIB89 has 6 chromosomes): ');

for a=1:x
    [myfilename, mypathname] = uigetfile('.fasta','Pick a file (ending with "_CDS.fasta"');
    myPathFile = [mypathname,myfilename];
    [header,sequence] =  fastaread(myPathFile);
    if ischar(sequence)==1
        seq=sequence;
        sequence=cell(1);
        sequence{1}=seq;
    end
    
    %Read the CDS sequences
    gene_cell=cell(length(sequence),1);
    for b=1:length(sequence)
        gene_cell{b}=sequence{1,b};
    end
    genes{a,1}=gene_cell;
    
    %Place CDS headers(names) into the cell
    name=cell(length(header),1);
    for b=1:length(header)
       if b<10
            name{b}=header{b}(44:56);
       elseif b>=10&&b<100
            name{b}=header{b}(45:57);
       elseif b>=100&&b<1000
            name{b}=header{b}(46:58);
       else
            name{b}=header{b}(47:59);
       end
    end
    genes{a,2}=name;
end
CDS = [genes{1,:};genes{2,:};genes{3,:};genes{4,:};genes{5,:};genes{6,:}];

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
