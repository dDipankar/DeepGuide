% This script imports a fastq file and stores it into the variables 'header'
% and 'sequence'. The DNA sequences/reads stored in the fastq file are then
% tabulated to provide every unique read and the associated count or the
% number of times the read appears in the fastq file.

% Inputs: fastq file to import (here portrayed by
% 'Examplefile1.fastqsanger') *Note*: This source of this fastq file is
% likley the output of using Trimmomatic on the Illumina raw reads of a
% given sample to keep only the 25 nt sgRNA sequence. For more information
% regarding processing Illumina reads on Galaxy refer " ", and specifically
% Supplementary table " ".

% Outputs Variables:s
% UniqueSeqs: Lists every unique sequence/read in the fastq file.
% Counts: Provides number of times each of those unique read appeared in fastq file.
% Counts_of_all_unique_reads1.mat: saves the above two variables locally as
% a MATLAB variable file.

% *Note*: This script requires the use the fast and efficient custom
% 'count_unique' function written by the user Anthony Kendall. This
% function may be downlaoded at:
% https://www.mathworks.com/matlabcentral/fileexchange/23333-determine-and-count-unique-values-of-an-array


%Author: Adithya Ramesh
%PhD Candidate, Wheeldon Lab
%UC Riverside, 900 University Ave
%Riverside, CA-92507, USA
%Email: arame003@ucr.edu
%       wheeldon@ucr.edu
%% Import fastq files, and generate the counts of every unique sequence
tic
clear
for a=1:1
    [header,sequence] =  fastqread('Examplefile1.fastqsanger'); %imports header and sequence of a fastq reads file of a given sample
    if ischar(sequence)==1
        seq=sequence;
        sequence=cell(1);
        sequence{1}=seq;
    end
end
Seq=sequence';
[UniqueSeqs1,Counts1]=count_unique(Seq); %Custom function that goes through the list of all sequences and provides the counts of every unique sequence in the list
save Counts_of_all_unique_reads1.mat UniqueSeqs1 Counts1
toc