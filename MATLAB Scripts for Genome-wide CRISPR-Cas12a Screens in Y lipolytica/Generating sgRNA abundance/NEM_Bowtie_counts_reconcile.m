% A simple script that switches the Bowtie counts of an sgRNA in any
% sample, with the Naive Exact Matching counts provided the Bowtie counts
% are lesser than the exact maathcing counts for that sgRNA in that sample.

% Note: This is important because while Bowtie allows for inexact matching
% (accounting for sequencing/PCR errors etc.) and therefore generates more
% accurate counts in a majority of the sgRNA, Bowtie was also designed for
% mRNA-seq that have longer read lengths than the sgRNAs. This leads to
% errors where some sgRNAs are imporperly aligned and thus have inaccurate
% counts. Thus, this ensures that the lowest value of counts for any sgRNA
% is the NEM value.

% Input variables:
% 1. NEM_counts_1: The NEM counts for a given sample.
% 2. Bowtie_counts_1: Bowtie counts for the same sample.

% Output variables:
% Final_counts_1.mat: Variable file that contains the final counts for any
% given sample.


%Author: Adithya Ramesh
%PhD Candidate, Wheeldon Lab
%UC Riverside, 900 University Ave
%Riverside, CA-92507, USA
%Email: arame003@ucr.edu
%       wheeldon@ucr.edu
%% Switch Bowtie Counts with NEM counts

tic
for i=1:length(NEM_counts_1)
    if Bowtie_counts_1(i,1)<NEM_counts_1(i,1)
       Bowtie_counts_1(i,1)=NEM_counts_1(i,1);
    end
end
Final_counts_1=Bowtie_counts_1;
save Final_counts_1.mat Final_counts_1
toc