% This script runs a test for uniqueness for every sgRNA in For_Cas12a and
% Rev_Cas12a. First 14 nt of each sgRNA for each CDS (seed) is compared to
% the first 14 nt of every sgRNA in the genome. If there is more than one
% match or zero matches, the sgRNA is dropped as non-unique or not present
% respectively.

% Input variables (from previous scripts):
% For_Cas12a
% Rev_Cas12a
% All_Cas12a
% CDS
% Chromosome

% Important Output Variables:
% For_Cas12a_unique: All unique sgRNA in For_Cas12a.
% Rev_Cas12a_unique: All unique sgRNA in Rev_Cas12a.
% For_Cas12a_unique2: Compressed list of unique For_Cas12a sgRNA to remove
% empty spaces due to deletion of non-unique sgRNA.
% Rev_Cas12a_unique2: Compressed list of unique Rev_Cas12a sgRNA to remove
% empty spaces due to deletion of non-unique sgRNA.

% Note: This script will take a while to run and thus the parfor loops have
% been utilized to save time. This requires that the MATLAB version
% also has the Parallel Computing toolbox installed. This script may also be
% run on the cluster to speed it up.


%Author: Adithya Ramesh
%PhD Candidate, Wheeldon Lab
%UC Riverside, 900 University Ave
%Riverside, CA-92507, USA
%Email: arame003@ucr.edu
%       wheeldon@ucr.edu
%% Initialization
For_Cas12a_unique=For_Cas12a;
Rev_Cas12a_unique=Rev_Cas12a;

[rowf,colf]=size(For_Cas12a);
[rowr,colr]=size(Rev_Cas12a);

N=200; %A number greater than maximum number of columns in For_Cas12a and Rev_cas12a.

nCDS=length(CDS); %number of CDS
nChr=6; %number of Chromosomes in CLIB89
lChr=length(All_Cas12a); %max number of sgRNA in any chromosome

%% Test for uniqueness in Top strand
% The test for uniqueness here compares the first 14 nt of each sgRNA in
% For_Cas12a to every sgRNA in All_Cas12a. If there is more than
% one match or zero matches, the sgRNA is dropped as non-unique or not
% present respectively.

parfor i=1:nCDS
    for j=1:colf
        if ~isempty(For_Cas12a{i,j})
            flag=0;
            for k=1:nChr
                for l=1:lChr
                    temp1=For_Cas12a{i,j};
                    if ~isempty(All_Cas12a{k,l})
                        temp2=All_Cas12a{k,l};
               
                        if(strcmp(temp1(1:14),temp2(1:14)))
                            flag=flag+1;
                        end
                        if flag==2
                            break;
                        end       
                    end
                end
            end
            if(flag>1||flag==0)
                For_Cas12a_unique{i,j}={};
            end
        end
    end
end

%% Test for uniqueness in Bottom strand
% The test for uniqueness here compares the first 14 nt of each sgRNA in
% Rev_Cas12a to every sgRNA in All_Cas12a. If there is more than
% one match or zero matches, the sgRNA is dropped as non-unique or not
% present respectively.

parfor i=1:nCDS
    for j=1:colr
        if ~isempty(Rev_Cas12a{i,j})
            flag=0;
            for k=1:nChr
                for l=1:lChr
                    temp1=Rev_Cas12a{i,j};
                    if ~isempty(All_Cas12a{k,l})
                        temp2=All_Cas12a{k,l};
               
                        if(strcmp(temp1(1:14),temp2(1:14)))
                            flag=flag+1;
                        end
                        if flag==2
                            break;
                        end       
                    end
                end
            end
            if(flag>1||flag==0)
                Rev_Cas12a_unique{i,j}={};
            end
        end
    end
end

%% Compress list of unique sgRNA to remove empty spaces due to deletion of non-unique sgRNA

For_Cas12a_unique2=For_Cas12a_unique;
Rev_Cas12a_unique2=Rev_Cas12a_unique;
sz1=size(For_Cas12a_unique2);
sz2=size(Rev_Cas12a_unique2);

for k=1:500
    for i=1:7919
        for j=1:(sz1(1,2)-1)
            if(isequal(isempty(For_Cas12a_unique2{i,j}),1))
               For_Cas12a_unique2{i,j}=For_Cas12a_unique2{i,j+1};
               For_Cas12a_unique2{i,j+1}={};
            end            
        end
    end
end

for k=1:500
    for i=1:7919
        for j=1:(sz2(1,2)-1)            
            if (isequal(isempty(Rev_Cas12a_unique2{i,j}),1))
                Rev_Cas12a_unique2{i,j}=Rev_Cas12a_unique2{i,j+1};
                Rev_Cas12a_unique2{i,j+1}={};
            end
        end
    end
end
