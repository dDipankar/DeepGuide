% This script generates sgRNA from the top and bottome strands of both all
% the CDS, as well as from the entire genome of CLIB89.

% Inputs: Variables 'CDS' and 'Chromosome' which are outputs of the sript 
% Gene_Import_Cas12a.m

% Important Output variables:
% For_Cas12a: List of all Top strand sgRNA from the given list of CDS.
% Rev_Cas12a: List of all Bottom strand sgRNA from the given list of CDS.
% For_Cas12a_context: Each sgRNA in For_Cas12a is listed along with its PAM
% and a couple of nucleotides (nt) on each side to total length of 32 nt.
% Rev_Cas12a_context: Same as above for Rev_Cas12a.
% All_Cas12a: List of all sgRNA (top and bottom) for each full length
% chromosome in the genome.


%Author: Adithya Ramesh
%PhD Candidate, Wheeldon Lab
%UC Riverside, 900 University Ave
%Riverside, CA-92507, USA
%Email: arame003@ucr.edu
%       wheeldon@ucr.edu
%% Generation of all Top and Bottom strand sgRNA in each CDS
tic
for n=1:length(CDS)
    ct_f=0;
    ct_r=0;
    
    %Top strand
    temp_f=CDS{n,1};
    for i=1:length(temp_f)-29        
        if(((temp_f(i)=='T')&&(temp_f(i+1)=='T')&&(temp_f(i+2)=='T')&&(temp_f(i+3)=='A'))||((temp_f(i)=='T')&&(temp_f(i+1)=='T')&&(temp_f(i+2)=='T')&&(temp_f(i+3)=='G'))||((temp_f(i)=='T')&&(temp_f(i+1)=='T')&&(temp_f(i+2)=='T')&&(temp_f(i+3)=='C')))
            j=i+4;
            k=j+24;
            x=i-1;
            y=k+2;
            if(i<=length(temp_f))
                ct_f=ct_f+1;                
                For_Cas12a{n,ct_f}=temp_f(j:k);
            end
            if (y<=length(temp_f) && x>0)
                For_Cas12a_context{n,ct_f}=temp_f(x:y);
            end
        end        
    end    
    
    %Bottom Strand
    temp_r=seqcomplement(temp_f);
    for i=26:length(temp_r)-4        
        if(((temp_r(i)=='A')&&(temp_r(i+1)=='T')&&(temp_r(i+2)=='T')&&(temp_r(i+3)=='T'))||((temp_r(i)=='C')&&(temp_r(i+1)=='T')&&(temp_r(i+2)=='T')&&(temp_r(i+3)=='T'))||((temp_r(i)=='G')&&(temp_r(i+1)=='T')&&(temp_r(i+2)=='T')&&(temp_r(i+3)=='T')))
            j=i-1;
            k=j-24;
            x=i+1;
            y=k-2;
            if(i+3<=length(temp_r))
                ct_r=ct_r+1;                
                Rev_Cas12a{n,ct_r}=reverse(temp_r(k:j)); 
            end
            if (x<=length(temp_f) && y>0)
                Rev_Cas12a_context{n,ct_r}=reverse(temp_r(y:x));
            end
        end        
    end    

end
toc
%% Generation of all sgRNA in the entire genome
tic
for n=1:length(Chromosome)
    ct=0;    
    
    %Top Strand
    temp_f=cell2mat(Chromosome{n,1});
    for i=1:length(temp_f)-29        
        if(((temp_f(i)=='T')&&(temp_f(i+1)=='T')&&(temp_f(i+2)=='T')&&(temp_f(i+3)=='A'))||((temp_f(i)=='T')&&(temp_f(i+1)=='T')&&(temp_f(i+2)=='T')&&(temp_f(i+3)=='G'))||((temp_f(i)=='T')&&(temp_f(i+1)=='T')&&(temp_f(i+2)=='T')&&(temp_f(i+3)=='C'))||((temp_f(i)=='T')&&(temp_f(i+1)=='T')&&(temp_f(i+2)=='T')&&(temp_f(i+3)=='T')))
            j=i+4;
            k=j+24;            
            if(i<=length(temp_f))
                ct=ct+1;                
                All_Cas12a{n,ct}=temp_f(j:k);
            end
        end
    end
    
    %Bottom Strand
    temp_r=seqcomplement(temp_f);
    for i=26:length(temp_r)-4        
        if(((temp_r(i)=='A')&&(temp_r(i+1)=='T')&&(temp_r(i+2)=='T')&&(temp_r(i+3)=='T'))||((temp_r(i)=='C')&&(temp_r(i+1)=='T')&&(temp_r(i+2)=='T')&&(temp_r(i+3)=='T'))||((temp_r(i)=='G')&&(temp_r(i+1)=='T')&&(temp_r(i+2)=='T')&&(temp_r(i+3)=='T'))||((temp_r(i)=='T')&&(temp_r(i+1)=='T')&&(temp_r(i+2)=='T')&&(temp_r(i+3)=='T')))
            j=i-1;
            k=j-24;            
            if(i+3<=length(temp_r))
                ct=ct+1;                
                All_Cas12a{n,ct}=reverse(temp_r(k:j)); 
            end          
        end
    end 

end
toc
