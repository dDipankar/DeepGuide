%Picks sgRNA from both Top and Bottom strands in an unbiased manner (equal
%representation as much as possible), to make an 8-fold coverage sgRNA
%library (8 sgRNA/gene when possible).

%The way this is done is by:
% (i) Picking 4 sgRNA each from Top and Bottom strands if each of them
% individually have more than 4 sgRNA.
% (ii) If only one of either the Top or Bottom strands have less than 4 sgRNA,
% then all of those are picked, and the remaining sgRNA are picked from the
% other strand to make a total of 8.
% (iii) If both Top and bottom strands have less than 4 sgRNA each, then
% all of the sgRNA from both strands are picked.

% Input Variables:
% For_Cas12a_unique2
% Rev_Cas12a_unique2

% Important Output Variables:
% Cas12a_Library_final: Final 8-fold coverage Cas12a sgRNA library for
% Y. lipolytica CLIB89.


%Author: Adithya Ramesh
%PhD Candidate, Wheeldon Lab
%UC Riverside, 900 University Ave
%Riverside, CA-92507, USA
%Email: arame003@ucr.edu
%       wheeldon@ucr.edu
%% Design of unbiased 8-fold coverage sgRNA library

isCas12for = ~cellfun('isempty',For_Cas12a_unique2);
isCas12rev = ~cellfun('isempty',Rev_Cas12a_unique2);
csCas12for = cumsum(isCas12for,2);
csCas12rev = cumsum(isCas12rev,2);
maxCas12for = max(csCas12for,[],2);
maxCas12rev = max(csCas12rev,[],2);
or4 = xor(maxCas12for<4,maxCas12rev<4); %True when either condition is true and false when both conditions are true or false
R = numel(or4);
Cas12a_Library_final = cell(R,8);

for k = 1:R
    if maxCas12for(k)<4
        % bpfor < 4 elements
        V = [For_Cas12a_unique2(k,isCas12for(k,:)),...
             Rev_Cas12a_unique2(k,isCas12rev(k,:))];
    elseif maxCas12rev(k)<4
        % bprev < 4 elements
        V = [Rev_Cas12a_unique2(k,isCas12rev(k,:)),...
             For_Cas12a_unique2(k,isCas12for(k,:))];
    else     
        % "bpfor and bprev individually have more than 4 non-empty elements"
        % and also "both bpfor and bprev individually have less than 4 elements"
        V = [For_Cas12a_unique2(k,isCas12for(k,:)&csCas12for(k,:)<=4),...
             Rev_Cas12a_unique2(k,isCas12rev(k,:)&csCas12rev(k,:)<=4)];
    end
    idx = 1:min(size(Cas12a_Library_final,2),numel(V));
    Cas12a_Library_final(k,idx) = V(idx);
end