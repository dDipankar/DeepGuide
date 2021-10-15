# DeepGuide
Genome-wide functional screens enable the prediction of high activity CRISPR-Cas9 and -Cas12a guides in _Yarrowia lipolytica_
# Installation
## Prerequisite
- python 3.6.10
- keras 2.4.3
- tensorflow 2.2.0
- scipy 1.4.1
- scikit-learn 0.23.1
- numpy 1.18.5
- pandas 1.0.5
## Installation Guide
The easiest way to get the needed prerequisites to run DeepGuide is through Anaconda. If you have Anaconda installed already you can skip this step, otherwise go to [https://docs.anaconda.com/anaconda/install/index.html](https://docs.anaconda.com/anaconda/install/index.html) to learn how to install conda on your system. We have used Anaconda version 4.8.5. Higher version is expected to work. Once Anaconda is correctly installed, You can run the following command to install requirements for DeepGuide

```
conda create -n deepguide python=3.6.10 ipykernel matplotlib pandas=1.0.5 numpy=1.18.5 scipy=1.4.1 tensorflow=2.2.0 keras=2.4.3 scikit-learn=0.23.1 biopython=1.71
```

## Running the software
### For cas12a
Assuming you have installed the prerequisites in a conda environment called deepguide, you can run the software for cas12a guides using following command

```
git clone https://github.com/dDipankar/DeepGuide
cd DeepGuide/src/
conda activate deepguide
python DetectAndScore_cas12a.py path_of_fasta_file
```
See data/seq_sample.fasta for FASTA format. Just remember that the program needs at least 32 nucleotides to fit the full target.

```
context (1nt) -- PAM (4nt, TTTV) -- target (25nt) -- context (2nt)
```

You will get an output file called activity_score_cas12a.csv in data directory. This file will contain the predicted cutting score by DeepGuide for each guide.

### For cas9(Sequence + Nucleosome Occupancy as Input)
You can run DeepGuide to get the prediction scores for cas9 guides using the sequence of guides and Nucleotide Occupancy by the following command

```
git clone https://github.com/dDipankar/DeepGuide
cd DeepGuide/src/
conda activate deepguide
python DetectAndScore_cas9.py path_of_fasta_file path_of_NucleosomeOccupancy_file
```
See data/seq_sample.fasta for FASTA format. Just remember that the program needs at least 28 nucleotides to fit the full target. 
See data/nu_sample.csv for nucleosome occupancy file. Here each number in the file represents nucleosome occupancy for each nucleotide potition of the fasta file.
Remember that total number of nuclesome occupanies has to be equal the total number of nucleotides in the fasta file


```
context (2nt) -- target (20nt) -- PAM (3nt, NGG) --  context (3nt)
```


### For cas9(only Sequence as Input)
To get prediction scores for cas9 guides using sequence only use the following command

```
git clone https://github.com/dDipankar/DeepGuide
cd DeepGuide/src/
conda activate deepguide
python DetectAndScore_cas9.py path_of_fasta_file
```
## Example run
### For cas12a
```
git clone https://github.com/dDipankar/DeepGuide
cd DeepGuide/src/
conda activate deepguide
python DetectAndScore_cas12a.py ../data/seq_sample.fasta
```
Then you will get an output file called activity_score_cas12a.csv in data directory. The format of the output is bellow:

| Guide  | Score |
| ------------- | ------------- |
| CTACCCGATATCTGTCACAGTCGTT  | 1.9962865114212036  |
| ACGACCCCAAGCTGACCGATGACTC  | 1.2889682054519653  |

### For cas9(Sequence + Nucleosome Occupancy as Input)

```
git clone https://github.com/dDipankar/DeepGuide
cd DeepGuide/src/
conda activate deepguide
python DetectAndScore_cas9.py  ../data/seq_sample.fasta  ../data/nu_sample.csv
```
Then you will get an output file called activity_score_cas9.csv in data directory. The format of the output is same as before.
| Guide  | Score |
| ------------- | ------------- |
| TCAAACGATTACCCACCCTC  | 3.4991283416748047  |
| TTACCCACCCTCCGGGACTG  | 4.1370463371276855  |

### For cas9(only Sequence as Input)
```
git clone https://github.com/dDipankar/DeepGuide
cd DeepGuide/src/
conda activate deepguide
python DetectAndScore_cas9_seq.py  ../data/seq_sample.fasta
```
Then you will get an output file called activity_score_cas9_seq.csv in data directory. The format of the output is same as above
## Citation
If you have used this tool in your publication please cite this

**Genome-wide functional screens enable the prediction of high activity CRISPR-Cas9 and -Cas12a guides in *Yarrowia lipolytica***. Dipankar Baisya, Adithya Ramesh, Cory Schwartz, Stefano Lonardi, and Ian Wheeldon. Submitted
