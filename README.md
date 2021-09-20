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

## Runing the software
### For cas12a
Assuming you have installed the prerequisites in a conda environment called deepguide, you can run the software using for cas12a guides using following command

```
conda activate deepguide
python DetectAndScore_cas12a.py path_of_fasta_file
```
See data/seq_fa.fasta for FASTA format. Just remember that the program needs at least 32 nucleotides to fit the full target.

```
context (1nt) -- PAM (4nt, TTTV) -- target (25nt) -- context (2nt)
```

You will get an output file called activity_score_cas12a.csv in data directory. This file will contain the predicted cutting score by DeepGuide for each guide.

## Example run
### For cas12a
```
conda activate deepguide
python DetectAndScore_cas12a.py ../data/seq_fa.fasta
```
Then you will get an output file called activity_score_cas12a.csv in data directory.
