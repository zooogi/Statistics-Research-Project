# Bayesian Survival Project

This project accompanies my MSc dissertation on **model checking in Bayesian survival analysis**, using an employee turnover dataset as a motivating example. The work explores the role of the observation window as an estimable parameter, evaluates its impact on censoring and inference, and applies posterior predictive checks to assess model adequacy.

## Project Structure
```bash
bayesian-survival-project/
│
├── data/                  
│   ├── turnover.csv
│
├── scr/                 
│   └── ...
│
├── images/               # All graphs generated after running R
|
├── paper/                 # Overleaf（LaTeX）
│   ├── main.tex
│   ├── sections/
│   │   ├── introduction.tex
│   │   ├── methods.tex
│   │   ├── results.tex
│   │   ├── discussion.tex
│   |   └── ...
|
├── README.md             
└── .gitignore            
```
## How to Reproduce
1. Run scripts in scr/ to:

   -Process the dataset in data/turnover.csv

   -Fit Bayesian survival models

   -Generate posterior predictive checks

   -Export figures to images/

3. Compile the LaTeX files in paper/ (or open on Overleaf) to generate the thesis.

## Dataset
The employee turnover dataset is publicly available at https://www.kaggle.com/datasets/davinwijaya/employee-turnover/data

## Reproducibility
All code and instructions are provided to ensure full reproducibility. Running the scripts will regenerate every table and figure reported in the paper.
