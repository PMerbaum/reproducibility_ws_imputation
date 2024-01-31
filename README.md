# repeat_imputation

This project is a part of short tandem repeat (STR) imputation accuracy analysis.
Allele lengths are cruisial in some cases. For example, for patients with amyotrophic lateral sclerosis - a C9ortf72 gene that has longer than 30 STR's means shorter survival and more aggressive disease. Therefore, it is crucial to be able to predict (impute) the length of C9ORF72 gene in all samples. We established a method to calculate these lengths for SNP array data.
We tested our model with known repeat lengths calculated with ExpansionHunter - we masked C9ortf72 gene region, reimputed it and now we check, how accurate our imputation is with build_roc_imputation_accuracy.py. 

## Usage

Upload the build_roc_imputation_accuracy.py function and simulated data from data/raw folder. Data is formated as a dataframe and contains name of samples, phased imputed lengths that were extracted from imputed.vcf file, phased true lengths extracted from ExpansionHunter.vcf file and a dosage calulated during the imputation with Beagle software.

## Project Structure

The project structure distinguishes three kinds of folders:
- read-only (RO): not edited by either code or researcher
- human-writeable (HW): edited by the researcher only.
- project-generated (PG): folders generated when running the code; these folders can be deleted or emptied and will be completely reconstituted as the project is run.


```
.
├── .gitignore
├── LICENSE
├── README.md
├── requirements.txt
├── data               <- All project data, ignored by git
│   ├── raw            <- The original, immutable data dump. (RO)
├── docs               <- Documentation notebook for users (HW)
│   ├── manuscript     <- Manuscript source, e.g., LaTeX, Markdown, etc. (HW)
│   └── reports        <- Other project reports and notebooks (e.g. Jupyter, .Rmd) (HW)
├── results
│   ├── figures        <- Figures for the manuscript or reports (PG)
│   └── output         <- Other output for the manuscript or reports (PG)
└── src                <- Source code for this project (HW)

```

## Add a citation file
Create a citation file for your repository using [cffinit](https://citation-file-format.github.io/cff-initializer-javascript/#/)

## License

This project is licensed under the terms of the [MIT License](/LICENSE).
