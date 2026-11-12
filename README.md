# ACE Docking Comprehensive Analyzer

A robust Python utility to perform deep statistical evaluation, frequency analysis, correlation testing (Pearson & Chi-Square), and automated publication-quality figure generation (SEM error bars, correlation regression, interaction pie charts) for ACE-peptide docking screening results.

## Features

- Computes frequency and binomial Standard Error of the Mean (SEM) for receptor active-site residues.
- Performs Pearson correlation between docking scores and total interaction counts.
- Executes Chi-Square goodness-of-fit tests on interaction distributions.
- Ranks candidate peptides using a multi-criteria scoring system.
- Exports high-resolution charts (`300 DPI`) and a multi-sheet comprehensive Excel report.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

Run in the current directory:
```bash
python comprehensive_analysis.py
```

Or specify custom paths and filenames:
```bash
python comprehensive_analysis.py /path/to/root -i Selected_Models_Interactions_Cleaned.xlsx -o Custom_Analysis_Report.xlsx
```
