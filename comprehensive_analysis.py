# -*- coding: utf-8 -*-
import argparse
from pathlib import Path
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
from scipy import stats


def run_comprehensive_analysis(root_dir: Path, input_excel_name: str, output_excel_name: str):
    """Run comprehensive analysis on cleaned ACE docking interactions, calculate statistics, and generate plots."""
    root_dir = root_dir.resolve()
    input_excel = root_dir / input_excel_name
    output_excel = root_dir / output_excel_name
    figures_dir = root_dir / "Analysis_Figures"

    figures_dir.mkdir(exist_ok=True)

    if not input_excel.exists():
        raise FileNotFoundError(f"Input file not found: {input_excel}")

    print(f"Loading Excel sheets from: {input_excel}")
    df_summary = pd.read_excel(input_excel, sheet_name="Peptide_Summary")
    df_detailed = pd.read_excel(input_excel, sheet_name="Detailed_Interactions")

    # 2. Receptor key residue frequency and SEM calculation
    def extract_residues(text):
        """Extract residue list from comma-separated string."""
        if pd.isna(text) or text == "None" or text == "No Interactions Found":
            return []
        return [r.strip() for r in str(text).split(',') if r.strip()]

    all_receptor_residues = []
    df_summary['Receptor_Residues_ChainA'].apply(lambda x: all_receptor_residues.extend(extract_residues(x)))

    rec_counts = Counter(all_receptor_residues)
    df_rec_freq = pd.DataFrame(rec_counts.items(), columns=["Receptor_Residue", "Frequency"]).sort_values(by="Frequency", ascending=False)

    reference_residues = ["ALA354", "ALA356", "GLN281", "LYS511", "TYR520", "TYR523", "HIS353", "HIS513", "HIS383", "HIS387", "GLU384", "GLU411"]
    df_rec_freq['Is_Reference_Residue'] = df_rec_freq['Receptor_Residue'].apply(lambda x: "Yes" if x in reference_residues else "No")

    # Binomial SEM calculation
    total_models = len(df_summary)
    p_freq = df_rec_freq['Frequency'] / total_models
    df_rec_freq['SEM'] = np.sqrt(p_freq * (1 - p_freq) * total_models)

    # Chi-square test for residue frequencies
    if not df_rec_freq.empty and len(df_rec_freq['Frequency']) > 1:
        chi2_stat, chi2_p_val = stats.chisquare(df_rec_freq['Frequency'])
    else:
        chi2_stat, chi2_p_val = 0.0, 1.0

    # 3. Peptide residue pattern analysis
    all_peptide_residues = []
    df_summary['Peptide_Residues_ChainB'].apply(lambda x: all_peptide_residues.extend(extract_residues(x)))

    pep_counts = Counter(all_peptide_residues)
    df_pep_freq = pd.DataFrame(pep_counts.items(), columns=["Peptide_Residue_Position", "Frequency"]).sort_values(by="Frequency", ascending=False)

    # 4. Interaction type distribution
    if "Interaction_Type" in df_detailed.columns:
        df_valid_detailed = df_detailed[df_detailed["Interaction_Type"] != "-"]
        df_interaction_type = df_valid_detailed["Interaction_Type"].value_counts().reset_index()
        df_interaction_type.columns = ["Interaction_Type", "Count"]
        df_interaction_type["Percentage (%)"] = (df_interaction_type["Count"] / df_interaction_type["Count"].sum()) * 100
    else:
        df_interaction_type = pd.DataFrame(columns=["Interaction_Type", "Count", "Percentage (%)"])

    # 5. Correlation analysis: docking score vs total interactions
    def parse_score(val):
        """Parse numeric docking score values."""
        try:
            if isinstance(val, (int, float)):
                return float(val)
            val_str = str(val).strip()
            if '/' in val_str:
                parts = val_str.split('/')
                return float(parts[0]) + (float(parts[1]) / 1000.0 if float(parts[0]) < 0 else -float(parts[1]) / 1000.0)
            return float(val_str)
        except Exception:
            return np.nan

    df_summary['Score_Numeric'] = df_summary['Docking_Score'].apply(parse_score)

    df_corr_clean = df_summary.dropna(subset=['Score_Numeric', 'Total_Interactions'])
    if len(df_corr_clean) > 2:
        corr_r, p_val_corr = stats.pearsonr(df_corr_clean['Score_Numeric'], df_corr_clean['Total_Interactions'])
    else:
        corr_r, p_val_corr = 0.0, 1.0

    df_stats_summary = pd.DataFrame([
        {
            "Analysis": "Pearson Correlation (Score vs Interactions)", 
            "Statistic_r": corr_r, 
            "p_value": p_val_corr,
            "Is_Significant (p < 0.05)": "Yes" if p_val_corr < 0.05 else "No"
        },
        {
            "Analysis": "Chi-Square Test (Receptor Residues Frequency)", 
            "Statistic_r": chi2_stat, 
            "p_value": chi2_p_val,
            "Is_Significant (p < 0.05)": "Yes" if chi2_p_val < 0.05 else "No"
        }
    ])

    # 6. Candidate lead selection
    def calc_ref_coverage(text):
        """Count reference residue matches."""
        res_list = extract_residues(text)
        return sum(1 for r in res_list if r in reference_residues)

    df_summary['Ref_Residue_Matches'] = df_summary['Receptor_Residues_ChainA'].apply(calc_ref_coverage)

    df_summary['Rank_Score'] = (
        df_summary['Score_Numeric'].rank(ascending=True) * 0.4 +
        df_summary['Total_Interactions'].rank(ascending=False) * 0.3 +
        df_summary['Ref_Residue_Matches'].rank(ascending=False) * 0.3
    )
    df_leads = df_summary.sort_values(by="Rank_Score", ascending=True).drop(columns=['Rank_Score'])

    # 7. Generate and save figures
    sns.set_theme(style="whitegrid")

    # Figure 1: Receptor residue frequencies with SEM bars
    fig1_path = figures_dir / "Fig1_Receptor_Residues_Frequency.png"
    df_top12 = df_rec_freq.head(12)

    plt.figure(figsize=(10, 5))
    palette_colors = sns.color_palette("viridis", len(df_top12))

    plt.bar(
        df_top12["Receptor_Residue"], 
        df_top12["Frequency"], 
        yerr=df_top12["SEM"], 
        capsize=4, 
        color=palette_colors, 
        edgecolor='black', 
        linewidth=0.5,
        error_kw={'ecolor': '#333333', 'elinewidth': 1.2, 'capthick': 1.2}
    )

    plt.xlabel("Receptor Residue", fontsize=10)
    plt.ylabel("Frequency", fontsize=10)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(fig1_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Figure 1 saved: {fig1_path}")

    # Figure 2: Docking score vs interactions correlation regression
    fig2_path = figures_dir / "Fig2_Score_vs_Interactions_Correlation.png"
    plt.figure(figsize=(8, 5))
    ax = sns.regplot(data=df_summary, x="Score_Numeric", y="Total_Interactions", color="b", scatter_kws={'alpha': 0.7})

    p_text = "p < 0.001" if p_val_corr < 0.001 else f"p = {p_val_corr:.4f}"
    stats_annotation = f"r = {corr_r:.2f}\n{p_text}"
    ax.text(0.95, 0.95, stats_annotation, transform=ax.transAxes, fontsize=11,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8, edgecolor='gray'))

    plt.xlabel("Docking Score (More Negative = Better)", fontsize=10)
    plt.ylabel("Total Interactions Count", fontsize=10)
    plt.tight_layout()
    plt.savefig(fig2_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Figure 2 saved: {fig2_path}")

    # Figure 3: Interaction types distribution pie chart
    if not df_interaction_type.empty:
        fig3_path = figures_dir / "Fig3_Interaction_Types_Distribution.png"
        plt.figure(figsize=(7, 7))
        plt.pie(df_interaction_type["Count"], labels=df_interaction_type["Interaction_Type"], autopct='%1.1f%%', startangle=140, colors=sns.color_palette("pastel"))
        plt.tight_layout()
        plt.savefig(fig3_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Figure 3 saved: {fig3_path}")

    # 8. Save output Excel file
    output_excel_path = root_dir / output_excel
    print(f"Saving results to Excel: {output_excel_path}")
    with pd.ExcelWriter(output_excel_path, engine="openpyxl") as writer:
        df_leads.to_excel(writer, sheet_name="Top_Ranked_Peptides", index=False)
        df_rec_freq.to_excel(writer, sheet_name="Receptor_Residues_Freq", index=False)
        df_pep_freq.to_excel(writer, sheet_name="Peptide_Residues_Freq", index=False)
        df_interaction_type.to_excel(writer, sheet_name="Interaction_Types", index=False)
        df_stats_summary.to_excel(writer, sheet_name="Statistical_Analysis", index=False)

    print("\nComprehensive analysis completed successfully!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run comprehensive analysis and statistics on ACE docking interactions.")
    parser.add_argument(
        "root_dir",
        nargs="?",
        default=".",
        type=Path,
        help="Root directory containing cleaned interaction files (default: current directory)."
    )
    parser.add_argument(
        "-i", "--input",
        default="Selected_Models_Interactions_Cleaned.xlsx",
        help="Input cleaned Excel filename (default: Selected_Models_Interactions_Cleaned.xlsx)."
    )
    parser.add_argument(
        "-o", "--output",
        default="ACE_Docking_Comprehensive_Analysis.xlsx",
        help="Output comprehensive analysis Excel filename (default: ACE_Docking_Comprehensive_Analysis.xlsx)."
    )
    args = parser.parse_args()
    
    run_comprehensive_analysis(args.root_dir, args.input, args.output)
