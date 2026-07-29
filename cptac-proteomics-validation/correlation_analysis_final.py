#!/usr/bin/env python3
"""
Final Working Correlation Analysis for CPTAC RTK/NRTK Data
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
import os
warnings.filterwarnings('ignore')

print("="*80)
print("FINAL CORRELATION ANALYSIS: RTK vs NRTK")
print("="*80)

# Create output directory
os.makedirs('correlation_final_results', exist_ok=True)

# Load preprocessed data
print("\n📁 Loading preprocessed data...")
expr_data = pd.read_csv('preprocessing_results/expression_standardized.csv', index_col=0)
print(f"✅ Loaded: {expr_data.shape[0]} genes × {expr_data.shape[1]} samples")

# Load gene info
gene_stats = pd.read_csv('eda_results/gene_statistics.csv')
print(f"✅ Gene info: {len(gene_stats)} genes")

print("\n" + "="*80)
print("🔬 IDENTIFYING RTK AND NRTK GENES")
print("="*80)

# Simple identification from gene IDs
rtk_genes = []
nrtk_genes = []

# Common RTK/NRTK patterns
rtk_patterns = ['EGFR', 'ERBB', 'FGFR', 'IGF1R', 'MET', 'PDGFR', 'KIT', 'FLT3', 'ALK',
               'ROS1', 'RET', 'NTRK', 'AXL', 'MER', 'TYRO3', 'EPHA', 'EPHB', 'DDR']

nrtk_patterns = ['SRC', 'JAK', 'ABL', 'FYN', 'LYN', 'YES', 'FGR', 'HCK', 'LCK', 'BLK',
                'FRK', 'BRK', 'TEC', 'BTK', 'ITK', 'TXK', 'SYK', 'ZAP70', 'CSK', 'MATK']

for gene_id in expr_data.index:
    gene_str = str(gene_id).upper()
    
    is_rtk = any(pattern in gene_str for pattern in rtk_patterns)
    is_nrtk = any(pattern in gene_str for pattern in nrtk_patterns)
    
    if is_rtk:
        rtk_genes.append(gene_id)
    elif is_nrtk:
        nrtk_genes.append(gene_id)

print(f"\nClassification results:")
print(f"  • RTKs: {len(rtk_genes)}")
print(f"  • NRTKs: {len(nrtk_genes)}")
print(f"  • Other: {expr_data.shape[0] - len(rtk_genes) - len(nrtk_genes)}")

if rtk_genes:
    print(f"\nRTK examples (first 5):")
    for i, gene in enumerate(rtk_genes[:5], 1):
        print(f"  {i}. {str(gene)[:60]}")

if nrtk_genes:
    print(f"\nNRTK examples (first 5):")
    for i, gene in enumerate(nrtk_genes[:5], 1):
        print(f"  {i}. {str(gene)[:60]}")

# Filter to RTKs and NRTKs
analysis_genes = rtk_genes + nrtk_genes
expr_filtered = expr_data.loc[expr_data.index.isin(analysis_genes)].copy()

print(f"\n📊 Expression data for analysis:")
print(f"  • Genes: {expr_filtered.shape[0]}")
print(f"  • Samples: {expr_filtered.shape[1]}")

# Calculate correlations
print("\n" + "="*80)
print("📈 CALCULATING CORRELATIONS")
print("="*80)

corr_matrix = expr_filtered.T.corr(method='pearson')
print(f"Correlation matrix: {corr_matrix.shape[0]} × {corr_matrix.shape[1]}")

# Extract correlation pairs
corr_pairs = []
for i, gene1 in enumerate(corr_matrix.index):
    for j, gene2 in enumerate(corr_matrix.columns):
        if i < j:
            corr_value = corr_matrix.iloc[i, j]
            type1 = 'RTK' if gene1 in rtk_genes else 'NRTK'
            type2 = 'RTK' if gene2 in rtk_genes else 'NRTK'
            
            corr_pairs.append({
                'Gene1': gene1,
                'Gene2': gene2,
                'Type1': type1,
                'Type2': type2,
                'Correlation': corr_value,
                'AbsCorrelation': abs(corr_value)
            })

corr_df = pd.DataFrame(corr_pairs)

print(f"\nTotal correlation pairs: {len(corr_df)}")

# Find top correlations
print("\n" + "="*80)
print("🏆 TOP CORRELATIONS")
print("="*80)

top_positive = corr_df.sort_values('Correlation', ascending=False).head(15)
top_negative = corr_df.sort_values('Correlation', ascending=True).head(15)

print("\nTop 5 positive correlations:")
for i, row in top_positive.head(5).iterrows():
    print(f"  {row['Gene1'][:30]} - {row['Gene2'][:30]}: r = {row['Correlation']:.3f}")

print("\nTop 5 negative correlations:")
for i, row in top_negative.head(5).iterrows():
    print(f"  {row['Gene1'][:30]} - {row['Gene2'][:30]}: r = {row['Correlation']:.3f}")

# RTK-NRTK correlations
rtk_nrtk_df = corr_df[(corr_df['Type1'] != corr_df['Type2'])].copy()
if not rtk_nrtk_df.empty:
    top_rtk_nrtk = rtk_nrtk_df.sort_values('AbsCorrelation', ascending=False).head(15)
    
    print("\nTop 5 RTK-NRTK correlations:")
    for i, row in top_rtk_nrtk.head(5).iterrows():
        print(f"  {row['Gene1'][:30]} - {row['Gene2'][:30]}: r = {row['Correlation']:.3f}")

# Statistical significance
print("\n" + "="*80)
print("📊 STATISTICAL SIGNIFICANCE")
print("="*80)

significant_corrs = []
for _, row in corr_df.iterrows():
    gene1 = row['Gene1']
    gene2 = row['Gene2']
    
    vec1 = expr_filtered.loc[gene1].values
    vec2 = expr_filtered.loc[gene2].values
    
    # Remove NaN
    mask = ~(np.isnan(vec1) | np.isnan(vec2))
    vec1_clean = vec1[mask]
    vec2_clean = vec2[mask]
    
    if len(vec1_clean) > 2:
        try:
            corr_val, p_val = stats.pearsonr(vec1_clean, vec2_clean)
            if p_val < 0.05:
                significant_corrs.append({
                    'Gene1': row['Gene1'],
                    'Gene2': row['Gene2'],
                    'Type1': row['Type1'],
                    'Type2': row['Type2'],
                    'Correlation': corr_val,
                    'P-value': p_val
                })
        except:
            pass

if significant_corrs:
    sig_df = pd.DataFrame(significant_corrs)
    
    # FDR correction
    sig_df = sig_df.sort_values('P-value')
    m = len(sig_df)
    sig_df['Rank'] = range(1, m + 1)
    sig_df['FDR'] = sig_df['P-value'] * (m / sig_df['Rank'])
    
    # Ensure monotonic FDR
    for i in range(1, len(sig_df)):
        if sig_df.iloc[i]['FDR'] < sig_df.iloc[i-1]['FDR']:
            sig_df.at[sig_df.index[i], 'FDR'] = sig_df.iloc[i-1]['FDR']
    
    # Count by type
    sig_rtk_rtk = sig_df[(sig_df['Type1'] == 'RTK') & (sig_df['Type2'] == 'RTK')]
    sig_nrtk_nrtk = sig_df[(sig_df['Type1'] == 'NRTK') & (sig_df['Type2'] == 'NRTK')]
    sig_rtk_nrtk = sig_df[(sig_df['Type1'] != sig_df['Type2'])]
    
    print(f"\nSignificant correlations (p < 0.05): {len(sig_df)}")
    print(f"  • RTK-RTK: {len(sig_rtk_rtk)}")
    print(f"  • NRTK-NRTK: {len(sig_nrtk_nrtk)}")
    print(f"  • RTK-NRTK: {len(sig_rtk_nrtk)}")
    
    print("\nTop 5 most significant correlations:")
    for i, row in sig_df.head(5).iterrows():
        print(f"  {row['Gene1'][:30]} - {row['Gene2'][:30]}: r = {row['Correlation']:.3f}, p = {row['P-value']:.2e}, FDR = {row['FDR']:.2e}")
    
else:
    print("\nNo significant correlations found (p < 0.05)")

# Visualization
print("\n" + "="*80)
print("📊 CREATING VISUALIZATIONS")
print("="*80)

# 1. Correlation heatmap
plt.figure(figsize=(12, 10))
mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
sns.heatmap(corr_matrix, mask=mask, cmap='RdBu_r', center=0,
            square=True, cbar_kws={'shrink': 0.8})
plt.title('RTK/NRTK Correlation Matrix', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('correlation_final_results/correlation_heatmap.png', dpi=150)
print("✓ Saved: correlation_final_results/correlation_heatmap.png")

# 2. Distribution plots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Get correlations by type
rtk_rtk_corrs = corr_df[(corr_df['Type1'] == 'RTK') & (corr_df['Type2'] == 'RTK')]['Correlation'].values
nrtk_nrtk_corrs = corr_df[(corr_df['Type1'] == 'NRTK') & (corr_df['Type2'] == 'NRTK')]['Correlation'].values
rtk_nrtk_corrs = corr_df[(corr_df['Type1'] != corr_df['Type2'])]['Correlation'].values

# Plot RTK-RTK
if len(rtk_rtk_corrs) > 0:
    axes[0].hist(rtk_rtk_corrs, bins=20, edgecolor='black', alpha=0.7, color='#ff6b6b')
    mean_val = np.mean(rtk_rtk_corrs)
    axes[0].axvline(mean_val, color='#c92a2a', linestyle='--', linewidth=2,
                   label=f'Mean: {mean_val:.3f}')
    axes[0].set_xlabel('Correlation Coefficient')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title(f'RTK-RTK Correlations\n(n={len(rtk_rtk_corrs)} pairs)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

# Plot NRTK-NRTK
if len(nrtk_nrtk_corrs) > 0:
    axes[1].hist(nrtk_nrtk_corrs, bins=20, edgecolor='black', alpha=0.7, color='#4d96ff')
    mean_val = np.mean(nrtk_nrtk_corrs)
    axes[1].axvline(mean_val, color='#1e56a0', linestyle='--', linewidth=2,
                   label=f'Mean: {mean_val:.3f}')
    axes[1].set_xlabel('Correlation Coefficient')
    axes[1].set_ylabel('Frequency')
    axes[1].set_title(f'NRTK-NRTK Correlations\n(n={len(nrtk_nrtk_corrs)} pairs)')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

# Plot RTK-NRTK
if len(rtk_nrtk_corrs) > 0:
    axes[2].hist(rtk_nrtk_corrs, bins=20, edgecolor='black', alpha=0.7, color='#9d65c9')
    mean_val = np.mean(rtk_nrtk_corrs)
    axes[2].axvline(mean_val, color='#6a1b9a', linestyle='--', linewidth=2,
                   label=f'Mean: {mean_val:.3f}')
    axes[2].set_xlabel('Correlation Coefficient')
    axes[2].set_ylabel('Frequency')
    axes[2].set_title(f'RTK-NRTK Cross-Correlations\n(n={len(rtk_nrtk_corrs)} pairs)')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('correlation_final_results/correlation_distributions.png', dpi=150)
print("✓ Saved: correlation_final_results/correlation_distributions.png")

# 3. Scatter plot of top RTK-NRTK correlations
if 'top_rtk_nrtk' in locals() and len(top_rtk_nrtk) > 0:
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    for idx, (_, row) in enumerate(top_rtk_nrtk.head(6).iterrows()):
        gene1 = row['Gene1']
        gene2 = row['Gene2']
        
        x_vals = expr_filtered.loc[gene1].values
        y_vals = expr_filtered.loc[gene2].values
        
        # Remove NaN
        mask = ~(np.isnan(x_vals) | np.isnan(y_vals))
        x_clean = x_vals[mask]
        y_clean = y_vals[mask]
        
        if len(x_clean) > 2:
            # Calculate correlation
            corr_val, p_val = stats.pearsonr(x_clean, y_clean)
            
            # Plot
            axes[idx].scatter(x_clean, y_clean, alpha=0.6, s=30)
            
            # Add trend line
            z = np.polyfit(x_clean, y_clean, 1)
            p = np.poly1d(z)
            axes[idx].plot(x_clean, p(x_clean), "r-", linewidth=2)
            
            # Shorten gene names for title
            name1 = str(gene1).split('|')[-1] if '|' in str(gene1) else str(gene1)[:15]
            name2 = str(gene2).split('|')[-1] if '|' in str(gene2) else str(gene2)[:15]
            
            axes[idx].set_xlabel(f'{name1}')
            axes[idx].set_ylabel(f'{name2}')
            axes[idx].set_title(f'r = {corr_val:.3f}\np = {p_val:.2e}')
            axes[idx].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('correlation_final_results/top_rtk_nrtk_scatters.png', dpi=150)
    print("✓ Saved: correlation_final_results/top_rtk_nrtk_scatters.png")

# 4. Bar chart of top correlations
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Top positive
if len(top_positive) > 0:
    top_pos_plot = top_positive.head(10).copy()
    # Shorten gene names
    top_pos_plot['Label'] = top_pos_plot.apply(
        lambda row: f"{str(row['Gene1']).split('|')[-1][:10]}-{str(row['Gene2']).split('|')[-1][:10]}", 
        axis=1
    )
    
    axes[0].barh(range(len(top_pos_plot)), top_pos_plot['Correlation'].values, 
                color=['#ff6b6b' if t1=='RTK' and t2=='RTK' else 
                      '#4d96ff' if t1=='NRTK' and t2=='NRTK' else 
                      '#9d65c9' for t1, t2 in zip(top_pos_plot['Type1'], top_pos_plot['Type2'])])
    axes[0].set_yticks(range(len(top_pos_plot)))
    axes[0].set_yticklabels(top_pos_plot['Label'].values)
    axes[0].set_xlabel('Correlation Coefficient')
    axes[0].set_title('Top 10 Positive Correlations')
    axes[0].invert_yaxis()
    axes[0].grid(True, alpha=0.3, axis='x')

# Top negative
if len(top_negative) > 0:
    top_neg_plot = top_negative.head(10).copy()
    top_neg_plot['Label'] = top_neg_plot.apply(
        lambda row: f"{str(row['Gene1']).split('|')[-1][:10]}-{str(row['Gene2']).split('|')[-1][:10]}", 
        axis=1
    )
    
    axes[1].barh(range(len(top_neg_plot)), top_neg_plot['Correlation'].values,
                color=['#ff6b6b' if t1=='RTK' and t2=='RTK' else 
                      '#4d96ff' if t1=='NRTK' and t2=='NRTK' else 
                      '#9d65c9' for t1, t2 in zip(top_neg_plot['Type1'], top_neg_plot['Type2'])])
    axes[1].set_yticks(range(len(top_neg_plot)))
    axes[1].set_yticklabels(top_neg_plot['Label'].values)
    axes[1].set_xlabel('Correlation Coefficient')
    axes[1].set_title('Top 10 Negative Correlations')
    axes[1].invert_yaxis()
    axes[1].grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig('correlation_final_results/top_correlations_barchart.png', dpi=150)
print("✓ Saved: correlation_final_results/top_correlations_barchart.png")

# Save results
print("\n" + "="*80)
print("💾 SAVING RESULTS")
print("="*80)

# Save all files
corr_matrix.to_csv('correlation_final_results/full_correlation_matrix.csv')
corr_df.to_csv('correlation_final_results/all_correlation_pairs.csv', index=False)
top_positive.to_csv('correlation_final_results/top_positive_correlations.csv', index=False)
top_negative.to_csv('correlation_final_results/top_negative_correlations.csv', index=False)

if 'rtk_nrtk_df' in locals() and not rtk_nrtk_df.empty:
    top_rtk_nrtk.to_csv('correlation_final_results/top_rtk_nrtk_correlations.csv', index=False)

if 'sig_df' in locals():
    sig_df.to_csv('correlation_final_results/significant_correlations.csv', index=False)

print("✅ All results saved in 'correlation_final_results/' folder")

# Create final summary
summary = f"""
CPTAC BRCA RTK/NRTK CORRELATION ANALYSIS - FINAL RESULTS
========================================================
Analysis Date: {pd.Timestamp.now()}

DATA SUMMARY
------------
• Total genes analyzed: {expr_filtered.shape[0]}
• RTKs identified: {len(rtk_genes)}
• NRTKs identified: {len(nrtk_genes)}
• Breast cancer samples: {expr_filtered.shape[1]}
• Total correlation pairs analyzed: {len(corr_df)}

KEY FINDINGS
------------
• Significant correlations (p < 0.05): {len(sig_df) if 'sig_df' in locals() else 0}
  - RTK-RTK pairs: {len(sig_rtk_rtk) if 'sig_rtk_rtk' in locals() else 0}
  - NRTK-NRTK pairs: {len(sig_nrtk_nrtk) if 'sig_nrtk_nrtk' in locals() else 0}
  - RTK-NRTK pairs: {len(sig_rtk_nrtk) if 'sig_rtk_nrtk' in locals() else 0}

• Strongest positive correlation: {top_positive.iloc[0]['Gene1'][:30]} - {top_positive.iloc[0]['Gene2'][:30]}
  r = {top_positive.iloc[0]['Correlation']:.3f}

• Strongest negative correlation: {top_negative.iloc[0]['Gene1'][:30]} - {top_negative.iloc[0]['Gene2'][:30]}
  r = {top_negative.iloc[0]['Correlation']:.3f}

"""

if 'top_rtk_nrtk' in locals() and len(top_rtk_nrtk) > 0:
    summary += f"""• Strongest RTK-NRTK correlation: {top_rtk_nrtk.iloc[0]['Gene1'][:30]} - {top_rtk_nrtk.iloc[0]['Gene2'][:30]}
  r = {top_rtk_nrtk.iloc[0]['Correlation']:.3f}

"""

summary += f"""
BIOLOGICAL INTERPRETATION
-------------------------
1. RTK-RTK correlations suggest co-regulation of receptor tyrosine kinases
2. RTK-NRTK correlations indicate potential signaling interactions
3. Negative correlations may suggest compensatory mechanisms

FILES GENERATED
---------------
1. full_correlation_matrix.csv - Complete correlation matrix
2. all_correlation_pairs.csv - All pairwise correlations
3. significant_correlations.csv - Statistically significant pairs
4. Various visualization PNG files

NEXT STEPS
----------
1. Validate top correlations in independent datasets
2. Investigate biological pathways for top RTK-NRTK pairs
3. Explore clinical implications for breast cancer therapy
4. Consider experimental validation of novel interactions
"""

with open('correlation_final_results/analysis_summary.txt', 'w') as f:
    f.write(summary)

print("✓ Saved: correlation_final_results/analysis_summary.txt")

print("\n" + "="*80)
print("🎉 ANALYSIS COMPLETE! 🎉")
print("="*80)
print("\n📊 RESULTS SUMMARY:")
print(f"   • Analyzed {expr_filtered.shape[0]} kinases in {expr_filtered.shape[1]} breast cancer samples")
print(f"   • Found {len(sig_df) if 'sig_df' in locals() else 0} statistically significant correlations")
print(f"   • Discovered {len(sig_rtk_nrtk) if 'sig_rtk_nrtk' in locals() else 0} significant RTK-NRTK interactions")
print("\n📁 Check 'correlation_final_results/' folder for:")
print("   • Correlation matrices and tables")
print("   • Statistical significance results")
print("   • Publication-ready visualizations")
print("\n🔬 These findings could reveal novel therapeutic targets")
print("   for breast cancer treatment!")
