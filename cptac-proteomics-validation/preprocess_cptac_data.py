#!/usr/bin/env python3
"""
Fixed preprocessing pipeline for CPTAC BRCA RTK/NRTK data
"""

import pandas as pd
import numpy as np
from sklearn.impute import KNNImputer
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("FIXED CPTAC DATA PREPROCESSING PIPELINE")
print("="*80)

# Create output directory first
os.makedirs('preprocessing_results', exist_ok=True)

# Load the cleaned expression matrix from EDA
print("\n📁 Loading data from EDA results...")
try:
    expr_matrix = pd.read_csv('eda_results/expression_matrix.csv', index_col=0)
    print(f"✅ Loaded expression matrix: {expr_matrix.shape[0]} genes × {expr_matrix.shape[1]} samples")
except:
    print("❌ Could not load EDA results. Loading original data...")
    df = pd.read_csv('BRCA_UMICH_proteomics_RTK_NRTK.csv')
    
    sample_cols = []
    for col in df.columns:
        if isinstance(col, str):
            if ('BR' in col and any(c.isdigit() for c in col[:3])) or 'CPT' in col:
                sample_cols.append(col)
    
    expr_matrix = df[sample_cols].copy()
    expr_matrix.index = df[df.columns[0]]
    expr_matrix = expr_matrix.apply(pd.to_numeric, errors='coerce')
    print(f"Loaded: {expr_matrix.shape[0]} genes × {expr_matrix.shape[1]} samples")

# Step 1: Remove genes with >30% missing data
print("\n" + "="*80)
print("STEP 1: FILTER GENES WITH HIGH MISSING DATA (>30%)")
print("="*80)

missing_by_gene = (expr_matrix.isnull().sum(axis=1) / expr_matrix.shape[1]) * 100
genes_before = expr_matrix.shape[0]

filtered_expr = expr_matrix[missing_by_gene <= 30].copy()
genes_after = filtered_expr.shape[0]

print(f"\nGene filtering results:")
print(f"  • Before filtering: {genes_before} genes")
print(f"  • After filtering (≤30% missing): {genes_after} genes")
print(f"  • Removed: {genes_before - genes_after} genes ({((genes_before - genes_after)/genes_before)*100:.1f}%)")

# Step 2: Impute missing values
print("\n" + "="*80)
print("STEP 2: IMPUTE MISSING VALUES")
print("="*80)

print(f"\nMissing values before imputation: {filtered_expr.isnull().sum().sum():,}")

# Simple mean imputation first
expr_mean_imputed = filtered_expr.copy()
for gene in expr_mean_imputed.index:
    gene_mean = expr_mean_imputed.loc[gene].mean()
    expr_mean_imputed.loc[gene] = expr_mean_imputed.loc[gene].fillna(gene_mean)

print(f"Missing after mean imputation: {expr_mean_imputed.isnull().sum().sum():,}")

# Try KNN imputation if we have enough data
if filtered_expr.shape[0] > 5 and filtered_expr.shape[1] > 5:
    try:
        knn_data = filtered_expr.T
        knn_imputer = KNNImputer(n_neighbors=min(5, knn_data.shape[0]-1), weights='uniform')
        knn_imputed = knn_imputer.fit_transform(knn_data)
        expr_knn_imputed = pd.DataFrame(knn_imputed.T, 
                                       index=filtered_expr.index, 
                                       columns=filtered_expr.columns)
        print(f"Missing after KNN imputation: {expr_knn_imputed.isnull().sum().sum():,}")
        expr_imputed = expr_knn_imputed
        imputation_method = "KNN"
    except:
        expr_imputed = expr_mean_imputed
        imputation_method = "Gene mean"
else:
    expr_imputed = expr_mean_imputed
    imputation_method = "Gene mean"

print(f"\n✅ Selected imputation method: {imputation_method}")

# Step 3: Check skewness
print("\n" + "="*80)
print("STEP 3: CHECK DATA DISTRIBUTION")
print("="*80)

values_before = expr_imputed.values.flatten()
values_before = values_before[~np.isnan(values_before)]
skewness_before = pd.Series(values_before).skew()

print(f"\nDistribution before transformation:")
print(f"  • Skewness: {skewness_before:.3f}")
print(f"  • Data is {'approximately symmetric' if abs(skewness_before) <= 1 else 'skewed'}")

if abs(skewness_before) > 1:
    print("Applying log2 transformation...")
    min_value = expr_imputed.min().min()
    if min_value <= 0:
        constant = abs(min_value) + 1e-6
        expr_transformed = np.log2(expr_imputed + constant)
    else:
        expr_transformed = np.log2(expr_imputed)
    needs_log_transform = True
else:
    expr_transformed = expr_imputed.copy()
    needs_log_transform = False
    print("Skipping log transformation.")

# Step 4: Standardize data
print("\n" + "="*80)
print("STEP 4: STANDARDIZE DATA")
print("="*80)

expr_standardized = expr_transformed.copy()
for gene in expr_standardized.index:
    gene_values = expr_standardized.loc[gene]
    mean_val = gene_values.mean()
    std_val = gene_values.std()
    if std_val > 0:
        expr_standardized.loc[gene] = (gene_values - mean_val) / std_val

print(f"Standardization complete:")
print(f"  • Mean: {expr_standardized.values.mean():.6f}")
print(f"  • Std: {expr_standardized.values.std():.3f}")

# Step 5: Check for batch effects (simplified)
print("\n" + "="*80)
print("STEP 5: CHECK FOR BATCH EFFECTS")
print("="*80)

samples = expr_standardized.columns.tolist()
batches = {}
for sample in samples:
    if isinstance(sample, str):
        if 'BR' in sample:
            parts = sample.split('BR')
            if parts[0].isdigit():
                batch = int(parts[0])
                batches.setdefault(batch, []).append(sample)
        elif 'CPT' in sample:
            batches.setdefault('CPT', []).append(sample)

print(f"\nIdentified {len(batches)} potential batches:")
for batch in list(batches.keys())[:10]:  # Show first 10
    print(f"  • Batch {batch}: {len(batches[batch])} samples")

# Simple batch visualization
if len(batches) > 1:
    plt.figure(figsize=(10, 6))
    
    # Calculate mean expression per batch
    batch_means = []
    batch_labels = []
    
    for batch, batch_samples in batches.items():
        if len(batch_samples) > 1:
            batch_data = expr_standardized[batch_samples].mean(axis=1).mean()
            batch_means.append(batch_data)
            batch_labels.append(str(batch))
    
    plt.bar(range(len(batch_means)), batch_means)
    plt.xticks(range(len(batch_means)), batch_labels, rotation=45)
    plt.xlabel('Batch')
    plt.ylabel('Mean Expression (Z-score)')
    plt.title('Mean Expression by Batch')
    plt.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig('preprocessing_results/batch_effects.png', dpi=150)
    print("✓ Saved: preprocessing_results/batch_effects.png")

# Step 6: Save preprocessed data
print("\n" + "="*80)
print("STEP 6: SAVE PREPROCESSED DATA")
print("="*80)

# Save all versions
filtered_expr.to_csv('preprocessing_results/expression_filtered.csv')
expr_imputed.to_csv('preprocessing_results/expression_imputed.csv')
expr_standardized.to_csv('preprocessing_results/expression_standardized.csv')

if needs_log_transform:
    expr_transformed.to_csv('preprocessing_results/expression_log_transformed.csv')

print(f"\n✅ Saved preprocessed datasets:")
print(f"  • preprocessing_results/expression_filtered.csv")
print(f"  • preprocessing_results/expression_imputed.csv")
print(f"  • preprocessing_results/expression_standardized.csv")
if needs_log_transform:
    print(f"  • preprocessing_results/expression_log_transformed.csv")

# Create summary
summary = f"""
CPTAC DATA PREPROCESSING SUMMARY
================================
Generated: {pd.Timestamp.now()}

PREPROCESSING STEPS:
1. Gene filtering: Removed genes with >30% missing data
   - Before: {genes_before} genes
   - After: {genes_after} genes

2. Missing value imputation: {imputation_method}
   - Missing before: {filtered_expr.isnull().sum().sum():,}
   - Missing after: {expr_imputed.isnull().sum().sum():,}

3. Data transformation: {'log2 applied' if needs_log_transform else 'none'}
   - Skewness: {skewness_before:.3f}

4. Standardization: Z-score normalization
   - Final mean: {expr_standardized.values.mean():.6f}
   - Final std: {expr_standardized.values.std():.3f}

5. Batch effects: {len(batches)} batches identified

FINAL DATASET:
• Genes: {expr_standardized.shape[0]}
• Samples: {expr_standardized.shape[1]}

RECOMMENDED FILE:
• preprocessing_results/expression_standardized.csv
"""

with open('preprocessing_results/preprocessing_summary.txt', 'w') as f:
    f.write(summary)

print(f"\n✓ Saved: preprocessing_results/preprocessing_summary.txt")

# Quick visualization
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Before/after distribution
axes[0].hist(values_before, bins=30, alpha=0.7, label='Original', edgecolor='black')
axes[0].hist(expr_standardized.values.flatten(), bins=30, alpha=0.7, label='Standardized', edgecolor='black')
axes[0].set_xlabel('Expression Value')
axes[0].set_ylabel('Frequency')
axes[0].set_title('Distribution: Before vs After')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Heatmap of standardized data (subset)
if expr_standardized.shape[0] > 5:
    subset = expr_standardized.iloc[:10, :10]
    im = axes[1].imshow(subset, aspect='auto', cmap='coolwarm', vmin=-3, vmax=3)
    axes[1].set_xlabel('Samples')
    axes[1].set_ylabel('Genes')
    axes[1].set_title('Standardized Data (Subset)')
    plt.colorbar(im, ax=axes[1])

plt.tight_layout()
plt.savefig('preprocessing_results/preprocessing_overview.png', dpi=150)
print("✓ Saved: preprocessing_results/preprocessing_overview.png")

print("\n" + "="*80)
print("✅ PREPROCESSING COMPLETE!")
print("="*80)
print("\nRun correlation analysis:")
print("python correlation_analysis_final.py")
