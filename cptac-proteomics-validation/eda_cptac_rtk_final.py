#!/usr/bin/env python3
"""
Corrected EDA for CPTAC BRCA RTK/NRTK Proteomics Data
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import warnings
import os
warnings.filterwarnings('ignore')

# Set visualization style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")
plt.rcParams['figure.figsize'] = (12, 8)

print("="*80)
print("CORRECTED EDA: CPTAC BRCA RTK/NRTK PROTEOMICS DATA")
print("="*80)

# Step 1: Load the data
data_file = "BRCA_UMICH_proteomics_RTK_NRTK.csv"
print(f"\n📁 Loading data from: {data_file}")

try:
    df = pd.read_csv(data_file)
    print(f"✅ Successfully loaded: {df.shape[0]} rows × {df.shape[1]} columns")
except Exception as e:
    print(f"❌ Error loading file: {e}")
    exit(1)

# Display basic info
print("\n" + "="*80)
print("📊 DATA STRUCTURE ANALYSIS")
print("="*80)

print(f"\nDataset dimensions: {df.shape[0]} rows × {df.shape[1]} columns")

# Show column information
print("\nFirst 15 columns:")
for i, col in enumerate(df.columns[:15], 1):
    dtype = df[col].dtype
    print(f"{i:2d}. {col[:50]:50} {str(dtype):10}")

# Check for key columns
has_gene_symbol = 'GeneSymbol' in df.columns
has_kinase_class = 'KinaseClass' in df.columns
has_gene = 'Gene' in df.columns

print(f"\nKey metadata columns found:")
print(f"  • GeneSymbol: {'Yes' if has_gene_symbol else 'No'}")
print(f"  • KinaseClass: {'Yes' if has_kinase_class else 'No'}")
print(f"  • Gene: {'Yes' if has_gene else 'No'}")

# Identify sample columns (columns that look like sample IDs)
sample_cols = []
for col in df.columns:
    if isinstance(col, str):
        # Look for patterns like 11BR047, CPT0018460005, etc.
        if ('BR' in col and any(c.isdigit() for c in col[:3])) or 'CPT' in col:
            sample_cols.append(col)

print(f"\nIdentified {len(sample_cols)} sample columns")
print("Sample columns (first 10):", sample_cols[:10])

if len(sample_cols) < 5:
    print("\n⚠️  Few sample columns found. Checking alternative patterns...")
    # Look for any numeric columns that might be samples
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    print(f"Found {len(numeric_cols)} numeric columns")
    sample_cols = numeric_cols[:50]  # Take first 50 numeric columns as samples

# Create expression matrix
if sample_cols:
    print(f"\nCreating expression matrix with {len(sample_cols)} samples...")
    
    # Use the first column as gene IDs (usually contains ENSEMBL IDs)
    gene_id_col = df.columns[0]
    print(f"Using '{gene_id_col}' as gene identifier")
    
    # Create expression matrix
    expr_matrix = df[sample_cols].copy()
    expr_matrix.index = df[gene_id_col]
    
    # Add metadata if available
    metadata = {}
    if has_gene_symbol:
        metadata['GeneSymbol'] = df['GeneSymbol'].values
    if has_kinase_class:
        metadata['KinaseClass'] = df['KinaseClass'].values
    if has_gene:
        metadata['Gene'] = df['Gene'].values
    
    # Convert to numeric (coerce errors to NaN)
    print("Converting to numeric values...")
    expr_numeric = expr_matrix.apply(pd.to_numeric, errors='coerce')
    
    print(f"\nExpression matrix created:")
    print(f"  • Genes: {expr_numeric.shape[0]}")
    print(f"  • Samples: {expr_numeric.shape[1]}")
    
    # Calculate missing values
    total_cells = expr_numeric.shape[0] * expr_numeric.shape[1]
    missing_cells = expr_numeric.isnull().sum().sum()
    missing_pct = (missing_cells / total_cells) * 100
    
    print(f"  • Missing values: {missing_cells:,} ({missing_pct:.1f}%)")
    
    # Add metadata to a separate dataframe
    gene_info = pd.DataFrame({'GeneID': expr_numeric.index})
    for key, values in metadata.items():
        gene_info[key] = values
    
    # Basic statistics
    print("\n" + "="*80)
    print("📈 BASIC STATISTICS")
    print("="*80)
    
    # Remove NaN for statistics
    numeric_values = expr_numeric.values.flatten()
    numeric_values = numeric_values[~np.isnan(numeric_values)]
    
    if len(numeric_values) > 0:
        print(f"\nGlobal expression statistics:")
        print(f"  • Mean: {np.mean(numeric_values):.3f}")
        print(f"  • Median: {np.median(numeric_values):.3f}")
        print(f"  • Std: {np.std(numeric_values):.3f}")
        print(f"  • Min: {np.min(numeric_values):.3f}")
        print(f"  • Max: {np.max(numeric_values):.3f}")
        print(f"  • Skewness: {stats.skew(numeric_values):.3f}")
        print(f"  • Kurtosis: {stats.kurtosis(numeric_values):.3f}")
        
        # Gene-level statistics
        gene_stats = pd.DataFrame({
            'GeneID': expr_numeric.index,
            'Mean': expr_numeric.mean(axis=1),
            'Std': expr_numeric.std(axis=1),
            'Min': expr_numeric.min(axis=1),
            'Max': expr_numeric.max(axis=1),
            'Missing%': (expr_numeric.isnull().sum(axis=1) / expr_numeric.shape[1]) * 100
        })
        
        # Add metadata
        for col in gene_info.columns:
            if col != 'GeneID':
                gene_stats[col] = gene_info[col].values
        
        print(f"\nGene-level statistics (first 10 genes):")
        print(gene_stats.head(10).to_string())
        
        # Analyze by kinase class if available
        if has_kinase_class:
            print("\n" + "="*80)
            print("🔬 KINASE CLASS ANALYSIS")
            print("="*80)
            
            # Count by kinase class
            kinase_counts = gene_stats['KinaseClass'].value_counts()
            print("\nKinase class distribution:")
            for kin_class, count in kinase_counts.items():
                if pd.isna(kin_class):
                    print(f"  • Unknown/Other: {count}")
                else:
                    print(f"  • {kin_class}: {count}")
            
            # Analyze RTKs and NRTKs
            rtks = gene_stats[gene_stats['KinaseClass'] == 'RTK']
            nrtks = gene_stats[gene_stats['KinaseClass'] == 'NRTK']
            
            print(f"\nReceptor Tyrosine Kinases (RTKs): {len(rtks)} found")
            if not rtks.empty:
                print("Top 10 RTKs by mean expression:")
                top_rtks = rtks.sort_values('Mean', ascending=False).head(10)
                for i, (_, row) in enumerate(top_rtks.iterrows(), 1):
                    gene_name = row.get('GeneSymbol', row.get('Gene', row['GeneID'][:30]))
                    print(f"  {i:2d}. {gene_name[:30]:30} Mean: {row['Mean']:.3f}")
            
            print(f"\nNon-Receptor Tyrosine Kinases (NRTKs): {len(nrtks)} found")
            if not nrtks.empty:
                print("Top 10 NRTKs by mean expression:")
                top_nrtks = nrtks.sort_values('Mean', ascending=False).head(10)
                for i, (_, row) in enumerate(top_nrtks.iterrows(), 1):
                    gene_name = row.get('GeneSymbol', row.get('Gene', row['GeneID'][:30]))
                    print(f"  {i:2d}. {gene_name[:30]:30} Mean: {row['Mean']:.3f}")
        
        # Visualization
        print("\n" + "="*80)
        print("📊 CREATING VISUALIZATIONS")
        print("="*80)
        
        os.makedirs('eda_results', exist_ok=True)
        
        # Plot 1: Distribution of expression values
        plt.figure(figsize=(10, 6))
        plt.hist(numeric_values, bins=50, edgecolor='black', alpha=0.7)
        plt.axvline(np.mean(numeric_values), color='red', linestyle='--', 
                   label=f'Mean: {np.mean(numeric_values):.3f}')
        plt.axvline(np.median(numeric_values), color='green', linestyle='--', 
                   label=f'Median: {np.median(numeric_values):.3f}')
        plt.xlabel('Expression Value')
        plt.ylabel('Frequency')
        plt.title('Distribution of All Expression Values')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('eda_results/expression_distribution.png', dpi=150)
        print("✓ Saved: eda_results/expression_distribution.png")
        
        # Plot 2: Missing data analysis
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Missing by gene
        missing_by_gene = (expr_numeric.isnull().sum(axis=1) / expr_numeric.shape[1]) * 100
        axes[0].hist(missing_by_gene, bins=30, edgecolor='black', alpha=0.7, color='coral')
        axes[0].axvline(20, color='red', linestyle='--', label='20% threshold')
        axes[0].set_xlabel('Missing Data (%) per Gene')
        axes[0].set_ylabel('Number of Genes')
        axes[0].set_title('Missing Data Distribution - Genes')
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # Missing by sample
        missing_by_sample = (expr_numeric.isnull().sum(axis=0) / expr_numeric.shape[0]) * 100
        axes[1].bar(range(len(missing_by_sample)), missing_by_sample.values)
        axes[1].set_xlabel('Sample Index')
        axes[1].set_ylabel('Missing Data (%)')
        axes[1].set_title('Missing Data Distribution - Samples')
        axes[1].grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig('eda_results/missing_data_analysis.png', dpi=150)
        print("✓ Saved: eda_results/missing_data_analysis.png")
        
        # Plot 3: Top expressed genes
        plt.figure(figsize=(12, 8))
        top_genes_all = gene_stats.sort_values('Mean', ascending=False).head(20)
        
        # Create labels with gene symbols if available
        labels = []
        for _, row in top_genes_all.iterrows():
            if 'GeneSymbol' in row and not pd.isna(row['GeneSymbol']):
                label = row['GeneSymbol']
            elif 'Gene' in row and not pd.isna(row['Gene']):
                label = row['Gene']
            else:
                label = str(row['GeneID'])[:20]
            
            if 'KinaseClass' in row and not pd.isna(row['KinaseClass']):
                label += f" ({row['KinaseClass']})"
            labels.append(label)
        
        y_pos = np.arange(len(top_genes_all))
        plt.barh(y_pos, top_genes_all['Mean'].values)
        plt.yticks(y_pos, labels, fontsize=9)
        plt.xlabel('Mean Expression')
        plt.title('Top 20 Most Expressed Genes')
        plt.gca().invert_yaxis()
        plt.grid(True, alpha=0.3, axis='x')
        plt.tight_layout()
        plt.savefig('eda_results/top_expressed_genes.png', dpi=150)
        print("✓ Saved: eda_results/top_expressed_genes.png")
        
        # Plot 4: Kinase class comparison (if available)
        if has_kinase_class:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            
            # Pie chart of kinase classes
            kinase_counts = gene_stats['KinaseClass'].value_counts()
            axes[0].pie(kinase_counts.values, labels=kinase_counts.index, autopct='%1.1f%%')
            axes[0].set_title('Kinase Class Distribution')
            
            # Box plot of expression by kinase class
            kin_classes = []
            expr_by_class = []
            
            for kin_class in ['RTK', 'NRTK']:
                if kin_class in gene_stats['KinaseClass'].values:
                    class_data = gene_stats[gene_stats['KinaseClass'] == kin_class]['Mean']
                    kin_classes.append(kin_class)
                    expr_by_class.append(class_data.values)
            
            if expr_by_class:
                axes[1].boxplot(expr_by_class, labels=kin_classes)
                axes[1].set_ylabel('Mean Expression')
                axes[1].set_title('Expression by Kinase Class')
                axes[1].grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig('eda_results/kinase_class_analysis.png', dpi=150)
            print("✓ Saved: eda_results/kinase_class_analysis.png")
        
        # Save data files
        print("\n" + "="*80)
        print("💾 SAVING DATA FILES")
        print("="*80)
        
        # Save expression matrix
        expr_numeric.to_csv('eda_results/expression_matrix.csv')
        print("✓ Saved: eda_results/expression_matrix.csv")
        
        # Save gene statistics
        gene_stats.to_csv('eda_results/gene_statistics.csv', index=False)
        print("✓ Saved: eda_results/gene_statistics.csv")
        
        # Create summary report
        report = f"""
CPTAC BRCA PROTEOMICS - EDA SUMMARY REPORT
Generated: {pd.Timestamp.now()}
Data file: {data_file}

DATA OVERVIEW
=============
• Total genes/proteins: {expr_numeric.shape[0]}
• Total samples: {expr_numeric.shape[1]}
• Total data points: {total_cells:,}
• Missing data: {missing_cells:,} ({missing_pct:.1f}%)

EXPRESSION STATISTICS
=====================
• Global mean: {np.mean(numeric_values):.3f}
• Global median: {np.median(numeric_values):.3f}
• Global std: {np.std(numeric_values):.3f}
• Range: [{np.min(numeric_values):.3f}, {np.max(numeric_values):.3f}]
• Skewness: {stats.skew(numeric_values):.3f}
• Kurtosis: {stats.kurtosis(numeric_values):.3f}

DATA QUALITY
============
• Genes with >20% missing: {(missing_by_gene > 20).sum()}
• Samples with >20% missing: {(missing_by_sample > 20).sum()}
• Complete genes (0% missing): {(missing_by_gene == 0).sum()}
• Complete samples (0% missing): {(missing_by_sample == 0).sum()}
"""
        
        if has_kinase_class:
            report += f"""
KINASE ANALYSIS
===============
• Total kinases identified: {kinase_counts.sum() if 'kinase_counts' in locals() else 'N/A'}
• Receptor Tyrosine Kinases (RTKs): {len(rtks) if 'rtks' in locals() else 0}
• Non-Receptor Tyrosine Kinases (NRTKs): {len(nrtks) if 'nrtks' in locals() else 0}

TOP EXPRESSED RTKs:
"""
            if 'rtks' in locals() and not rtks.empty:
                top_5_rtks = rtks.sort_values('Mean', ascending=False).head(5)
                for i, (_, row) in enumerate(top_5_rtks.iterrows(), 1):
                    gene_name = row.get('GeneSymbol', row.get('Gene', row['GeneID'][:30]))
                    report += f"{i}. {gene_name}: Mean={row['Mean']:.3f}, Missing={row['Missing%']:.1f}%\n"
            
            report += "\nTOP EXPRESSED NRTKs:"
            if 'nrtks' in locals() and not nrtks.empty:
                top_5_nrtks = nrtks.sort_values('Mean', ascending=False).head(5)
                for i, (_, row) in enumerate(top_5_nrtks.iterrows(), 1):
                    gene_name = row.get('GeneSymbol', row.get('Gene', row['GeneID'][:30]))
                    report += f"{i}. {gene_name}: Mean={row['Mean']:.3f}, Missing={row['Missing%']:.1f}%\n"
        
        report += f"""
FILES GENERATED
===============
1. eda_results/expression_distribution.png - Expression value distribution
2. eda_results/missing_data_analysis.png - Missing data analysis
3. eda_results/top_expressed_genes.png - Top expressed genes
4. eda_results/kinase_class_analysis.png - Kinase class analysis (if available)
5. eda_results/expression_matrix.csv - Clean expression matrix
6. eda_results/gene_statistics.csv - Gene-level statistics

RECOMMENDATIONS
===============
1. Consider removing genes with >30% missing data
2. Impute missing values using appropriate methods (kNN, mean, etc.)
3. Standardize data if comparing across samples
4. Check for batch effects among samples
5. Log-transform if data is highly skewed
"""

        with open('eda_results/eda_summary_report.txt', 'w') as f:
            f.write(report)
        
        print("✓ Saved: eda_results/eda_summary_report.txt")
        
        print("\n" + "="*80)
        print("✅ EDA COMPLETE SUCCESSFULLY!")
        print("="*80)
        print("\nAll results saved in 'eda_results/' folder")
        print("Check the summary report: eda_results/eda_summary_report.txt")
        
    else:
        print("\n❌ No numeric values found in expression matrix!")
        
else:
    print("\n❌ Could not identify sample columns for expression matrix!")
    print("Please check the data format.")
