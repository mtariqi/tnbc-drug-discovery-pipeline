#!/usr/bin/env python3
"""
Check data structure and fix column names
"""

import pandas as pd
import numpy as np

print("Checking data structure...")

# Load your correlation data
print("\n1. Checking correlation_final_results/all_correlation_pairs.csv")
corr_df = pd.read_csv('correlation_final_results/all_correlation_pairs.csv')
print(f"Columns: {corr_df.columns.tolist()}")
print(f"First few rows:")
print(corr_df.head())

print("\n2. Checking correlation_final_results/significant_correlations.csv")
try:
    sig_df = pd.read_csv('correlation_final_results/significant_correlations.csv')
    print(f"Columns: {sig_df.columns.tolist()}")
    print(f"First few rows:")
    print(sig_df.head())
except:
    print("File not found or error loading")

print("\n3. Checking preprocessing_results/expression_standardized.csv")
expr_data = pd.read_csv('preprocessing_results/expression_standardized.csv', index_col=0)
print(f"Shape: {expr_data.shape}")
print(f"First few gene IDs:")
print(expr_data.index[:5].tolist())
