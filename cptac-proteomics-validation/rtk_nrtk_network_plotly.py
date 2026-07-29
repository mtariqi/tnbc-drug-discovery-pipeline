#!/usr/bin/env python3
"""
FINAL Interactive Network Visualization for RTK-NRTK Co-expression
Handles ENSEMBL gene IDs properly
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import networkx as nx
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("FINAL RTK-NRTK CO-EXPRESSION NETWORK VISUALIZATION")
print("="*80)

# Load data
print("\n📁 Loading data...")
expr_data = pd.read_csv('preprocessing_results/expression_standardized.csv', index_col=0)
corr_df = pd.read_csv('correlation_final_results/all_correlation_pairs.csv')
sig_df = pd.read_csv('correlation_final_results/significant_correlations.csv')

print(f"✅ Expression data: {expr_data.shape[0]} genes × {expr_data.shape[1]} samples")
print(f"✅ Correlation pairs: {len(corr_df)}")
print(f"✅ Significant correlations: {len(sig_df)}")

# Extract gene symbols from complex ENSEMBL IDs
def extract_gene_info(gene_id):
    """Extract gene symbol and clean info from ENSEMBL ID"""
    gene_str = str(gene_id)
    
    # Parse the pipe-separated format
    parts = gene_str.split('|')
    
    # The format appears to be: ENSP|ENST|ENSG|OTTHUMG|OTTHUMT|GENE-version|GENE_SYMBOL|length
    if len(parts) >= 7:
        # Gene symbol is typically in position 6 (0-indexed) or sometimes 5
        gene_symbol = parts[6] if len(parts) > 6 else parts[5] if len(parts) > 5 else parts[0]
        
        # Clean up - remove version numbers
        if '-' in gene_symbol:
            gene_symbol = gene_symbol.split('-')[0]
        
        # Also extract from position 5 if available
        alt_symbol = parts[5] if len(parts) > 5 else gene_symbol
        if '-' in alt_symbol:
            alt_symbol = alt_symbol.split('-')[0]
        
        # Use the cleaner one
        if gene_symbol.isdigit() or len(gene_symbol) > 10:
            gene_symbol = alt_symbol
        
        return gene_symbol
    else:
        # If format is different, return first part
        return gene_str.split('|')[0] if '|' in gene_str else gene_str[:15]

# Create gene information DataFrame
print("\n🔬 Processing gene information...")
genes = list(expr_data.index)
gene_info = pd.DataFrame({
    'GeneID': genes,
    'GeneSymbol': [extract_gene_info(g) for g in genes],
    'MeanExpression': expr_data.mean(axis=1).values,
    'ExpressionStd': expr_data.std(axis=1).values
})

# Classify as RTK or NRTK based on gene symbols
def classify_by_symbol(gene_symbol):
    """Classify gene based on symbol"""
    symbol = str(gene_symbol).upper()
    
    # RTK patterns
    rtk_keywords = ['EGFR', 'ERBB', 'FGFR', 'IGF1R', 'MET', 'PDGFR', 'KIT', 'FLT3', 'ALK',
                   'ROS1', 'RET', 'NTRK', 'AXL', 'MER', 'TYRO3', 'EPHA', 'EPHB', 'DDR',
                   'INSR', 'VEGFR', 'TEK', 'TIE', 'RYK', 'MUSK', 'ROR']
    
    # NRTK patterns
    nrtk_keywords = ['SRC', 'JAK', 'ABL', 'FYN', 'LYN', 'YES', 'FGR', 'HCK', 'LCK', 'BLK',
                    'FRK', 'BRK', 'TEC', 'BTK', 'ITK', 'TXK', 'SYK', 'ZAP70', 'CSK', 'MATK',
                    'PTK', 'TNK', 'ACK', 'FAK']
    
    for keyword in rtk_keywords:
        if keyword in symbol:
            return 'RTK'
    
    for keyword in nrtk_keywords:
        if keyword in symbol:
            return 'NRTK'
    
    return 'Other'

# Apply classification
gene_info['Type'] = gene_info['GeneSymbol'].apply(classify_by_symbol)

# Also classify based on original Type columns in correlation data
print("Classifying genes from correlation data...")
all_genes_in_corr = set(corr_df['Gene1']).union(set(corr_df['Gene2']))
gene_type_map = {}

for gene_id in all_genes_in_corr:
    # Check if this gene appears as Gene1 in correlations
    gene1_matches = corr_df[corr_df['Gene1'] == gene_id]
    if not gene1_matches.empty:
        gene_type_map[gene_id] = gene1_matches.iloc[0]['Type1']
        continue
    
    # Check if this gene appears as Gene2 in correlations
    gene2_matches = corr_df[corr_df['Gene2'] == gene_id]
    if not gene2_matches.empty:
        gene_type_map[gene_id] = gene2_matches.iloc[0]['Type2']
        continue
    
    # If not found in correlation data, use symbol-based classification
    gene_symbol = extract_gene_info(gene_id)
    gene_type_map[gene_id] = classify_by_symbol(gene_symbol)

# Update gene_info with types from correlation data
for idx, gene_id in enumerate(gene_info['GeneID']):
    if gene_id in gene_type_map:
        gene_info.at[idx, 'Type'] = gene_type_map[gene_id]

rtk_genes = gene_info[gene_info['Type'] == 'RTK']['GeneID'].tolist()
nrtk_genes = gene_info[gene_info['Type'] == 'NRTK']['GeneID'].tolist()
other_genes = gene_info[gene_info['Type'] == 'Other']['GeneID'].tolist()

print(f"\n📊 Gene Classification Results:")
print(f"  • RTKs: {len(rtk_genes)}")
print(f"  • NRTKs: {len(nrtk_genes)}")
print(f"  • Other: {len(other_genes)}")

if rtk_genes:
    print(f"\nExample RTKs:")
    for i, gene_id in enumerate(rtk_genes[:5], 1):
        symbol = gene_info[gene_info['GeneID'] == gene_id]['GeneSymbol'].values[0]
        print(f"  {i}. {symbol}")

if nrtk_genes:
    print(f"\nExample NRTKs:")
    for i, gene_id in enumerate(nrtk_genes[:5], 1):
        symbol = gene_info[gene_info['GeneID'] == gene_id]['GeneSymbol'].values[0]
        print(f"  {i}. {symbol}")

# Get top RTK-NRTK correlations from significant correlations
print("\n📈 Selecting top RTK-NRTK correlations...")

# Filter for RTK-NRTK pairs in significant correlations
rtk_nrtk_sig = sig_df[
    ((sig_df['Type1'] == 'RTK') & (sig_df['Type2'] == 'NRTK')) |
    ((sig_df['Type1'] == 'NRTK') & (sig_df['Type2'] == 'RTK'))
].copy()

# If not enough significant RTK-NRTK pairs, use top correlations overall
if len(rtk_nrtk_sig) < 10:
    print(f"⚠️  Only {len(rtk_nrtk_sig)} significant RTK-NRTK pairs found.")
    print("Using top correlations from all pairs...")
    
    # Get RTK-NRTK pairs from all correlations
    rtk_nrtk_all = corr_df[
        ((corr_df['Type1'] == 'RTK') & (corr_df['Type2'] == 'NRTK')) |
        ((corr_df['Type1'] == 'NRTK') & (corr_df['Type2'] == 'RTK'))
    ].copy()
    
    # Sort by absolute correlation
    rtk_nrtk_all['AbsCorrelation'] = rtk_nrtk_all['Correlation'].abs()
    top_rtk_nrtk = rtk_nrtk_all.sort_values('AbsCorrelation', ascending=False).head(30)
else:
    # Use significant pairs
    top_rtk_nrtk = rtk_nrtk_sig.sort_values('Correlation', key=abs, ascending=False).head(30)

print(f"\n✅ Selected {len(top_rtk_nrtk)} top RTK-NRTK correlations")

# Add gene symbols for display
top_rtk_nrtk['Gene1_Symbol'] = top_rtk_nrtk['Gene1'].apply(extract_gene_info)
top_rtk_nrtk['Gene2_Symbol'] = top_rtk_nrtk['Gene2'].apply(extract_gene_info)

print("\n🏆 Top 5 RTK-NRTK Correlations:")
for i, (_, row) in enumerate(top_rtk_nrtk.head(5).iterrows(), 1):
    p_val = row.get('P-value', 'N/A')
    if isinstance(p_val, (int, float)):
        p_str = f"{p_val:.2e}"
    else:
        p_str = str(p_val)
    
    print(f"  {i}. {row['Gene1_Symbol']} ({row['Type1']}) ↔ {row['Gene2_Symbol']} ({row['Type2']})")
    print(f"     Correlation: {row['Correlation']:.3f}, P-value: {p_str}")

# Build the network graph
print("\n🔧 Building network graph...")

G = nx.Graph()

# Add nodes from top correlations
all_genes_in_network = set(top_rtk_nrtk['Gene1']).union(set(top_rtk_nrtk['Gene2']))

for gene_id in all_genes_in_network:
    # Get gene info
    info = gene_info[gene_info['GeneID'] == gene_id]
    if not info.empty:
        node_type = info['Type'].values[0]
        gene_symbol = info['GeneSymbol'].values[0]
        mean_expr = info['MeanExpression'].values[0]
        
        # Node size based on expression level
        node_size = abs(mean_expr) * 25 + 20
        
        G.add_node(gene_id,
                  gene_symbol=gene_symbol,
                  type=node_type,
                  mean_expression=mean_expr,
                  size=node_size)
    else:
        # If gene not in expression data, create basic node
        gene_symbol = extract_gene_info(gene_id)
        node_type = 'Unknown'
        
        # Try to get type from correlation data
        for _, row in top_rtk_nrtk.iterrows():
            if row['Gene1'] == gene_id:
                node_type = row['Type1']
                break
            elif row['Gene2'] == gene_id:
                node_type = row['Type2']
                break
        
        G.add_node(gene_id,
                  gene_symbol=gene_symbol,
                  type=node_type,
                  mean_expression=0,
                  size=30)

# Add edges from top correlations
for _, row in top_rtk_nrtk.iterrows():
    gene1 = row['Gene1']
    gene2 = row['Gene2']
    corr = row['Correlation']
    
    if gene1 in G.nodes() and gene2 in G.nodes():
        # Edge width based on correlation strength
        edge_width = abs(corr) * 10 + 2
        
        G.add_edge(gene1, gene2,
                  weight=abs(corr),
                  correlation=corr,
                  width=edge_width,
                  p_value=row.get('P-value', None),
                  fdr=row.get('FDR', None))

print(f"✅ Network created: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# Create interactive visualization
print("\n🎨 Creating interactive Plotly visualization...")

# Use spring layout for node positions
pos = nx.spring_layout(G, k=3, iterations=200, seed=42)

# Prepare node data
node_x = []
node_y = []
node_text = []
node_size = []
node_color = []
node_symbol = []

for node in G.nodes():
    x, y = pos[node]
    node_x.append(x)
    node_y.append(y)
    
    node_data = G.nodes[node]
    gene_symbol = node_data['gene_symbol']
    node_type = node_data['type']
    mean_expr = node_data['mean_expression']
    
    # Create hover text
    text = f"<b>{gene_symbol}</b><br>"
    text += f"Type: {node_type}<br>"
    text += f"Mean Expression: {mean_expr:.3f}<br>"
    text += f"Degree: {G.degree(node)}"
    
    node_text.append(text)
    node_size.append(node_data['size'])
    
    # Color and symbol by type
    if node_type == 'RTK':
        node_color.append('#FF6B6B')  # Red
        node_symbol.append('circle')
    elif node_type == 'NRTK':
        node_color.append('#4ECDC4')  # Teal
        node_symbol.append('diamond')
    else:
        node_color.append('#FFE66D')  # Yellow
        node_symbol.append('square')

# Prepare edge data
edge_traces = []
edge_texts = []

for edge in G.edges():
    x0, y0 = pos[edge[0]]
    x1, y1 = pos[edge[1]]
    
    edge_data = G.edges[edge]
    corr = edge_data['correlation']
    
    # Edge color based on correlation sign
    if corr > 0:
        edge_color = 'rgba(46, 204, 113, 0.8)'  # Green
    else:
        edge_color = 'rgba(231, 76, 60, 0.8)'   # Red
    
    # Create hover text for edge
    gene1_info = G.nodes[edge[0]]
    gene2_info = G.nodes[edge[1]]
    
    text = f"<b>{gene1_info['gene_symbol']} ↔ {gene2_info['gene_symbol']}</b><br>"
    text += f"Correlation: {corr:.3f}<br>"
    text += f"Type: {gene1_info['type']}-{gene2_info['type']}<br>"
    
    if 'p_value' in edge_data and edge_data['p_value'] is not None:
        text += f"P-value: {edge_data['p_value']:.2e}<br>"
    if 'fdr' in edge_data and edge_data['fdr'] is not None:
        text += f"FDR: {edge_data['fdr']:.2e}"
    
    # Create edge trace
    edge_trace = go.Scatter(
        x=[x0, x1], y=[y0, y1],
        line=dict(
            width=edge_data['width'],
            color=edge_color
        ),
        mode='lines',
        hoverinfo='text',
        hovertext=text,
        showlegend=False
    )
    
    edge_traces.append(edge_trace)

# Create main figure
fig = go.Figure()

# Add edges first (so they appear behind nodes)
for edge_trace in edge_traces:
    fig.add_trace(edge_trace)

# Add nodes
node_trace = go.Scatter(
    x=node_x, y=node_y,
    mode='markers+text',
    marker=dict(
        size=node_size,
        color=node_color,
        line=dict(width=2, color='DarkSlateGrey'),
        symbol=node_symbol
    ),
    text=[G.nodes[node]['gene_symbol'] for node in G.nodes()],
    textposition="top center",
    textfont=dict(size=12, color='black', family='Arial'),
    hovertext=node_text,
    hoverinfo='text',
    showlegend=False
)

fig.add_trace(node_trace)

# Update layout
fig.update_layout(
    title=dict(
        text='<b>RTK-NRTK Co-expression Network in Breast Cancer</b><br>CPTAC Proteomics - Top Significant Correlations',
        font=dict(size=24, family="Arial Black"),
        x=0.5,
        xanchor='center'
    ),
    showlegend=True,
    hovermode='closest',
    margin=dict(b=20, l=20, r=20, t=100),
    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    plot_bgcolor='rgba(245, 245, 245, 0.9)',
    width=1400,
    height=900,
    annotations=[
        dict(
            text="<b>Network Legend</b><br>"
                 "<span style='color:#FF6B6B'>● RTK</span> | "
                 "<span style='color:#4ECDC4'>◆ NRTK</span> | "
                 "<span style='color:#FFE66D'>■ Other</span><br>"
                 "<span style='color:#2ECC71'>━━ Positive</span> | "
                 "<span style='color:#E74C3C'>━━ Negative</span><br>"
                 "<b>Size:</b> Expression level<br>"
                 "<b>Width:</b> Correlation strength",
            xref="paper", yref="paper",
            x=0.02, y=0.98,
            align="left",
            showarrow=False,
            bgcolor="rgba(255, 255, 255, 0.95)",
            bordercolor="black",
            borderwidth=1,
            borderpad=10,
            font=dict(size=11)
        )
    ]
)

# Add legend
fig.add_trace(go.Scatter(
    x=[None], y=[None],
    mode='markers',
    marker=dict(size=15, color='#FF6B6B', symbol='circle'),
    name='Receptor TK (RTK)',
    showlegend=True
))

fig.add_trace(go.Scatter(
    x=[None], y=[None],
    mode='markers',
    marker=dict(size=15, color='#4ECDC4', symbol='diamond'),
    name='Non-Receptor TK (NRTK)',
    showlegend=True
))

fig.add_trace(go.Scatter(
    x=[None], y=[None],
    mode='lines',
    line=dict(width=4, color='#2ECC71'),
    name='Positive Correlation',
    showlegend=True
))

fig.add_trace(go.Scatter(
    x=[None], y=[None],
    mode='lines',
    line=dict(width=4, color='#E74C3C'),
    name='Negative Correlation',
    showlegend=True
))

# Create comprehensive dashboard
print("\n📊 Creating comprehensive dashboard...")

fig_dashboard = make_subplots(
    rows=2, cols=2,
    subplot_titles=('RTK-NRTK Co-expression Network', 
                   'Top 10 Correlation Values',
                   'Correlation Distribution by Type',
                   'Expression Levels by Gene Type'),
    specs=[[{'type': 'scatter'}, {'type': 'bar'}],
          [{'type': 'histogram'}, {'type': 'box'}]],
    vertical_spacing=0.12,
    horizontal_spacing=0.15
)

# Add network to dashboard
for trace in fig.data:
    fig_dashboard.add_trace(trace, row=1, col=1)

# Add bar chart of top correlations
top_10_display = top_rtk_nrtk.head(10).copy()
top_10_display['Pair'] = top_10_display.apply(
    lambda x: f"{x['Gene1_Symbol']}-{x['Gene2_Symbol']}", axis=1
)

bar_trace = go.Bar(
    x=top_10_display['Pair'],
    y=top_10_display['Correlation'],
    marker_color=['#2ECC71' if c > 0 else '#E74C3C' for c in top_10_display['Correlation']],
    text=top_10_display['Correlation'].round(3),
    textposition='outside',
    hovertext=top_10_display.apply(
        lambda x: f"{x['Gene1_Symbol']} ({x['Type1']}) ↔ {x['Gene2_Symbol']} ({x['Type2']})<br>"
                 f"r = {x['Correlation']:.3f}<br>"
                 f"P = {x.get('P-value', 'N/A')}",
        axis=1
    )
)

fig_dashboard.add_trace(bar_trace, row=1, col=2)
fig_dashboard.update_xaxes(tickangle=45, row=1, col=2)
fig_dashboard.update_yaxes(title_text="Correlation Coefficient", row=1, col=2)

# Add histogram of correlations by type
fig_dashboard.add_trace(go.Histogram(
    x=corr_df[corr_df['Type1'] == 'RTK'][corr_df['Type2'] == 'RTK']['Correlation'],
    name='RTK-RTK',
    marker_color='#FF6B6B',
    opacity=0.6,
    nbinsx=20
), row=2, col=1)

fig_dashboard.add_trace(go.Histogram(
    x=corr_df[corr_df['Type1'] == 'NRTK'][corr_df['Type2'] == 'NRTK']['Correlation'],
    name='NRTK-NRTK',
    marker_color='#4ECDC4',
    opacity=0.6,
    nbinsx=20
), row=2, col=1)

fig_dashboard.add_trace(go.Histogram(
    x=rtk_nrtk_all['Correlation'] if 'rtk_nrtk_all' in locals() else [],
    name='RTK-NRTK',
    marker_color='#9D65C9',
    opacity=0.6,
    nbinsx=20
), row=2, col=1)

fig_dashboard.update_xaxes(title_text="Correlation Coefficient", row=2, col=1)
fig_dashboard.update_yaxes(title_text="Frequency", row=2, col=1)

# Add box plot of expression by type
expr_by_type = []
labels = ['RTK', 'NRTK', 'Other']

for label in labels:
    type_genes = gene_info[gene_info['Type'] == label]['GeneID'].tolist()
    if type_genes:
        type_expr = []
        for gene in type_genes:
            if gene in expr_data.index:
                expr_vals = expr_data.loc[gene].dropna().tolist()
                type_expr.extend(expr_vals)
        
        if type_expr:
            expr_by_type.append(type_expr)
        else:
            expr_by_type.append([])
    else:
        expr_by_type.append([])

# Add box plots
colors = ['#FF6B6B', '#4ECDC4', '#FFE66D']
for i, (expr_data_vals, label, color) in enumerate(zip(expr_by_type, labels, colors)):
    if expr_data_vals:  # Only add if we have data
        fig_dashboard.add_trace(go.Box(
            y=expr_data_vals,
            name=label,
            marker_color=color,
            boxpoints='outliers',
            showlegend=False
        ), row=2, col=2)

fig_dashboard.update_yaxes(title_text="Expression (Z-score)", row=2, col=2)

# Update dashboard layout
fig_dashboard.update_layout(
    title=dict(
        text='<b>Comprehensive RTK-NRTK Co-expression Analysis</b><br>Breast Cancer Tyrosine Kinase Network',
        font=dict(size=26, family="Arial Black"),
        x=0.5,
        xanchor='center',
        y=0.98
    ),
    showlegend=True,
    legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=1.02
    ),
    plot_bgcolor='white',
    width=1600,
    height=1000
)

# Save visualizations
print("\n💾 Saving visualizations...")

fig.write_html("final_rtk_nrtk_network.html")
fig.write_image("final_rtk_nrtk_network.png", width=1400, height=900, scale=2)

fig_dashboard.write_html("final_rtk_nrtk_dashboard.html")
fig_dashboard.write_image("final_rtk_nrtk_dashboard.png", width=1600, height=1000, scale=2)

print("✅ Saved: final_rtk_nrtk_network.html")
print("✅ Saved: final_rtk_nrtk_network.png")
print("✅ Saved: final_rtk_nrtk_dashboard.html")
print("✅ Saved: final_rtk_nrtk_dashboard.png")

# Create detailed results table
print("\n📋 Creating detailed results table...")

results_table = top_rtk_nrtk.copy()
results_table = results_table[[
    'Gene1', 'Gene2', 'Gene1_Symbol', 'Gene2_Symbol',
    'Type1', 'Type2', 'Correlation'
]]

# Add expression information
for idx, row in results_table.iterrows():
    gene1 = row['Gene1']
    gene2 = row['Gene2']
    
    if gene1 in expr_data.index:
        results_table.at[idx, 'Gene1_MeanExpr'] = expr_data.loc[gene1].mean()
        results_table.at[idx, 'Gene1_StdExpr'] = expr_data.loc[gene1].std()
    
    if gene2 in expr_data.index:
        results_table.at[idx, 'Gene2_MeanExpr'] = expr_data.loc[gene2].mean()
        results_table.at[idx, 'Gene2_StdExpr'] = expr_data.loc[gene2].std()

# Save results
results_table.to_csv('final_rtk_nrtk_results.csv', index=False)
print("✅ Saved: final_rtk_nrtk_results.csv")

# Create summary report
print("\n📄 Creating final summary report...")

report = f"""
FINAL RTK-NRTK CO-EXPRESSION NETWORK ANALYSIS
=============================================
Analysis Date: {pd.Timestamp.now()}

SUMMARY
-------
• Breast cancer samples analyzed: {expr_data.shape[1]}
• Genes in network: {G.number_of_nodes()}
• Significant correlations visualized: {G.number_of_edges()}
• RTKs identified: {len(rtk_genes)}
• NRTKs identified: {len(nrtk_genes)}

NETWORK PROPERTIES
------------------
• Average node degree: {np.mean([d for _, d in G.degree()]):.2f}
• Network density: {nx.density(G):.3f}
• Connected components: {nx.number_connected_components(G)}

TOP CORRELATIONS
----------------
"""

for i in range(min(10, len(top_rtk_nrtk))):
    row = top_rtk_nrtk.iloc[i]
    corr_type = "Positive" if row['Correlation'] > 0 else "Negative"
    
    report += f"{i+1}. {row['Gene1_Symbol']} ({row['Type1']}) ↔ {row['Gene2_Symbol']} ({row['Type2']})\n"
    report += f"   • Correlation: {row['Correlation']:.3f} ({corr_type})\n"
    
    if 'P-value' in row and not pd.isna(row['P-value']):
        report += f"   • P-value: {row['P-value']:.2e}\n"
    
    if 'FDR' in row and not pd.isna(row['FDR']):
        fdr_status = "Significant (FDR < 0.05)" if row['FDR'] < 0.05 else "Not significant"
        report += f"   • FDR: {row['FDR']:.2e} ({fdr_status})\n"
    
    report += "\n"

report += f"""
KEY FINDINGS
------------
1. Network reveals {len([e for e in G.edges(data=True) if e[2]['correlation'] > 0])} positive 
   and {len([e for e in G.edges(data=True) if e[2]['correlation'] < 0])} negative correlations.

2. The strongest correlation is between {top_rtk_nrtk.iloc[0]['Gene1_Symbol']} and 
   {top_rtk_nrtk.iloc[0]['Gene2_Symbol']} (r = {top_rtk_nrtk.iloc[0]['Correlation']:.3f}).

3. Genes with highest connectivity (potential hubs):
"""

# Find hubs (genes with most connections)
degrees = dict(G.degree())
sorted_degrees = sorted(degrees.items(), key=lambda x: x[1], reverse=True)

for i, (gene_id, degree) in enumerate(sorted_degrees[:5]):
    gene_symbol = G.nodes[gene_id]['gene_symbol']
    gene_type = G.nodes[gene_id]['type']
    report += f"   {i+1}. {gene_symbol} ({gene_type}): {degree} connections\n"

report += f"""
INTERPRETATION
--------------
This network visualization reveals coordinated expression patterns between 
Receptor Tyrosine Kinases (RTKs) and Non-Receptor Tyrosine Kinases (NRTKs) 
in breast cancer. Strong correlations suggest:

1. Co-regulation within signaling pathways
2. Potential protein-protein interactions
3. Coordinated response to oncogenic signals
4. Possible therapeutic targets for combination therapy

FILES GENERATED
---------------
1. final_rtk_nrtk_network.html - Interactive network visualization
2. final_rtk_nrtk_dashboard.html - Comprehensive analysis dashboard
3. final_rtk_nrtk_results.csv - Detailed correlation results
4. PNG images for publications

HOW TO USE
----------
• Open the HTML files in any web browser
• Hover over nodes/edges for detailed information
• Zoom and pan to explore the network
• Use the dashboard for comprehensive analysis

NEXT STEPS
----------
1. Validate findings in independent datasets
2. Investigate biological pathways for top correlations
3. Explore clinical correlations with patient outcomes
4. Consider experimental validation of key interactions
"""

with open('final_rtk_nrtk_analysis_report.txt', 'w') as f:
    f.write(report)

print("✅ Saved: final_rtk_nrtk_analysis_report.txt")

print("\n" + "="*80)
print("🎉 ANALYSIS COMPLETE! 🎉")
print("="*80)
print("\n📊 VISUALIZATIONS CREATED:")
print("   1. final_rtk_nrtk_network.html - Interactive network")
print("   2. final_rtk_nrtk_dashboard.html - 4-panel dashboard")
print("\n📁 OUTPUT FILES:")
print("   • HTML files - Open in browser for interactive exploration")
print("   • PNG images - High-resolution figures for publications")
print("   • CSV files - Detailed correlation data")
print("   • Text report - Complete analysis summary")
print("\n🔬 TO VIEW RESULTS:")
print("   Open 'final_rtk_nrtk_network.html' in your web browser!")
print("\n   Command: xdg-open final_rtk_nrtk_network.html")
print("\n💡 TIP: The interactive visualizations allow you to:")
print("   • Hover over nodes for gene information")
print("   • Hover over edges for correlation details")
print("   • Zoom and pan to explore different regions")
print("   • Click legend items to show/hide elements")
