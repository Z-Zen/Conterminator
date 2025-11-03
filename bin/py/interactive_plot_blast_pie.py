#!/usr/bin/env python3
"""
Create an interactive pie chart using Plotly
Usage: python plot_blast_interactive.py <blast_results.tsv> [sample_name]
"""

import sys
import plotly.graph_objects as go
from collections import Counter
import re
import os

def extract_genus(description):
    """Extract genus name from BLAST description"""
    match = re.search(r'^([A-Z][a-z]+)\s+[a-z]+', description)
    if match:
        return match.group(1)
    return description.split()[0] if description else "Unknown"

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_blast_interactive.py <blast_results.tsv> [sample_name]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    # Get sample name from argument or filename
    if len(sys.argv) >= 3:
        sample_name = sys.argv[2]
    else:
        # Extract from filename: "43_high_gc_blast_results.tsv" -> "43"
        sample_name = os.path.basename(input_file).split('_')[0]
    
    # Read and count genera
    genus_counts = Counter()
    
    with open(input_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 13:
                description = fields[12]
                genus = extract_genus(description)
                genus_counts[genus] += 1
    
    # Prepare data
    top_n = 10
    total_hits = sum(genus_counts.values())
    sorted_genera = genus_counts.most_common()
    top_genera = sorted_genera[:top_n]
    others_count = sum(count for _, count in sorted_genera[top_n:])
    
    labels = [genus for genus, _ in top_genera]
    values = [count for _, count in top_genera]
    
    if others_count > 0:
        labels.append('Others')
        values.append(others_count)
    
    # Create interactive pie chart
    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=values,
        hovertemplate='<b>%{label}</b><br>' +
                      'Hits: %{value:,}<br>' +
                      'Percentage: %{percent}<br>' +
                      '<extra></extra>',
        textinfo='label+percent',
        textposition='auto',
    )])
    
    fig.update_layout(
        title=f'Contamination Profile by Genus - Sample: {sample_name}<br>(Total BLAST hits: {total_hits:,})',
        title_font_size=16,
        showlegend=True,
        height=700,
        width=1000
    )
    
    # Save as HTML
    output_file = input_file.replace('.tsv', '_interactive.html')
    fig.write_html(output_file)
    print(f"Interactive chart saved to: {output_file}")
    print(f"Sample: {sample_name}")
    print(f"Open in browser: file://{output_file}")

if __name__ == "__main__":
    main()
