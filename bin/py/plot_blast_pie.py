#!/usr/bin/env python3
"""
Create a pie chart of major species/genera from BLAST contamination results
Usage:
  python plot_blast_pie.py <blast_results.tsv> [sample_name]
  python plot_blast_pie.py <blast_results.tsv> [sample_name] --genus <genus_name>
"""

import sys
import matplotlib.pyplot as plt
from collections import Counter
import re
import os
import argparse


def extract_genus(description):
    """Extract and normalize genus name from BLAST description."""
    clean_desc = re.sub(r'^PREDICTED:\s+', '', description)
    clean_desc = re.sub(r'^\S+:\s+', '', clean_desc)

    if clean_desc.lower().startswith('mouse'):
        return 'Mus'

    match = re.search(r'^([A-Z][a-z]+)\s+[a-z]+', clean_desc)
    if match:
        genus = match.group(1)
        if genus.lower() == 'mouse':
            return 'Mus'
        return genus

    return clean_desc.split()[0] if clean_desc else "Unknown"

def extract_species(description):
    """Extract species/subspecies name from BLAST description."""
    clean_desc = re.sub(r'^PREDICTED:\s+', '', description)
    clean_desc = re.sub(r'^\S+:\s+', '', clean_desc)
    
    # Match "Genus species" or "Genus species subspecies"
    match = re.search(r'^([A-Z][a-z]+(?:\s+[a-z]+){1,2})', clean_desc)
    if match:
        return match.group(1)
        
    words = clean_desc.split()
    if len(words) >= 2:
        return f"{words[0]} {words[1]}"
    return words[0] if words else "Unknown"

def normalize_description(desc):
    """Normalize descriptions to group similar entries, preserving all details."""
    species = extract_species(desc)
    base_label = species

    # More robust regex to find strain or isolate, handling different formats
    strain_match = re.search(r'(?:strain|isolate)[\s:]+([^\s,]+)', desc, re.IGNORECASE)
    if strain_match:
        strain_name = strain_match.group(1).strip()
        # Avoid adding generic words as the strain name
        if strain_name.lower() not in ['complete', 'mitochondrion', 'genome']:
            base_label += f" strain {strain_name}"

    # Append context like (mitochondrion) or (genome) instead of replacing the label
    context = ""
    if 'mitochondrion' in desc.lower() or 'mitochondrial' in desc.lower():
        context = " (mitochondrion)"
    elif 'chromosome' in desc.lower() or 'scaffold' in desc.lower() or 'unlocalized' in desc.lower() or 'unplaced' in desc.lower():
        context = " (genome)"
        
    return base_label + context

def plot_genus_overview(genus_counts, sample_name, input_file, total_hits):
    """Plot pie chart of all genera, with percentages in the legend."""
    top_n = 10
    sorted_genera = genus_counts.most_common()
    top_genera = sorted_genera[:top_n]
    others_count = sum(count for _, count in sorted_genera[top_n:])

    labels = [genus for genus, _ in top_genera]
    sizes = [count for _, count in top_genera]

    if others_count > 0:
        labels.append('Others')
        sizes.append(others_count)

    percentages = [(s / total_hits) * 100 for s in sizes]
    display_labels = [label if pct >= 2 else '' for label, pct in zip(labels, percentages)]

    fig, ax = plt.subplots(figsize=(12, 8))
    colors = plt.cm.Set3(range(len(labels)))

    wedges, texts, autotexts = ax.pie(
        sizes, labels=display_labels, colors=colors,
        autopct=lambda p: f'{p:.1f}%' if p >= 2 else '',
        startangle=90, textprops={'fontsize': 10}
    )

    for autotext in autotexts:
        autotext.set_color('black')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(9)

    ax.set_title(f'Contamination Profile by Genus - Sample: {sample_name}\n(Total BLAST hits: {total_hits:,})',
                 fontsize=14, fontweight='bold', pad=20)

    legend_labels = [f'{label}: {size:,} hits ({pct:.1f}%)' for label, size, pct in zip(labels, sizes, percentages)]
    ax.legend(wedges, legend_labels, loc='center left', bbox_to_anchor=(1, 0, 0.5, 1), fontsize=9)
    plt.tight_layout()

    output_file = input_file.replace('.tsv', '_pie_chart.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')

    print(f"Pie chart saved to: {output_file}")
    print(f"\nSample: {sample_name}")
    print(f"Top {top_n} Genera:")
    print("-" * 50)
    for genus, count in top_genera:
        pct = (count/total_hits)*100
        print(f"{genus:20s}: {count:6,} hits ({pct:5.1f}%)")
    if others_count > 0:
        pct = (others_count/total_hits)*100
        print(f"{'Others':20s}: {others_count:6,} hits ({pct:5.1f}%)")

def plot_genus_breakdown(descriptions, genus_name, sample_name, input_file, total_hits):
    """Plot pie chart of species/strains within a specific genus"""
    desc_counts = Counter(normalize_description(d) for d in descriptions if extract_genus(d) == genus_name)

    if not desc_counts:
        print(f"Error: No hits found for genus '{genus_name}'")
        return

    genus_total = sum(desc_counts.values())
    top_n = 15
    sorted_descs = desc_counts.most_common()
    top_descs = sorted_descs[:top_n]
    others_count = sum(count for _, count in sorted_descs[top_n:])

    labels = [desc for desc, _ in top_descs]
    sizes = [count for _, count in top_descs]

    if others_count > 0:
        labels.append('Others')
        sizes.append(others_count)

    percentages = [(s / genus_total) * 100 for s in sizes]
    display_labels = [label if pct >= 4 else '' for label, pct in zip(labels, percentages)]

    fig, ax = plt.subplots(figsize=(14, 10))
    colors = plt.cm.Paired(range(len(labels)))

    wedges, texts, autotexts = ax.pie(
        sizes, labels=display_labels, colors=colors,
        autopct=lambda p: f'{p:.1f}%' if p >= 4 else '',
        startangle=90, textprops={'fontsize': 8}
    )

    for autotext in autotexts:
        autotext.set_color('black')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(8)

    genus_pct_total = (genus_total / total_hits) * 100
    ax.set_title(f'Breakdown for Genus: {genus_name} - Sample: {sample_name}\n' +
                 f'({genus_name} hits: {genus_total:,} of {total_hits:,} total = {genus_pct_total:.1f}%)',
                 fontsize=14, fontweight='bold', pad=20)

    legend_labels = [f'{label}: {size:,} hits ({pct:.1f}%)' for label, size, pct in zip(labels, sizes, percentages)]
    ax.legend(wedges, legend_labels, loc='center left', bbox_to_anchor=(1, 0, 0.5, 1), fontsize=9)
    plt.tight_layout()

    output_file = input_file.replace('.tsv', f'_{genus_name}_breakdown.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')

    print(f"\nGenus breakdown chart saved to: {output_file}")
    print(f"\nBreakdown for {genus_name}:")
    print("-" * 70)
    for (desc, count), pct in zip(top_descs, percentages):
        print(f"{desc[:60]:60s}: {count:6,} hits ({pct:5.1f}%)")
    if others_count > 0:
        pct = (others_count / genus_total) * 100
        print(f"{'Others':60s}: {others_count:6,} hits ({pct:5.1f}%)")

def main():
    parser = argparse.ArgumentParser(
        description='Create pie charts from BLAST contamination results',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('input_file', help='BLAST results TSV file')
    parser.add_argument('sample_name', nargs='?', help='Sample name (optional, derived from filename if not provided)')
    parser.add_argument('--genus', '-g', help='Generate species breakdown for specific genus')
    args = parser.parse_args()

    if args.sample_name:
        sample_name = args.sample_name
    else:
        basename = os.path.basename(args.input_file)
        # A more robust way to get a sample name from filenames like 'blast_results_43_new.tsv'
        sample_name = os.path.splitext(basename)[0].replace('blast_results_', '').replace('_new', '')

    genus_counts = Counter()
    descriptions = []

    try:
        with open(args.input_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                # stitle is the last column, a format "6" output has 12 columns by default
                if len(fields) >= 13: 
                    description = fields[12]
                    descriptions.append(description)
                    genus = extract_genus(description)
                    genus_counts[genus] += 1
    except FileNotFoundError:
        print(f"Error: Input file not found at '{args.input_file}'")
        sys.exit(1)

    total_hits = sum(genus_counts.values())

    if not total_hits:
        print("Error: No valid BLAST hits found in the input file.")
        sys.exit(1)

    if args.genus:
        plot_genus_breakdown(descriptions, args.genus, sample_name, args.input_file, total_hits)
    else:
        plot_genus_overview(genus_counts, sample_name, args.input_file, total_hits)

if __name__ == "__main__":
    main()
