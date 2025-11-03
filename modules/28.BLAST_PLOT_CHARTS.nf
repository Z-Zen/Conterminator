process BLAST_PLOT_CHARTS {
    tag "${sample}_${db_name}"
    publishDir "${params.outdir}/Output/contamination_check/${sample}/${reads_type}/${db_name}/plots", mode: 'copy'
    cpus 1
    memory '4 GB'
    errorStrategy 'ignore'

    input:
    tuple val(sample), val(db_name), path(blast_tsv), path(blast_summary), val(reads_type)

    output:
    path "*.png", optional: true, emit: plots
    path "*.html", optional: true, emit: html

    when:
    params.run_contamination_check

    script:
    """
    set -euo pipefail
    
    # Check if BLAST results have hits
    if [ ! -s ${blast_tsv} ]; then
        echo "No BLAST hits found, skipping plots"
        exit 0
    fi
    
    # Count number of hits
    HITS=\$(wc -l < ${blast_tsv})
    
    if [ "\$HITS" -lt 1 ]; then
        echo "No BLAST hits found, skipping plots"
        exit 0
    fi
    
    echo "Generating plots for ${sample} (${db_name}): \$HITS hits"
    
    # Generate matplotlib pie chart (overview)
    python3 ${params.blast_plot} ${blast_tsv} ${sample}_${db_name}
    
    # Generate interactive plotly chart
    python3 ${params.blast_plot_interac} ${blast_tsv} ${sample}_${db_name}
    
    # Auto-detect top genus and create breakdown
    TOP_GENUS=\$(python3 -c "
import sys
from collections import Counter
import re

def extract_genus(desc):
    clean = re.sub(r'^PREDICTED:\\s+', '', desc)
    clean = re.sub(r'^\\S+:\\s+', '', clean)
    if clean.lower().startswith('mouse'):
        return 'Mus'
    match = re.search(r'^([A-Z][a-z]+)\\s+[a-z]+', clean)
    return match.group(1) if match else clean.split()[0]

genus_counts = Counter()
with open('${blast_tsv}') as f:
    for line in f:
        fields = line.strip().split('\\t')
        if len(fields) >= 13:
            genus_counts[extract_genus(fields[12])] += 1

if genus_counts:
    print(genus_counts.most_common(1)[0][0])
" || echo "")
    
    # Generate genus breakdown if top genus found and has enough hits
    if [ -n "\$TOP_GENUS" ] && [ "\$HITS" -ge 10 ]; then
        echo "Generating breakdown for genus: \$TOP_GENUS"
        python3 ${params.blast_plot} ${blast_tsv} ${sample}_${db_name} --genus "\$TOP_GENUS" || true
    fi
    
    echo "Plot generation complete"
    ls -lh *.png *.html 2>/dev/null || echo "No plot files found"
    """
}