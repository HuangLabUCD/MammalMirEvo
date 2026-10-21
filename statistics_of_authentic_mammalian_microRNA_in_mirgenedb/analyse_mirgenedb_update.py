#!/usr/bin/env python3
"""
Analyze mammalian MirGeneDB miRNA precursor sequences
Generates statistics and distribution plots for key features
"""

import RNA
import Bio.SeqIO
import sys
import re
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
from pathlib import Path

# Set style for better looking plots
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

def analyze_sequence(record):
    """
    Fold and analyze a single sequence using BOTH MFE and centroid structures
    """
    seq = str(record.seq).upper()
    seq_length = len(seq)
    
    # Initialize default values
    mfe_struct = ""
    centroid_struct = ""
    mfe = 0
    struct_distance = 0.0
    bp_distance = 0.0
    
    try:
        md = RNA.md()
        md.temperature = 37.0
        md.dangles = 2
        md.noLP = True
        fc = RNA.fold_compound(seq, md)
        mfe_struct, mfe = fc.mfe()
        
        if mfe_struct is None:
            mfe_struct = ""
        elif not isinstance(mfe_struct, str):
            mfe_struct = str(mfe_struct) if mfe_struct else ""
        
        try:
            fc.pf()
            centroid_result = fc.centroid()
            
            if centroid_result is None:
                centroid_struct = mfe_struct
            elif isinstance(centroid_result, str):
                centroid_struct = centroid_result
            elif isinstance(centroid_result, (list, tuple)) and len(centroid_result) > 0:
                if isinstance(centroid_result[0], str):
                    centroid_struct = centroid_result[0]
                else:
                    centroid_struct = str(centroid_result[0]) if centroid_result[0] else mfe_struct
            else:
                centroid_struct = mfe_struct
        except:
            centroid_struct = mfe_struct
        
        if (mfe_struct and centroid_struct and 
            isinstance(mfe_struct, str) and isinstance(centroid_struct, str) and
            len(mfe_struct) == len(centroid_struct) and len(mfe_struct) > 0):
            struct_distance = calculate_structure_distance(mfe_struct, centroid_struct)
            bp_distance = calculate_base_pair_distance(mfe_struct, centroid_struct)
        
    except Exception as e:
        print(f"  Warning: Folding failed for {record.id}: {e}")
        mfe_struct = ""
        centroid_struct = ""
    
    # Analyze structures
    mfe_analysis = analyze_structure(mfe_struct, seq, mfe)
    centroid_analysis = analyze_structure(centroid_struct, seq, mfe) if centroid_struct else mfe_analysis.copy()
    
    return {
        'id': record.id,
        'sequence': seq,
        'length': seq_length,
        'loop_count': mfe_analysis['loop_count'],
        'terminal_loop_length': mfe_analysis['terminal_loop_length'],
        'pairing_percent': mfe_analysis['pairing_percent'],
        'mfei': mfe_analysis['mfei'],
        'stem_segments': mfe_analysis['stem_segments'],
        'structure_distance': struct_distance,
        'base_pair_distance': bp_distance,
        'mfe': mfe,
        'centroid_loop_count': centroid_analysis['loop_count'],
        'centroid_terminal_loop_length': centroid_analysis['terminal_loop_length'],
    }

def analyze_structure(structure, seq, mfe=None):
    """Extract features from a structure"""
    if not structure or not isinstance(structure, str) or structure == "":
        return {
            'loop_count': 0,
            'terminal_loop_length': 0,
            'pairing_percent': 0,
            'mfei': 0,
            'stem_segments': 0
        }
    
    seq_length = len(seq)
    
    # Loop analysis
    loop_count = count_loops_simple(structure)
    terminal_loop_length = get_terminal_loop_length(structure)
    terminal_loop_sequence = get_terminal_loop_sequence(seq, structure)
    
    # Stem analysis
    paired_bases = structure.count('(') + structure.count(')')
    pairing_percent = (paired_bases / seq_length) * 100 if seq_length > 0 else 0
    
    # Count stem segments
    in_stem = False
    stem_segments = 0
    for char in structure:
        if char == '(' and not in_stem:
            stem_segments += 1
            in_stem = True
        elif char == ')':
            in_stem = False
        elif char == '.':
            in_stem = False
    
    # Thermodynamic features
    gc_content = calculate_gc_content(seq)
    mfei = calculate_mfei(mfe, seq_length, gc_content) if mfe and gc_content > 0 else 0
    
    return {
        'loop_count': loop_count,
        'terminal_loop_length': terminal_loop_length,
        'pairing_percent': pairing_percent,
        'mfei': mfei,
        'stem_segments': stem_segments
    }

def calculate_structure_distance(struct1, struct2):
    if not struct1 or not struct2 or len(struct1) != len(struct2):
        return 1.0
    differences = sum(1 for i in range(len(struct1)) if struct1[i] != struct2[i])
    return differences / len(struct1)

def calculate_base_pair_distance(struct1, struct2):
    def get_base_pairs(structure):
        pairs = set()
        stack = []
        for i, char in enumerate(structure):
            if char == '(':
                stack.append(i)
            elif char == ')' and stack:
                j = stack.pop()
                pairs.add((j, i))
        return pairs
    
    pairs1 = get_base_pairs(struct1)
    pairs2 = get_base_pairs(struct2)
    if len(pairs1) == 0 and len(pairs2) == 0:
        return 0.0
    diff = len(pairs1.symmetric_difference(pairs2))
    total = max(len(pairs1), len(pairs2))
    return diff / total if total > 0 else 0

def count_loops_simple(structure):
    if not structure or not isinstance(structure, str):
        return 0
    
    loop_count = 0
    for i in range(len(structure) - 1):
        if structure[i] == ')' and structure[i+1] == '(':
            loop_count += 1
    
    try:
        hairpins = re.findall(r'\(+\.+\)+', structure)
        loop_count += len(hairpins)
    except:
        pass
    
    if loop_count == 0 and ('(' in structure or ')' in structure):
        loop_count = 1
    
    return loop_count

def get_terminal_loop_length(structure):
    if not structure or not isinstance(structure, str):
        return 0
    try:
        matches = re.findall(r'\(+\.+\)+', structure)
        if matches:
            loop_region = matches[-1]
            return loop_region.count('.')
    except:
        pass
    return 0

def get_terminal_loop_sequence(seq, structure):
    if not structure or not seq or not isinstance(structure, str):
        return ""
    try:
        match = re.search(r'\(+\.+\)+', structure)
        if match:
            start = match.start()
            end = match.end()
            loop_start = start + structure[start:end].find('.')
            loop_end = loop_start + structure[start:end].count('.')
            if loop_end <= len(seq):
                return seq[loop_start:loop_end]
    except:
        pass
    return ""

def calculate_gc_content(seq):
    if not seq:
        return 0
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100

def calculate_mfei(mfe, length, gc_content):
    if gc_content == 0 or mfe is None:
        return 0
    mfe_per_nt = abs(mfe) / length
    return mfe_per_nt / (gc_content / 100)

def generate_stats_and_plots(data, output_dir):
    """
    Generate statistics and plots for all features
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # REMOVED gc_content and stem_segments from features
    features = {
        'loop_count': {'name': 'Loop Count', 'xlabel': 'Number of loops', 'dtype': 'int'},
        'terminal_loop_length': {'name': 'Terminal Loop Length', 'xlabel': 'Length (nt)', 'dtype': 'int'},
        'pairing_percent': {'name': 'Pairing Percentage', 'xlabel': 'Percentage (%)', 'dtype': 'float'},
        'mfei': {'name': 'MFEI (Minimum Free Energy Index)', 'xlabel': 'MFEI value', 'dtype': 'float'},
        'length': {'name': 'Sequence Length', 'xlabel': 'Length (nt)', 'dtype': 'int'},
        'structure_distance': {'name': 'Structure Distance (MFE vs Centroid)', 'xlabel': 'Distance', 'dtype': 'float'}
    }
    
    # Calculate statistics
    stats_file = os.path.join(output_dir, "statistics.txt")
    with open(stats_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("MAMMALIAN MirGeneDB miRNA PRECURSOR ANALYSIS\n")
        f.write("="*80 + "\n\n")
        f.write(f"Total sequences analyzed: {len(data)}\n\n")
        
        f.write("-"*80 + "\n")
        f.write("FEATURE STATISTICS\n")
        f.write("-"*80 + "\n\n")
        
        for feature_key, feature_info in features.items():
            values = [d[feature_key] for d in data if d[feature_key] is not None]
            values = np.array(values)
            
            f.write(f"{feature_info['name']}:\n")
            f.write(f"  Count: {len(values)}\n")
            f.write(f"  Mean: {np.mean(values):.2f}\n")
            f.write(f"  Median: {np.median(values):.2f}\n")
            f.write(f"  Std Dev: {np.std(values):.2f}\n")
            f.write(f"  Min: {np.min(values):.2f}\n")
            f.write(f"  Max: {np.max(values):.2f}\n")
            
            # Percentiles
            for p in [1, 5, 10, 25, 50, 75, 90, 95, 99]:
                f.write(f"  {p}th percentile: {np.percentile(values, p):.2f}\n")
            
            f.write("\n")
    
    # Generate plots
    print("\nGenerating plots...")
    
    # Create figure with subplots (2 rows x 3 columns for 6 features)
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    for idx, (feature_key, feature_info) in enumerate(features.items()):
        ax = axes[idx]
        values = [d[feature_key] for d in data if d[feature_key] is not None]
        
        if feature_info['dtype'] == 'int':
            # For integer features, use bar plot
            unique_vals = sorted(set(values))
            counts = [values.count(v) for v in unique_vals]
            
            # Special handling for sequence length to avoid overlapping x-labels
            if feature_key == 'length':
                # Calculate optimal number of ticks based on data range
                min_val = min(unique_vals)
                max_val = max(unique_vals)
                val_range = max_val - min_val
                
                # Determine tick step (aim for ~10 ticks or less)
                if val_range <= 20:
                    step = 2
                elif val_range <= 50:
                    step = 5
                elif val_range <= 100:
                    step = 10
                else:
                    step = 20
                
                # Generate tick positions
                tick_positions = list(range(min_val, max_val + 1, step))
                
                # Create bar plot
                ax.bar(unique_vals, counts, edgecolor='black', alpha=0.7, color='skyblue', width=0.8)
                ax.set_xticks(tick_positions)
                # REMOVED rotation for x labels
                ax.set_xticklabels(tick_positions, fontsize=8)
                
                # Optionally, add a note if step > 1
                if step > 2:
                    ax.set_xlabel(f"{feature_info['xlabel']} (tick every {step} nt)", fontsize=10)
                else:
                    ax.set_xlabel(feature_info['xlabel'], fontsize=10)
            else:
                ax.bar(unique_vals, counts, edgecolor='black', alpha=0.7, color='skyblue')
                ax.set_xticks(unique_vals)
                ax.set_xlabel(feature_info['xlabel'], fontsize=10)
        else:
            # For continuous features, use histogram
            ax.hist(values, bins=30, edgecolor='black', alpha=0.7, color='skyblue')
            ax.set_xlabel(feature_info['xlabel'], fontsize=10)
        
        # Only add y-label "Frequency" for first plot in each row (indices 0 and 3)
        if idx == 0 or idx == 3:
            ax.set_ylabel('Frequency', fontsize=10)
        else:
            ax.set_ylabel('')  # Empty y-label for other plots
        
        # REMOVED (n=xxxx) from title
        ax.set_title(f"{feature_info['name']}", fontsize=12)
        ax.grid(False)
    
    # Hide any unused subplots (if less than 6 features)
    for idx in range(len(features), len(axes)):
        axes[idx].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "all_features_distributions.png"), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Individual plots for each feature (higher resolution)
    for feature_key, feature_info in features.items():
        fig, ax = plt.subplots(figsize=(10, 6))
        
        values = [d[feature_key] for d in data if d[feature_key] is not None]
        
        if feature_info['dtype'] == 'int':
            unique_vals = sorted(set(values))
            counts = [values.count(v) for v in unique_vals]
            
            # Special handling for sequence length to avoid overlapping x-labels
            if feature_key == 'length':
                # Calculate optimal number of ticks based on data range
                min_val = min(unique_vals)
                max_val = max(unique_vals)
                val_range = max_val - min_val
                
                # Determine tick step (aim for ~10-15 ticks maximum)
                if val_range <= 20:
                    step = 2
                elif val_range <= 50:
                    step = 5
                elif val_range <= 100:
                    step = 10
                elif val_range <= 200:
                    step = 20
                else:
                    step = 50
                
                # Generate tick positions
                tick_positions = list(range(min_val, max_val + 1, step))
                
                # Create bar plot
                ax.bar(unique_vals, counts, edgecolor='black', alpha=0.7, color='skyblue', width=0.8)
                ax.set_xticks(tick_positions)
                # REMOVED rotation for x labels
                ax.set_xticklabels(tick_positions, fontsize=10)
                
                # Add informative x-label
                if step > 2:
                    ax.set_xlabel(f"{feature_info['xlabel']} (tick every {step} nt)", fontsize=12)
                else:
                    ax.set_xlabel(feature_info['xlabel'], fontsize=12)
                
                # Optionally, add minor ticks for all values
                ax.set_xticks(unique_vals, minor=True)
                ax.tick_params(axis='x', which='minor', length=2, labelsize=0)
            else:
                ax.bar(unique_vals, counts, edgecolor='black', alpha=0.7, color='skyblue')
                ax.set_xticks(unique_vals)
                ax.set_xlabel(feature_info['xlabel'], fontsize=12)
        else:
            ax.hist(values, bins=30, edgecolor='black', alpha=0.7, color='skyblue')
            ax.set_xlabel(feature_info['xlabel'], fontsize=12)
        
        ax.set_ylabel('Frequency', fontsize=12)
        # REMOVED (n=xxxx) from title
        ax.set_title(f"{feature_info['name']}", fontsize=14)
        ax.grid(False)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{feature_key}_distribution.png"), dpi=150, bbox_inches='tight')
        plt.close()
    
    # Generate box plot summary (updated with remaining features)
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Prepare data for box plots (normalize for visualization)
    plot_data = []
    labels = []
    for feature_key, feature_info in features.items():
        values = [d[feature_key] for d in data if d[feature_key] is not None]
        # Normalize to 0-1 range for comparison
        min_val = np.min(values)
        max_val = np.max(values)
        if max_val > min_val:
            normalized = [(v - min_val) / (max_val - min_val) for v in values]
        else:
            normalized = values
        plot_data.append(normalized)
        labels.append(feature_info['name'])
    
    bp = ax.boxplot(plot_data, labels=labels, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('skyblue')
        patch.set_alpha(0.7)
    
    ax.set_ylabel('Normalized Value (0-1 scale)', fontsize=12)
    ax.set_title('Feature Distribution Comparison (Normalized)', fontsize=14)
    ax.tick_params(axis='x', rotation=45)
    ax.grid(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "features_boxplot.png"), dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"\nStatistics and plots saved to: {output_dir}")

def generate_correlation_matrix(data, output_dir):
    """
    Generate correlation matrix between features
    """
    # REMOVED gc_content and stem_segments from features
    features = ['loop_count', 'terminal_loop_length', 'pairing_percent', 
                'mfei', 'length', 'structure_distance']
    
    # Create correlation matrix
    corr_data = []
    for d in data:
        row = [d[f] for f in features]
        corr_data.append(row)
    
    corr_matrix = np.corrcoef(np.array(corr_data).T)
    
    # Plot correlation matrix
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(corr_matrix, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    
    # Add colorbar
    plt.colorbar(im, ax=ax)
    
    # Add labels
    ax.set_xticks(range(len(features)))
    ax.set_yticks(range(len(features)))
    ax.set_xticklabels(features, rotation=45, ha='right')
    ax.set_yticklabels(features)
    
    # Add correlation values
    for i in range(len(features)):
        for j in range(len(features)):
            ax.text(j, i, f'{corr_matrix[i, j]:.2f}',
                    ha="center", va="center", color="black", fontsize=8)
    
    ax.set_title('Feature Correlation Matrix', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "correlation_matrix.png"), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Save correlation matrix as text
    corr_file = os.path.join(output_dir, "correlation_matrix.txt")
    with open(corr_file, 'w') as f:
        f.write("Feature Correlation Matrix\n")
        f.write("="*60 + "\n\n")
        f.write("\t" + "\t".join(features) + "\n")
        for i, feat_i in enumerate(features):
            f.write(f"{feat_i}\t")
            for j, feat_j in enumerate(features):
                f.write(f"{corr_matrix[i, j]:.3f}\t")
            f.write("\n")
    
    print(f"Correlation matrix saved to: {corr_file}")

def generate_loop_count_detailed(data, output_dir):
    """
    Generate detailed analysis of loop count patterns
    """
    loop_counts = [d['loop_count'] for d in data]
    centroid_loops = [d['centroid_loop_count'] for d in data]
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # MFE loop count distribution
    unique_loops = sorted(set(loop_counts))
    counts = [loop_counts.count(l) for l in unique_loops]
    axes[0].bar(unique_loops, counts, edgecolor='black', alpha=0.7, color='skyblue')
    axes[0].set_xlabel('Number of loops (MFE structure)', fontsize=12)
    axes[0].set_ylabel('Frequency', fontsize=12)
    axes[0].set_title('MFE Loop Count Distribution', fontsize=12)
    axes[0].grid(False)
    
    # Add percentage labels
    for i, (l, c) in enumerate(zip(unique_loops, counts)):
        pct = c / len(loop_counts) * 100
        axes[0].text(l, c, f'{pct:.1f}%', ha='center', va='bottom')
    
    # Centroid vs MFE loop count agreement
    agreement = [1 if l == c else 0 for l, c in zip(loop_counts, centroid_loops)]
    agree_count = sum(agreement)
    disagree_count = len(agreement) - agree_count
    
    axes[1].bar(['Agree', 'Disagree'], [agree_count, disagree_count], 
                edgecolor='black', alpha=0.7, color=['skyblue', 'lightcoral'])
    axes[1].set_ylabel('Frequency', fontsize=12)
    axes[1].set_title(f'MFE vs Centroid Loop Count Agreement\n(Agree: {agree_count/len(agreement)*100:.1f}%)', fontsize=12)
    axes[1].grid(False)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "loop_count_detailed.png"), dpi=150, bbox_inches='tight')
    plt.close()
    
    # Save agreement statistics
    agree_file = os.path.join(output_dir, "loop_agreement.txt")
    with open(agree_file, 'w') as f:
        f.write("Loop Count Agreement Analysis\n")
        f.write("="*60 + "\n\n")
        f.write(f"Total sequences: {len(loop_counts)}\n")
        f.write(f"MFE and Centroid agree: {agree_count} ({agree_count/len(loop_counts)*100:.1f}%)\n")
        f.write(f"MFE and Centroid disagree: {disagree_count} ({disagree_count/len(loop_counts)*100:.1f}%)\n\n")
        
        # Show disagreement patterns
        f.write("Disagreement patterns:\n")
        for l, c in zip(loop_counts, centroid_loops):
            if l != c:
                f.write(f"  MFE: {l} loop(s) → Centroid: {c} loop(s)\n")

def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_mirgenedb.py <input_fasta> [output_dir]")
        print("\nExample:")
        print("  python analyze_mirgenedb.py mammalian_mirnas.fa")
        print("  python analyze_mirgenedb.py mammalian_mirnas.fa my_analysis/")
        sys.exit(1)
    
    input_file = sys.argv[1]
    
    if len(sys.argv) > 2:
        output_dir = sys.argv[2]
    else:
        output_dir = "mirgenedb_analysis"
    
    print(f"\n{'='*80}")
    print("MAMMALIAN MirGeneDB miRNA PRECURSOR ANALYSIS")
    print(f"{'='*80}")
    print(f"Input file: {input_file}")
    print(f"Output directory: {output_dir}")
    print(f"{'='*80}\n")
    
    # Read and analyze sequences
    print("Analyzing sequences (this may take a while)...")
    data = []
    
    for i, record in enumerate(Bio.SeqIO.parse(input_file, "fasta"), 1):
        if i % 100 == 0:
            print(f"  Processed {i} sequences...")
        
        analysis = analyze_sequence(record)
        data.append(analysis)
    
    print(f"\n✓ Analyzed {len(data)} sequences")
    
    # Generate statistics and plots
    generate_stats_and_plots(data, output_dir)
    generate_correlation_matrix(data, output_dir)
    generate_loop_count_detailed(data, output_dir)
    
    # Print quick summary to screen
    print(f"\n{'='*80}")
    print("QUICK SUMMARY")
    print(f"{'='*80}")
    
    loop_counts = [d['loop_count'] for d in data]
    print(f"\nLoop count distribution:")
    for l in sorted(set(loop_counts)):
        count = loop_counts.count(l)
        pct = count / len(loop_counts) * 100
        print(f"  {l} loop(s): {count} ({pct:.1f}%)")
    
    print(f"\nTerminal loop length:")
    tl = [d['terminal_loop_length'] for d in data if d['terminal_loop_length'] > 0]
    print(f"  Range: {min(tl)} - {max(tl)} nt")
    print(f"  Mean: {np.mean(tl):.1f} nt")
    print(f"  Median: {np.median(tl):.1f} nt")
    
    print(f"\nPairing percentage:")
    pp = [d['pairing_percent'] for d in data]
    print(f"  Range: {min(pp):.1f} - {max(pp):.1f}%")
    print(f"  Mean: {np.mean(pp):.1f}%")
    print(f"  Median: {np.median(pp):.1f}%")
    
    print(f"\nMFEI:")
    mfei = [d['mfei'] for d in data if d['mfei'] > 0]
    print(f"  Range: {min(mfei):.2f} - {max(mfei):.2f}")
    print(f"  Mean: {np.mean(mfei):.2f}")
    print(f"  Median: {np.median(mfei):.2f}")
    
    print(f"\nSequence length:")
    lengths = [d['length'] for d in data]
    print(f"  Range: {min(lengths)} - {max(lengths)} nt")
    print(f"  Mean: {np.mean(lengths):.1f} nt")
    print(f"  Median: {np.median(lengths):.1f} nt")
    
    print(f"\nStructure distance (MFE vs Centroid):")
    sd = [d['structure_distance'] for d in data]
    print(f"  Range: {min(sd):.3f} - {max(sd):.3f}")
    print(f"  Mean: {np.mean(sd):.3f}")
    print(f"  Median: {np.median(sd):.3f}")
    
    print(f"\n{'='*80}")
    print("✓ Analysis complete!")
    print(f"  Results saved to: {output_dir}")
    print("  Files generated:")
    print("    - statistics.txt (full statistics)")
    print("    - all_features_distributions.png (summary plot)")
    print("    - [feature]_distribution.png (individual plots)")
    print("    - correlation_matrix.png (feature correlations)")
    print("    - loop_count_detailed.png (loop analysis)")
    print(f"{'='*80}\n")

if __name__ == "__main__":
    import Bio.SeqIO
    main()