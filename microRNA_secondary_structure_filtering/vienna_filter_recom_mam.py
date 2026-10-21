#!/usr/bin/env python3
"""
Batch filter for miRNA precursors - MAMMALIAN OPTIMIZED VERSION
Based on empirical thresholds from 7,556 mammalian MirGeneDB sequences
Matches online RNAfold command: RNAfold -p -d2 --noLP
Compatible with ViennaRNA 2.7.2

Thresholds:
  Strict mode:  loop_count=1, term_loop=3-15, pairing≥55%, MFEI≥0.60, length=54-86
  Permissive mode: loop_count≤2, term_loop≥2, pairing≥50%, MFEI≥0.50, length=50-86
"""

import RNA
import Bio.SeqIO
import sys
import re
import os
import glob
from collections import defaultdict
from pathlib import Path

# ============================================================================
# THRESHOLDS BASED ON MAMMALIAN MirGeneDB ANALYSIS
# ============================================================================

STRICT_THRESHOLDS = {
    # Hard binary rules
    'loop_count': 1,                    # 99.4% of mammalian miRNAs
    'max_centroid_loops': 2,            # 99.6% have ≤2
    'require_centroid_match': True,     # MFE and centroid loop counts must match
    
    # Range-based filters
    'terminal_loop_min': 3,
    'terminal_loop_max': 15,
    'pairing_percent': 55,              # Minimum 55%
    'mfei': 0.60,                       # Minimum 0.60
    'length_min': 54,
    'length_max': 86,
    'structure_distance': 0.15,         # Maximum 0.15
    'max_loop_ratio': 3.0,              # Centroid/MFE terminal loop ratio
}

PERMISSIVE_THRESHOLDS = {
    # Hard binary rules (relaxed)
    'loop_count': 2,                    # Allow ≤2 loops (captures let-7)
    'max_centroid_loops': 2,
    'require_centroid_match': True,
    
    # Range-based filters (relaxed)
    'terminal_loop_min': 2,             # Allow let-7's 2-nt loop
    'terminal_loop_max': 22,            # Full range from data
    'pairing_percent': 50,              # Minimum 50%
    'mfei': 0.50,                       # Minimum 0.50
    'length_min': 50,
    'length_max': 86,
    'structure_distance': 0.25,         # More tolerant
    'max_loop_ratio': 5.0,              # Allow let-7's large ratio
}

# ============================================================================
# CORE FUNCTIONS
# ============================================================================

def analyze_sequence(record):
    """
    Fold and analyze a single sequence using BOTH MFE and centroid structures
    Matches online RNAfold -p -d2 --noLP parameters for ViennaRNA 2.7.2
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
        try:
            fc = RNA.fold_compound(seq)
            mfe_struct, mfe = fc.mfe()
            if mfe_struct is None:
                mfe_struct = ""
            elif not isinstance(mfe_struct, str):
                mfe_struct = str(mfe_struct) if mfe_struct else ""
            
            try:
                fc.pf()
                centroid_result = fc.centroid()
                if isinstance(centroid_result, str):
                    centroid_struct = centroid_result
                elif isinstance(centroid_result, (list, tuple)) and len(centroid_result) > 0:
                    centroid_struct = str(centroid_result[0]) if centroid_result[0] else mfe_struct
                else:
                    centroid_struct = mfe_struct
            except:
                centroid_struct = mfe_struct
                
            if mfe_struct and centroid_struct and len(mfe_struct) == len(centroid_struct):
                struct_distance = calculate_structure_distance(mfe_struct, centroid_struct)
                bp_distance = calculate_base_pair_distance(mfe_struct, centroid_struct)
        except Exception as e2:
            print(f"  Warning: Folding failed for {record.id}: {e2}")
            mfe_struct = ""
            centroid_struct = ""
    
    # Analyze structures
    mfe_analysis = analyze_structure(mfe_struct, seq, mfe)
    centroid_analysis = analyze_structure(centroid_struct, seq, mfe) if centroid_struct else mfe_analysis.copy()
    
    # Determine confidence level
    if struct_distance == 0:
        confidence = "IDENTICAL"
    elif struct_distance < 0.1:
        confidence = "HIGH"
    elif struct_distance < 0.2:
        confidence = "MEDIUM"
    else:
        confidence = "LOW"
    
    return {
        'id': record.id,
        'sequence': seq,
        'mfe_structure': mfe_struct,
        'centroid_structure': centroid_struct,
        'structure_distance': struct_distance,
        'base_pair_distance': bp_distance,
        'confidence': confidence,
        'mfe': mfe if mfe else 0,
        'length': seq_length,
        'loop_count': mfe_analysis['loop_count'],
        'terminal_loop_length': mfe_analysis['terminal_loop_length'],
        'terminal_loop_gc': mfe_analysis['terminal_loop_gc'],
        'stem_segments': mfe_analysis['stem_segments'],
        'pairing_percent': mfe_analysis['pairing_percent'],
        'centroid_loop_count': centroid_analysis['loop_count'],
        'centroid_terminal_loop_length': centroid_analysis['terminal_loop_length'],
        'centroid_pairing_percent': centroid_analysis['pairing_percent'],
        'centroid_stem_segments': centroid_analysis['stem_segments'],
        'gc_content': mfe_analysis['gc_content'],
        'mfei': mfe_analysis['mfei']
    }

def analyze_structure(structure, seq, mfe=None):
    """Extract features from a single structure"""
    if not structure or not isinstance(structure, str) or structure == "":
        return {
            'structure': "",
            'loop_count': 0,
            'terminal_loop_length': 0,
            'terminal_loop_gc': 0,
            'stem_segments': 0,
            'pairing_percent': 0,
            'gc_content': calculate_gc_content(seq),
            'mfei': 0
        }
    
    seq_length = len(seq)
    
    # Loop analysis
    loop_count = count_loops_simple(structure)
    terminal_loop_length = get_terminal_loop_length(structure)
    terminal_loop_sequence = get_terminal_loop_sequence(seq, structure)
    terminal_loop_gc = calculate_gc_content(terminal_loop_sequence) if terminal_loop_sequence else 0
    
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
        'structure': structure,
        'loop_count': loop_count,
        'terminal_loop_length': terminal_loop_length,
        'terminal_loop_gc': terminal_loop_gc,
        'stem_segments': stem_segments,
        'pairing_percent': pairing_percent,
        'gc_content': gc_content,
        'mfei': mfei
    }

def calculate_structure_distance(struct1, struct2):
    """Calculate the structural distance between two dot-bracket strings"""
    if not struct1 or not struct2:
        return 0.0
    if not isinstance(struct1, str) or not isinstance(struct2, str):
        return 0.0
    if len(struct1) != len(struct2):
        return 1.0
    
    differences = sum(1 for i in range(len(struct1)) if struct1[i] != struct2[i])
    return differences / len(struct1)

def calculate_base_pair_distance(struct1, struct2):
    """Calculate base pair distance between two structures"""
    def get_base_pairs(structure):
        pairs = set()
        stack = []
        if not isinstance(structure, str):
            return pairs
        for i, char in enumerate(structure):
            if char == '(':
                stack.append(i)
            elif char == ')' and stack:
                j = stack.pop()
                pairs.add((j, i))
        return pairs
    
    if not struct1 or not struct2:
        return 0.0
    if not isinstance(struct1, str) or not isinstance(struct2, str):
        return 0.0
    
    pairs1 = get_base_pairs(struct1)
    pairs2 = get_base_pairs(struct2)
    
    if len(pairs1) == 0 and len(pairs2) == 0:
        return 0.0
    
    diff = len(pairs1.symmetric_difference(pairs2))
    total = max(len(pairs1), len(pairs2))
    return diff / total if total > 0 else 0

def count_loops_simple(structure):
    """Count loops in dot-bracket structure"""
    if not structure or not isinstance(structure, str):
        return 0
    
    loop_count = 0
    
    # Count transitions from stem to loop to stem
    for i in range(len(structure) - 1):
        if structure[i] == ')' and structure[i+1] == '(':
            loop_count += 1
    
    # Count hairpin loops
    try:
        hairpins = re.findall(r'\(+\.+\)+', structure)
        loop_count += len(hairpins)
    except:
        pass
    
    # If no loops found but there are stems, it's a simple hairpin (1 loop)
    if loop_count == 0 and ('(' in structure or ')' in structure):
        loop_count = 1
    
    return loop_count

def get_terminal_loop_length(structure):
    """Extract terminal loop length"""
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
    """Extract the sequence of the terminal loop"""
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
    """Calculate GC percentage"""
    if not seq:
        return 0
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100

def calculate_mfei(mfe, length, gc_content):
    """Calculate MFE Index"""
    if gc_content == 0 or mfe is None:
        return 0
    mfe_per_nt = abs(mfe) / length
    return mfe_per_nt / (gc_content / 100)

def is_mirna_candidate(analysis, thresholds, mode='strict'):
    """
    Apply filtering thresholds based on mammalian MirGeneDB analysis
    
    Args:
        analysis: Dictionary with analysis results
        thresholds: Threshold dictionary (STRICT_THRESHOLDS or PERMISSIVE_THRESHOLDS)
        mode: 'strict' or 'permissive' (for reporting only)
    
    Returns:
        (is_valid, reasons, warnings)
    """
    reasons = []
    warnings = []
    
    # ===== HARD BINARY RULES (Immediate reject if failed) =====
    
    # 1. Loop count
    if analysis['loop_count'] > thresholds['loop_count']:
        reasons.append(f"loops={analysis['loop_count']} (max {thresholds['loop_count']})")
        return False, reasons, warnings
    
    # 2. Centroid loop count must match MFE loop count
    if thresholds.get('require_centroid_match', True):
        if analysis['centroid_loop_count'] != analysis['loop_count']:
            reasons.append(f"loop_mismatch({analysis['loop_count']}→{analysis['centroid_loop_count']})")
            return False, reasons, warnings
    
    # 3. Centroid cannot have too many loops
    if analysis['centroid_loop_count'] > thresholds['max_centroid_loops']:
        reasons.append(f"cent_loops={analysis['centroid_loop_count']} (>{thresholds['max_centroid_loops']})")
        return False, reasons, warnings
    
    # ===== RANGE-BASED CHECKS (Record warnings but don't reject unless in strict mode) =====
    
    # Terminal loop length
    tl = analysis['terminal_loop_length']
    if tl < thresholds['terminal_loop_min']:
        warnings.append(f"term_loop={tl} (<{thresholds['terminal_loop_min']})")
    elif tl > thresholds['terminal_loop_max']:
        warnings.append(f"term_loop={tl} (>{thresholds['terminal_loop_max']})")
    
    # Pairing percentage
    pp = analysis['pairing_percent']
    if pp < thresholds['pairing_percent']:
        warnings.append(f"pairing={pp:.0f}% (<{thresholds['pairing_percent']}%)")
    
    # MFEI
    mfei_val = analysis['mfei']
    if mfei_val < thresholds['mfei']:
        warnings.append(f"MFEI={mfei_val:.2f} (<{thresholds['mfei']})")
    
    # Length
    length = analysis['length']
    if length < thresholds['length_min']:
        warnings.append(f"length={length} (<{thresholds['length_min']})")
    elif length > thresholds['length_max']:
        warnings.append(f"length={length} (>{thresholds['length_max']})")
    
    # Structure distance
    sd = analysis['structure_distance']
    if sd > thresholds['structure_distance']:
        warnings.append(f"dist={sd:.3f} (>{thresholds['structure_distance']})")
    
    # Loop ratio (centroid vs MFE terminal loop)
    if analysis['centroid_terminal_loop_length'] > 0 and analysis['terminal_loop_length'] > 0:
        ratio = analysis['centroid_terminal_loop_length'] / analysis['terminal_loop_length']
        if ratio > thresholds['max_loop_ratio']:
            warnings.append(f"loop_ratio={ratio:.2f} (>{thresholds['max_loop_ratio']})")
        elif ratio < (1/thresholds['max_loop_ratio']):
            warnings.append(f"loop_ratio={ratio:.2f} (<{1/thresholds['max_loop_ratio']})")
    
    # Decision: In strict mode, warnings cause rejection. In permissive, warnings are allowed.
    if mode == 'strict':
        return len(warnings) == 0, reasons + warnings, warnings
    else:
        return True, reasons, warnings

def is_fasta_file(filepath):
    """Check if a file appears to be in FASTA format"""
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    return line.startswith('>')
    except:
        pass
    return False

def process_single_file(input_file, output_file, thresholds, mode='strict', verbose=True):
    """
    Process a single FASTA file and keep only sequences passing all thresholds
    """
    results = []
    kept = 0
    filtered = 0
    errors = 0
    warnings_count = 0
    
    loop_distribution = defaultdict(int)
    confidence_distribution = defaultdict(int)
    filter_reasons = defaultdict(int)
    
    species_name = Path(input_file).name
    
    if verbose:
        print(f"\n  Analyzing {species_name}...")
    
    try:
        for record in Bio.SeqIO.parse(input_file, "fasta"):
            try:
                analysis = analyze_sequence(record)
                results.append(analysis)
                
                # Update distributions
                loop_distribution[analysis['loop_count']] += 1
                confidence_distribution[analysis['confidence']] += 1
                
                # Apply filters
                is_valid, reasons, warnings = is_mirna_candidate(analysis, thresholds, mode)
                
                if is_valid:
                    kept += 1
                    if warnings:
                        warnings_count += 1
                    if verbose and kept <= 10:  # Show first 10 kept
                        warn_str = f" [warnings: {', '.join(warnings[:2])}]" if warnings else ""
                        print(f"    ✓ KEPT: {analysis['id'][:50]} | loops={analysis['loop_count']} | "
                              f"term_loop={analysis['terminal_loop_length']} | "
                              f"pair={analysis['pairing_percent']:.0f}% | "
                              f"MFEI={analysis['mfei']:.2f}{warn_str}")
                else:
                    filtered += 1
                    reason_str = ", ".join(reasons[:3])
                    if len(reasons) > 3:
                        reason_str += f" (+{len(reasons)-3})"
                    filter_reasons[reason_str] += 1
                    if verbose and filtered <= 10:  # Show first 10 filtered
                        print(f"    ✗ FILTERED: {analysis['id'][:50]} | {reason_str}")
            except Exception as e:
                errors += 1
                if verbose:
                    print(f"    ERROR: {record.id[:50]} | {str(e)[:50]}")
    except Exception as e:
        print(f"    ERROR: Cannot parse {species_name} as FASTA: {e}")
        return None
    
    # Write output file (only kept sequences)
    with open(output_file, 'w') as out:
        for analysis in results:
            is_valid, _, _ = is_mirna_candidate(analysis, thresholds, mode)
            if is_valid:
                out.write(f">{analysis['id']}\n")
                out.write(f"{analysis['sequence']}\n")
    
    total = len(results)
    return {
        'species': species_name,
        'input_file': input_file,
        'output_file': output_file,
        'total': total,
        'kept': kept,
        'filtered': filtered,
        'errors': errors,
        'warnings': warnings_count,
        'loop_distribution': dict(loop_distribution),
        'confidence_distribution': dict(confidence_distribution),
        'filter_reasons': dict(filter_reasons)
    }

def process_folder(input_folder, output_folder, thresholds, mode='strict', verbose=True, recursive=False):
    """
    Process ALL files in a folder that appear to be in FASTA format
    """
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    
    # Get all files
    if recursive:
        all_files = []
        for root, dirs, files in os.walk(input_folder):
            for file in files:
                all_files.append(os.path.join(root, file))
    else:
        all_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) 
                    if os.path.isfile(os.path.join(input_folder, f))]
    
    # Filter for FASTA files
    fasta_files = []
    skipped_files = []
    
    for filepath in all_files:
        if is_fasta_file(filepath):
            fasta_files.append(filepath)
        else:
            skipped_files.append(filepath)
    
    if not fasta_files:
        print(f"Error: No FASTA files found in {input_folder}")
        return None
    
    # Print header
    mode_name = "STRICT" if mode == 'strict' else "PERMISSIVE"
    print(f"\n{'='*80}")
    print(f"miRNA FILTERING - MAMMALIAN OPTIMIZED ({mode_name} MODE)")
    print(f"{'='*80}")
    print(f"Input folder:  {input_folder}")
    print(f"Output folder: {output_folder}")
    print(f"Filter mode:   {mode_name}")
    print(f"Thresholds:    {mode_name} thresholds from 7,556 mammalian miRNAs")
    print(f"Files found:   {len(fasta_files)} FASTA files")
    if skipped_files:
        print(f"Skipped:       {len(skipped_files)} non-FASTA files")
    print(f"{'='*80}\n")
    
    all_stats = []
    total_sequences = 0
    total_kept = 0
    total_warnings = 0
    
    # Process each file
    for fasta_file in sorted(fasta_files):
        input_path = Path(fasta_file)
        output_path = Path(output_folder) / input_path.name
        
        stats = process_single_file(fasta_file, output_path, thresholds, mode, verbose)
        if stats:
            all_stats.append(stats)
            total_sequences += stats['total']
            total_kept += stats['kept']
            total_warnings += stats['warnings']
            
            if verbose:
                kept_pct = (stats['kept'] / stats['total'] * 100) if stats['total'] > 0 else 0
                print(f"    → Saved {stats['kept']}/{stats['total']} ({kept_pct:.1f}%) to {output_path.name}")
    
    # Print summary
    print(f"\n{'='*80}")
    print(f"SUMMARY - {mode_name} MODE")
    print(f"{'='*80}")
    print(f"{'Species':<30} {'Total':<8} {'Kept':<8} {'Kept %':<8} {'Warnings':<8}")
    print(f"{'-'*70}")
    
    for stats in all_stats:
        kept_percent = (stats['kept'] / stats['total'] * 100) if stats['total'] > 0 else 0
        print(f"{stats['species'][:30]:<30} {stats['total']:<8} {stats['kept']:<8} {kept_percent:<7.1f}% {stats['warnings']:<8}")
    
    print(f"{'-'*70}")
    kept_percent_total = (total_kept / total_sequences * 100) if total_sequences > 0 else 0
    print(f"{'TOTAL':<30} {total_sequences:<8} {total_kept:<8} {kept_percent_total:<7.1f}% {total_warnings:<8}")
    
    # Additional statistics
    if total_sequences > 0:
        print(f"\n{'='*80}")
        print("AGGREGATE STATISTICS")
        print(f"{'='*80}")
        
        # Loop count distribution
        all_loops = defaultdict(int)
        for stats in all_stats:
            for loops, count in stats['loop_distribution'].items():
                all_loops[loops] += count
        
        print(f"\nLoop count distribution (all species combined):")
        for loops in sorted(all_loops.keys()):
            count = all_loops[loops]
            percent = (count / total_sequences) * 100
            bar = '█' * int((count / max(all_loops.values())) * 40) if all_loops else ''
            print(f"  {loops} loop(s): {count:6} ({percent:5.1f}%) {bar}")
        
        # Top filter reasons
        all_filters = defaultdict(int)
        for stats in all_stats:
            for reason, count in stats['filter_reasons'].items():
                all_filters[reason] += count
        
        if all_filters:
            print(f"\nTop filter reasons:")
            for reason, count in sorted(all_filters.items(), key=lambda x: -x[1])[:10]:
                print(f"  {reason}: {count}")
    
    print(f"\n{'='*80}")
    print(f"✓ Processing complete!")
    print(f"  Mode: {mode_name}")
    print(f"  Input files: {len(fasta_files)} FASTA files processed")
    print(f"  Total sequences: {total_sequences}")
    print(f"  Kept: {total_kept} ({kept_percent_total:.1f}%)")
    print(f"  Output folder: {output_folder}")
    print(f"{'='*80}\n")
    
    return all_stats

def main():
    # Parse command line arguments
    if len(sys.argv) < 3:
        print("Usage: python filter_mirna_mammalian.py <input_folder> <output_folder> [options]")
        print("\nDESCRIPTION:")
        print("  Filter miRNA precursors using thresholds derived from 7,556 mammalian")
        print("  MirGeneDB sequences. Processes all FASTA files in a folder.")
        print("\nOPTIONS:")
        print("  --permissive      Use permissive mode (captures let-7 and edge cases)")
        print("  --recursive       Search subdirectories recursively")
        print("  --quiet           Reduce output verbosity")
        print("\nTHRESHOLDS:")
        print("  STRICT mode (default):")
        print("    - Loop count: = 1")
        print("    - Terminal loop: 3-15 nt")
        print("    - Pairing: ≥55%")
        print("    - MFEI: ≥0.60")
        print("    - Length: 54-86 nt")
        print("    - Structure distance: ≤0.15")
        print("\n  PERMISSIVE mode (--permissive):")
        print("    - Loop count: ≤2")
        print("    - Terminal loop: ≥2 nt (allows let-7)")
        print("    - Pairing: ≥50%")
        print("    - MFEI: ≥0.50")
        print("    - Length: 50-86 nt")
        print("    - Structure distance: ≤0.25")
        print("\nEXAMPLES:")
        print("  # Strict mode (high specificity)")
        print("  python filter_mirna_mammalian.py precursors/ filtered/")
        print()
        print("  # Permissive mode (high sensitivity, includes let-7)")
        print("  python filter_mirna_mammalian.py precursors/ filtered/ --permissive")
        print()
        print("  # Process recursively")
        print("  python filter_mirna_mammalian.py precursors/ filtered/ --permissive --recursive")
        sys.exit(1)
    
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    
    # Parse options
    permissive = '--permissive' in sys.argv
    recursive = '--recursive' in sys.argv
    verbose = '--quiet' not in sys.argv
    
    # Select thresholds based on mode
    if permissive:
        thresholds = PERMISSIVE_THRESHOLDS
        mode = 'permissive'
    else:
        thresholds = STRICT_THRESHOLDS
        mode = 'strict'
    
    # Validate input folder
    if not os.path.isdir(input_folder):
        print(f"Error: Input folder '{input_folder}' does not exist.")
        sys.exit(1)
    
    # Process all files
    process_folder(input_folder, output_folder, thresholds, mode, verbose, recursive)

if __name__ == "__main__":
    main()