Description

This Python script (vienna_filter_recom_mam.py) folds and filters mammalian miRNA precursor sequences in batches. It predicts both minimum free energy (MFE) and centroid RNA structures using ViennaRNA, extracts structural and thermodynamic features, and retains sequences that satisfy mammalian miRNA precursor criteria derived from 7,556 MirGeneDB sequences.

The folding configuration corresponds to RNAfold -p -d2 --noLP at 37 degrees Celsius and is designed for ViennaRNA 2.7.2.

Features

RNA Structure Prediction:
1) Predicts the MFE structure
2) Calculates the partition function and centroid structure
3) Uses dangle model 2 and disables isolated base pairs
4) Falls back to default folding settings if the configured fold fails

Feature Extraction:
1) MFE and centroid dot-bracket structures
2) MFE and centroid loop counts
3) Terminal-loop length and terminal-loop GC content
4) Pairing percentage
5) Number of stem segments
6) Sequence length and GC content
7) Minimum Free Energy Index (MFEI)
8) Normalized character-level structure distance
9) Normalized base-pair distance
10) Structure-confidence category based on MFE-versus-centroid distance

Filtering Modes

Strict Mode (default):
1) MFE loop count must equal 1
2) MFE and centroid loop counts must match
3) Centroid loop count must not exceed 2
4) Terminal-loop length must be 3-15 nucleotides
5) Pairing percentage must be at least 55%
6) MFEI must be at least 0.60
7) Sequence length must be 54-86 nucleotides
8) Structure distance must not exceed 0.15
9) Centroid-to-MFE terminal-loop ratio must remain within the configured threefold range

All strict-mode criteria must pass for a sequence to be retained.

Permissive Mode:
1) MFE loop count must not exceed 2
2) MFE and centroid loop counts must match
3) Centroid loop count must not exceed 2
4) Terminal-loop length target is 2-22 nucleotides
5) Pairing-percentage target is at least 50%
6) MFEI target is at least 0.50
7) Sequence-length target is 50-86 nucleotides
8) Structure-distance target is no greater than 0.25
9) The configured centroid-to-MFE terminal-loop ratio range is fivefold

In permissive mode, the loop-count rules remain mandatory. Failures of the range-based criteria are reported as warnings but do not cause rejection.

Input

The input is a folder containing one or more FASTA files. A file is recognized as FASTA when its first non-empty line begins with the > character; the filename extension is not used for detection.

Each FASTA record should contain:
1) A unique sequence identifier in the header
2) An RNA or DNA nucleotide sequence representing a candidate miRNA precursor

The script converts sequences to uppercase before analysis.

Output

The script creates the requested output folder and writes one filtered FASTA file for each input FASTA file, preserving the original filename. Each output file contains only the sequences retained under the selected filtering mode.

The terminal report includes:
1) Total, retained, filtered, and error counts for each input file
2) Retention percentages
3) Counts of retained sequences with warnings
4) Aggregate loop-count distribution
5) The ten most frequent filtering-reason combinations

The script does not write the detailed feature measurements or summary statistics to a separate table.

Dependencies

1) Python 3
2) ViennaRNA Python bindings, designed for version 2.7.2
3) Biopython

Example installation with Conda:

conda install -c bioconda viennarna biopython

Usage

Strict mode:

python vienna_filter_recom_mam.py <input_folder> <output_folder>

Permissive mode:

python vienna_filter_recom_mam.py <input_folder> <output_folder> --permissive

Recursive folder search:

python vienna_filter_recom_mam.py <input_folder> <output_folder> --recursive

Reduced terminal output:

python vienna_filter_recom_mam.py <input_folder> <output_folder> --quiet

Options can be combined:

python vienna_filter_recom_mam.py <input_folder> <output_folder> --permissive --recursive --quiet

Interpretation Notes

1) The thresholds are empirical screening criteria rather than proof that a retained sequence is an authentic miRNA precursor.
2) Permissive mode retains sequences that pass the mandatory loop-count rules even when other measurements fall outside their target ranges.
3) A structure distance of zero is classified as IDENTICAL, below 0.1 as HIGH confidence, below 0.2 as MEDIUM confidence, and 0.2 or higher as LOW confidence.
4) If centroid prediction fails, the script substitutes the MFE structure, which produces an apparent structure distance of zero.
5) When recursive processing is enabled, outputs from different subdirectories are written into the same output folder. Input files with identical basenames may overwrite one another.
6) Existing files with the same names in the output folder are overwritten.
