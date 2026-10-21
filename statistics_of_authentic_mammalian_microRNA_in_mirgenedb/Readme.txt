Description

This Python script (analyse_mirgenedb.py) performs comprehensive structural analysis of mammalian microRNA (miRNA) precursor sequences from MirGeneDB. It folds sequences using both Minimum Free Energy (MFE) and centroid structure prediction methods, extracting key features and generating statistical summaries and visualizations.

Features

Dual Structure Prediction: Analyzes sequences using both MFE and centroid folding algorithms for robust structural characterization

Comprehensive Feature Extraction:
1) Apical loop count and apical loop length
2) Pairing percentage
3) Minimum Free Energy Index (MFEI)
4) Precursor sequence length
5) Structure distance between MFE and centroid predictions
6) Base pair distance

Statistical Analysis: Generates detailed statistics including mean, median, standard deviation, and percentiles

Visualization Outputs:
1) Distribution plots for all features (individual and combined)
2) Correlation matrix showing feature relationships
3) Loop count agreement analysis between MFE and centroid structures
4) Detailed Loop Analysis: Compares loop counts between MFE and centroid structures with agreement/disagreement metrics

Dependencies:
1) ViennaRNA 
2) Python 3.6+

Usage:
python analyse_mirgenedb.py <input fasta> [output_dir]   

#input: mammalian microRNA precursor sequences curated in MirGeneDB (v3) in FASTA format