# ACCESSIBILITY
#LNCRNA
./accessibility_processing.sh ../data/latest1000/functional-lncrna-exon1-dataset-features.csv functional-lncrna-exon1-accessibility-feature
./accessibility_processing.sh ../data/latest1000/functional-lncrna-exon2-dataset-features.csv functional-lncrna-exon2-accessibility-feature

./accessibility_processing.sh ../data/latest1000/lncrna-exon1-negative-control-dataset-features.csv lncrna-exon1-negative-control-accessibility-feature
./accessibility_processing.sh ../data/latest1000/lncrna-exon2-negative-control-dataset-features.csv lncrna-exon2-negative-control-accessibility-feature

#Short NCRNA
./accessibility_processing.sh ../data/latest1000/functional-short-ncrna-dataset-features.csv functional-short-ncrna-accessibility-feature

./accessibility_processing.sh ../data/latest1000/short-ncrna-negative-control-dataset-features.csv short-ncrna-negative-control-accessibility-feature

#Protein Coding
./accessibility_processing.sh ../data/latest1000/functional-protein-exon2-dataset-features.csv functional-protein-exon2-accessibility-feature
./accessibility_processing.sh ../data/latest1000/functional-protein-exon3-dataset-features.csv functional-protein-exon3-accessibility-feature

./accessibility_processing.sh ../data/latest1000/protein-exon2-negative-control-dataset-features.csv protein-exon2-negative-control-accessibility-feature
./accessibility_processing.sh ../data/latest1000/protein-exon3-negative-control-dataset-features.csv protein-exon3-negative-control-accessibility-feature
