# Phylliroe Placement manuscript code repository

Citation: Goodheart, Jessica A., Wägele Heike. 2020. Swimming behavior in nudibranch gastropods facilitated the transition to a pelagic lifestyle in Phylliroe. Organisms, Diversity, and Evolution, 20: 657–667. doi: 10.1007/s13127-020-00458-9.

## Scripts Included

### calculate_basic_denovo_transcriptome_assembly_statistics.pl
* Generates output statistics (# of transcripts (TFs), # of bases, N50, L50) for a transcriptome in FASTA format
* Usage: ``` perl calculate_basic_denovo_transcriptome_assembly_statistics.pl Trinity.fasta ```

### create_data_matrix_gastropoda.pl
* Aligns protein sequences in HaMStR assigned orthologous groups using MAFFT and converts the protein alignments to corresponding nucleotide alignments.
* Usage: ``` perl create_data_matrix_gastropoda.pl con/rep [OUTPUT DIRECTORY] ``` 
* NOTE: con = consensus sequences, rep = representative sequences
