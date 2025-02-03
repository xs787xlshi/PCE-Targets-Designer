# Programmable Chromosome Engineering (PCE) system: Designing pegRNAs for Various Large-Scale DNA Manipulations

Given the complexity of designing pegRNAs for the PCE and RePCE systems, we developed a web-based tool, 
PCE Targets Designer v1.0 (http://www.engineeringchromosome.net/), to streamline and simplify the design
 of target sites for diverse chromosome editing applications. These applications include integration, de
letion, inversion, replacement, translocation, and scarless editing. This tool aims to provide researche
rs with a convenient and efficient platform for utilizing our systems.
When designing target sites, the tool considers several critical factors:
    • The PAM types associated with the target site,
    • The distance between cleavage sites for the dual-pegRNA strategy,
    • The secondary structure formed between the target site sequence and the RT sequence, and
    • The potential off-target effects of the target site sequence.
The single pegRNA designer part of  PCE Targets Designer was developed by modifying CRISPR-GE (Xie et al
., 2017). The output includes a user-friendly list of target sites based on the input genomic DNA sequen
ces, ranked according to the factors mentioned above. For subsequent pegRNA primer design, users can ref
er to the PlantPegDesigner tool (http://www.plantgenomeediting.net/).
Currently, PCE Targets Designer v1.0 is specifically tailored for designing target sites in the rice gen
ome. We plan to enhance the tool in future versions by adding advanced functionalities, improving design
 accuracy, and expanding support for additional species’ genomes.
# PCEDesign_multiple.py (Design pegRNA for integration, deletion, inversion, replace ment and translocat
ion manipulations)
Usage: python PCEDesign_multiple.py ./example/test_PCE.fa
Input file: ./example/test_PCE.fa
Output file: ./result_PCE/results/*.TwinPair.txt
# PCEDesignScarless_multiple.py (Design pegRNA for removing scar)
Usage: python PCEDesignScarless_multiple.py ./example/test_scarless.fa
Input file: ./example/test_scarless.fa
Output files: ./ result_scarless/results/*.TwinPair.txt

Xie, X., Ma, X., Zhu, Q., Zeng, D., Li, G., Liu, YG., (2017). CRISPR-GE: A Convenient Software Toolkit f
or CRISPR-Based Genome Editing. Mol. Plant 10, 1246–1249.
