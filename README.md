Built on Anaconda 2.6.0, Python 3.11.7

This FASTA Parser tools script looks for text starting with ">" and grabs the amino acid sequence below it until two line breaks are detected (rather than just looking for the next ">".

The script then uses BioPython to calculate various biochemical properties and writes them to a csv file along with the sequence name, original sequence, and the sequence with an appended C-terminus sequence (used for biochemical calculation)
% Charged Amino Acids
% Hydrophobic Amino Acids
Charge at pH 7
Molecular Extinction Coefficient Oxidized
Abs 0.1% Oxidized (Molecular Extinction Coefficient Oxidized / Molecular Weight)
Isoelectric Point
Molecular Weight in daltons
