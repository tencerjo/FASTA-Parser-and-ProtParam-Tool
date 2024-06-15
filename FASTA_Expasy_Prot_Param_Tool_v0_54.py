# v0_54
# csv output does not need to be specified, script will save to same location as file input.
# File path 'not' entered by user

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv
import os

# To Add Fc at C-terminus or another sequence define append_str with the sequence in quotes.
append_str = "EPKSSDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLYITREPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSC SVMHEALHNHYTQKSLSLSPGK"

def calculate_charge_at_pH7(sequence):
    # Calculate charge at pH 7
    analysed_seq = ProteinAnalysis(sequence)
    charge_at_pH7 = analysed_seq.charge_at_pH(7.0)
    return charge_at_pH7

def calculate_percent_charged(sequence):
    # Calculate percentage of charged amino acids (DEHKR)
    charged_aa = sum(sequence.count(aa) for aa in 'DEHKR')
    total_aa = len(sequence)
    percent_charged = (charged_aa / total_aa) * 100
    return percent_charged

def calculate_percent_hydrophobic(sequence):
    # Calculate percentage of hydrophobic amino acids (AILMFWYV)
    hydrophobic_aa = sum(sequence.count(aa) for aa in 'AILMFWYV')
    total_aa = len(sequence)
    percent_hydrophobic = (hydrophobic_aa / total_aa) * 100
    return percent_hydrophobic

def calculate_molar_extinction_coefficient_oxidized(sequence):
    # Calculate molar extinction coefficient for oxidized cysteines
    analysed_seq = ProteinAnalysis(sequence)
    epsilon_prot = analysed_seq.molar_extinction_coefficient()  # [reduced, oxidized]
    return epsilon_prot[1]  # Return the oxidized form coefficient

def calculate_abs_0_1_percent_oxidized(sequence):
    # Calculate Abs 0.1% (oxidized)
    extinction_coefficient = calculate_molar_extinction_coefficient_oxidized(sequence)
    molecular_weight = ProteinAnalysis(sequence).molecular_weight()
    abs_0_1_percent = extinction_coefficient / molecular_weight
    return abs_0_1_percent

def parse_fasta_and_export_with_protparam(input_file, append_str=""):
    sequences = []
    
    with open(input_file, 'r') as file:
        sequence_id = None
        sequence_lines = []
        line_count = 0  # Initialize line counter
        
        for line in file:
            line = line.strip()
            
            if line.startswith('>'):
                # Save the current sequence if there is one
                if sequence_id is not None and sequence_lines:
                    sequence = ''.join(sequence_lines)
                    sequence_with_append = sequence + append_str
                    
                    # Calculate Protein Analysis Parameters using the appended sequence
                    charge_at_pH7 = calculate_charge_at_pH7(sequence_with_append)
                    percent_charged = calculate_percent_charged(sequence_with_append)
                    percent_hydrophobic = calculate_percent_hydrophobic(sequence_with_append)
                    extinction_coefficient = calculate_molar_extinction_coefficient_oxidized(sequence_with_append)
                    abs_0_1_percent = calculate_abs_0_1_percent_oxidized(sequence_with_append)
                    isoelectric_point = ProteinAnalysis(sequence_with_append).isoelectric_point()
                    molecular_weight = ProteinAnalysis(sequence_with_append).molecular_weight()
                    
                    # Append data to sequences list
                    sequences.append((sequence_id, percent_charged, percent_hydrophobic, charge_at_pH7, isoelectric_point, extinction_coefficient, molecular_weight, abs_0_1_percent, sequence, sequence_with_append))
                
                # Start a new sequence
                sequence_id = line[1:]
                sequence_lines = []
                line_count = 0  # Reset line counter
            else:
                # Collect up to two lines of the sequence
                if line_count < 2:
                    sequence_lines.append(line)
                    line_count += 1
        
        # Save the last sequence if it is valid
        if sequence_id is not None and sequence_lines:
            sequence = ''.join(sequence_lines)
            sequence_with_append = sequence + append_str
            
            # Calculate Protein Analysis Parameters using the appended sequence
            charge_at_pH7 = calculate_charge_at_pH7(sequence_with_append)
            percent_charged = calculate_percent_charged(sequence_with_append)
            percent_hydrophobic = calculate_percent_hydrophobic(sequence_with_append)
            extinction_coefficient = calculate_molar_extinction_coefficient_oxidized(sequence_with_append)
            abs_0_1_percent = calculate_abs_0_1_percent_oxidized(sequence_with_append)
            isoelectric_point = ProteinAnalysis(sequence_with_append).isoelectric_point()
            molecular_weight = ProteinAnalysis(sequence_with_append).molecular_weight()
            
            # Append data to sequences list
            sequences.append((sequence_id, percent_charged, percent_hydrophobic, charge_at_pH7, isoelectric_point, extinction_coefficient, molecular_weight, abs_0_1_percent, sequence, sequence_with_append))
    
    # Get the directory and base name of the input file
    output_dir = os.path.dirname(input_file)
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    # Set the output CSV file path
    output_csv = os.path.join(output_dir, f'{base_name}_ProtParam_Results.csv')
    
    # Write to CSV
    with open(output_csv, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['Sequence Name', '% Charged Amino Acids', '% Hydrophobic Amino Acids', 'Charge at pH 7', 'Isoelectric Point', 'Extinction Coefficient M-1cm-1 Oxidized', 'MW in daltons', 'Abs 0.1% (1g/L) Oxidized', 'Sequence', 'Sequence Appended'])
        for seq_id, percent_charged, percent_hydrophobic, charge_pH7, pI, extinction_coefficient, mw, abs_0_1_percent, seq, seq_with_append in sequences:
            writer.writerow([seq_id, percent_charged, percent_hydrophobic, charge_pH7, pI, extinction_coefficient, mw, abs_0_1_percent, seq, seq_with_append])

    print(f"Sequences exported to '{output_csv}'")

def find_fasta_file(root_folder):
    for root, dirs, files in os.walk(root_folder):
        for file in files:
            if file.endswith('.fasta'):
                return os.path.join(root, file)
    return None

# Usage
root_folder = input("Please enter the path to the root folder: ")
input_file = find_fasta_file(root_folder)

if input_file:
    print(f"Found FASTA file: {input_file}")
    parse_fasta_and_export_with_protparam(input_file, append_str)
else:
    print("No FASTA file found in the specified root folder.")
