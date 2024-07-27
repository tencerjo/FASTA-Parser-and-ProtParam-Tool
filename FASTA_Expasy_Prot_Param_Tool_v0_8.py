# v 0.8
# Adjusted to place N-terminus at the beginning and C-terminus at the end
# Added the option to enter "Fc" for the C-terminus

import os
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import csv

Fc_sequence = ("EPKSSDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLYITREPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSC"
              "SVMHEALHNHYTQKSLSLSPGK")

def calculate_charge_at_pH7(sequence):
    analysed_seq = ProteinAnalysis(sequence)
    return analysed_seq.charge_at_pH(7.0)

def calculate_percent_charged(sequence):
    charged_aa = sum(sequence.count(aa) for aa in 'DEHKR')
    total_aa = len(sequence)
    return (charged_aa / total_aa) * 100

def calculate_percent_hydrophobic(sequence):
    hydrophobic_aa = sum(sequence.count(aa) for aa in 'AILMFWYV')
    total_aa = len(sequence)
    return (hydrophobic_aa / total_aa) * 100

def calculate_molar_extinction_coefficient_oxidized(sequence):
    analysed_seq = ProteinAnalysis(sequence)
    epsilon_prot = analysed_seq.molar_extinction_coefficient()
    return epsilon_prot[1]

def calculate_abs_0_1_percent_oxidized(sequence):
    extinction_coefficient = calculate_molar_extinction_coefficient_oxidized(sequence)
    molecular_weight = ProteinAnalysis(sequence).molecular_weight()
    return extinction_coefficient / molecular_weight

def parse_fasta_and_export_with_protparam(input_file, c_terminus="", n_terminus="", multimer_count=1):
    sequences = []
    
    with open(input_file, 'r') as file:
        sequence_id = None
        sequence_lines = []
        
        for line in file:
            line = line.strip()
            
            if line.startswith('>'):
                if sequence_id is not None and sequence_lines:
                    sequence = ''.join(sequence_lines)
                    if c_terminus == "Fc":
                        c_terminus_seq = Fc_sequence
                    else:
                        c_terminus_seq = c_terminus
                    sequence_with_modifications = n_terminus + (sequence * multimer_count) + c_terminus_seq
                    
                    charge_at_pH7 = calculate_charge_at_pH7(sequence_with_modifications)
                    percent_charged = calculate_percent_charged(sequence_with_modifications)
                    percent_hydrophobic = calculate_percent_hydrophobic(sequence_with_modifications)
                    extinction_coefficient = calculate_molar_extinction_coefficient_oxidized(sequence_with_modifications)
                    abs_0_1_percent = calculate_abs_0_1_percent_oxidized(sequence_with_modifications)
                    isoelectric_point = ProteinAnalysis(sequence_with_modifications).isoelectric_point()
                    molecular_weight = ProteinAnalysis(sequence_with_modifications).molecular_weight()
                    
                    sequences.append((sequence_id, percent_charged, percent_hydrophobic, charge_at_pH7, isoelectric_point, extinction_coefficient, molecular_weight, abs_0_1_percent, sequence, c_terminus_seq, n_terminus, multimer_count, sequence_with_modifications))
                
                sequence_id = line[1:]
                sequence_lines = []
            else:
                sequence_lines.append(line)
        
        if sequence_id is not None and sequence_lines:
            sequence = ''.join(sequence_lines)
            if c_terminus == "Fc":
                c_terminus_seq = Fc_sequence
            else:
                c_terminus_seq = c_terminus
            sequence_with_modifications = n_terminus + (sequence * multimer_count) + c_terminus_seq
            
            charge_at_pH7 = calculate_charge_at_pH7(sequence_with_modifications)
            percent_charged = calculate_percent_charged(sequence_with_modifications)
            percent_hydrophobic = calculate_percent_hydrophobic(sequence_with_modifications)
            extinction_coefficient = calculate_molar_extinction_coefficient_oxidized(sequence_with_modifications)
            abs_0_1_percent = calculate_abs_0_1_percent_oxidized(sequence_with_modifications)
            isoelectric_point = ProteinAnalysis(sequence_with_modifications).isoelectric_point()
            molecular_weight = ProteinAnalysis(sequence_with_modifications).molecular_weight()
            
            sequences.append((sequence_id, percent_charged, percent_hydrophobic, charge_at_pH7, isoelectric_point, extinction_coefficient, molecular_weight, abs_0_1_percent, sequence, c_terminus_seq, n_terminus, multimer_count, sequence_with_modifications))
    
    output_dir = os.path.dirname(input_file)
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_csv = os.path.join(output_dir, f'{base_name}_ProtParam_Results.csv')
    
    with open(output_csv, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['Sequence Name', '% Charged Amino Acids', '% Hydrophobic Amino Acids', 'Charge at pH 7', 'Isoelectric Point', 'Extinction Coefficient M-1cm-1 Oxidized', 'MW in daltons', 'Abs 0.1% (1g/L) Oxidized', 'Sequence', 'C-Terminus', 'N-Terminus', 'Multimer Count', 'Analyzed Sequence'])
        for seq_id, percent_charged, percent_hydrophobic, charge_pH7, pI, extinction_coefficient, mw, abs_0_1_percent, seq, c_term, n_term, multimer, seq_with_modifications in sequences:
            c_term_display = c_term if c_term else "No C-terminus added"
            n_term_display = n_term if n_term else "No N-terminus added"
            writer.writerow([seq_id, percent_charged, percent_hydrophobic, charge_pH7, pI, extinction_coefficient, mw, abs_0_1_percent, seq, c_term_display, n_term_display, multimer, seq_with_modifications])
    
    print(f"Sequences exported to '{output_csv}'")
    input("Press the <ENTER> key to continue...")

def find_fasta_file(root_folder):
    for root, dirs, files in os.walk(root_folder):
        for file in files:
            if file.endswith('.fasta'):
                return os.path.join(root, file)
    return None

print("The directory must contain a .fasta file to be used.")
root_folder = input("Please enter the path to the root folder: ")
input_file = find_fasta_file(root_folder)

if input_file:
    print(f"Found FASTA file: {input_file}")
    
    multimer_count_input = input("Enter the multimer count (Example: dimer = 2): ")
    multimer_count = int(multimer_count_input) if multimer_count_input else 1
    n_terminus = input("Enter the N-Terminus to add at the beginning (Press enter for no N-Terminus): ")
    c_terminus = input("Enter the C-Terminus to add at the end (Press enter for no C-Terminus, or type 'Fc' for Fc sequence): ")
    
    parse_fasta_and_export_with_protparam(input_file, c_terminus, n_terminus, multimer_count)
else:
    print("No FASTA file found in the specified root folder.")
