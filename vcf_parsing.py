import os


def parse_vcf(path: str) -> None
"""
Parse VCF files with the suffix '_rare.vcf' after snakemake pipeline in the specified directory.

Args:
path (str): The path to the directory containing VCF files.

Returns:
None: The function does not return a value but writes processed data to new files.

Processed Data:
- Extracts specific information from VCF files and writes it to a new text file.
- Information extracted includes position, reference allele, alternate allele, and variant frequency.
"""
path = path

for file in os.listdir(path):
    if file.endswith('_rare.vcf'):
        file_path = os.path.join(path, file)
        with open (file_path, mode='r') as rare_vcf:
            rare_var = rare_vcf.readlines()
            destination = []
            for line in rare_var:
                dest_str = []
                if not line.startswith('#'):
                    line = line.split('\t')
                    dest_str.append(line[1])
                    dest_str.append(line[3])
                    dest_str.append(line[4])
                    subline1 = line[8].split(':')
                    subline2 = line[9].split(':')
                    for i in range(len(subline1)):
                        if subline1[i] == 'FREQ':
                            dest_str.append(subline2[i])
                    destination.append('\t'.join(dest_str))
            
            output_file_name = file.replace('_rare.vcf', '_processed_rare.txt')
            output_file_path = os.path.join(path, output_file_name)
            with open(output_file_path, mode='a') as output_file:
                output_file.write('\n'.join(destination))