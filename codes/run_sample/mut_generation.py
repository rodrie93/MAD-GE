import sys

def revcom(string):
    mapping = {'A': 'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    com = ""
    for i in range(len(string)):
        com += mapping.get(string.upper()[i], 'N')
    return com[::-1]

def filter_indels(bed_in, bed_out):
    print(f"Filtering InDels from {bed_in} (<30nt) into {bed_out}...")
    with open(bed_in) as fh, open(bed_out, "w") as ofh:
        for line in fh:
            row = line.split()
            # Fixed Indices: 5 is REF, 6 is ALT
            if len(row) > 6 and len(row[5]) < 30 and len(row[6]) < 30:
                ofh.write(line)

def revcom_negative_exons(fasta_in, fasta_out):
    print(f"Applying revcom to negative strand exons from {fasta_in}...")
    with open(fasta_in, 'r') as f_in, open(fasta_out, 'w') as f_out:
        header = None
        sequence = ''
        for line in f_in:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    if header.endswith('-'): 
                        f_out.write(header + '\n' + revcom(sequence) + '\n')
                    else:
                        f_out.write(header + '\n' + sequence + '\n')
                header = line
                sequence = ''
            else:
                sequence += line
                
        if header is not None:
            if header.endswith('-'):
                f_out.write(header + '\n' + revcom(sequence) + '\n')
            else:
                f_out.write(header + '\n' + sequence + '\n')

def read_fasta(fasta_file):
    seqs = {}
    with open(fasta_file) as fh:
        header = None
        seq = []
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(seq)
                header = line[1:] 
                seq = []
            else:
                seq.append(line)
        if header:
            seqs[header] = "".join(seq)
    return seqs

def mutate_exons(mut_gff, filtered_bed, exons_fasta, out_fasta):
    print("Parsing reference exome and somatic calls to create mutant exons...")
    
    wt_exons = read_fasta(exons_fasta)
    
    mutations = []
    with open(filtered_bed) as fh:
        for line in fh:
            row = line.split()
            chrom, start, end = row[0], int(row[1]), int(row[2])
            
            # FIXED INDICES: Extract actual sequences instead of "PASS"
            ref_allele, alt_allele = row[5], row[6]
            
            mutations.append({"chrom": chrom, "start": start, "ref": ref_allele, "alt": alt_allele})
            
    mutated_seqs = {}
    with open(mut_gff) as fh:
        for line in fh:
            row = line.strip().split("\t")
            chrom = row[0]
            exon_id = row[2]
            exon_start = int(row[3]) - 1 # 0-based for list indexing
            exon_strand = row[6]
            
            fasta_key = f"{exon_id}{exon_strand}"
            if fasta_key not in wt_exons:
                continue
                
            seq = list(wt_exons[fasta_key])
            
            # Find mutations inside this exon and sort in reverse order 
            # to prevent sequence shifting issues when replacing InDels
            exon_muts = [m for m in mutations if m["chrom"] == chrom and m["start"] >= exon_start and m["start"] < (exon_start + len(seq))]
            exon_muts.sort(key=lambda x: x["start"], reverse=True)
            
            for mut in exon_muts:
                local_pos = mut["start"] - exon_start
                ref_len = len(mut["ref"])
                seq[local_pos:local_pos+ref_len] = list(mut["alt"])
                
            mutated_seqs[f"{exon_id}_mut{exon_strand}"] = "".join(seq)

    with open(out_fasta, "w") as ofh:
        for mut_id, seq in mutated_seqs.items():
            ofh.write(f">{mut_id}\n{seq}\n")
    print(f"Successfully wrote mutated exons to {out_fasta}")

def build_isoforms(tx_exon_map, wild_fasta, mut_fasta, out_fasta):
    print("Building patient-specific mutant isoforms from map...")
    
    wt_exons = read_fasta(wild_fasta)
    mut_exons = read_fasta(mut_fasta)
    
    mut_base_ids = {}
    for header, seq in mut_exons.items():
        base_id = header[:-1] # Remove strand (+ or -)
        mut_base_ids[base_id] = seq

    with open(tx_exon_map) as f_map, open(out_fasta, "w") as f_out:
        for line in f_map:
            row = line.strip().split("\t")
            enst = row[0]
            
            exons_string = row[1]
            exons_list = exons_string.split(",") if "," in exons_string else exons_string.split(";")
            
            ensg = row[2] if len(row) > 2 else ""
            
            isoform_seq = ""
            is_mutant = False
            
            for ex in exons_list:
                exon_id = ex.rsplit(":", 1)[0]
                mut_id = f"{exon_id}_mut"
                
                if mut_id in mut_base_ids:
                    isoform_seq += mut_base_ids[mut_id]
                    is_mutant = True
                else:
                    if f"{exon_id}+" in wt_exons:
                        isoform_seq += wt_exons[f"{exon_id}+"]
                    elif f"{exon_id}-" in wt_exons:
                        isoform_seq += wt_exons[f"{exon_id}-"]
            
            if is_mutant:
                f_out.write(f">{enst}_mut {ensg}\n{isoform_seq}\n")
                
    print(f"Successfully wrote mutant isoforms to {out_fasta}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python mut_generation.py [filter_indels|revcom|mutate_exons|build_isoforms] <args>")
        sys.exit(1)

    step = sys.argv[1]
    
    if step == "filter_indels":
        filter_indels(sys.argv[2], sys.argv[3])
    elif step == "revcom":
        revcom_negative_exons(sys.argv[2], sys.argv[3])
    elif step == "mutate_exons":
        mutate_exons(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    elif step == "build_isoforms":
        build_isoforms(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
