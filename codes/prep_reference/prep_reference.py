import sys

def parse_gff_and_map(gff_in, gff_out, map_out):
    tx_exons = {}
    tx_gene = {}
    unique_exons = {}

    print("Parsing GFF and building custom annotations and map...")
    with open(gff_in) as f:
        for line in f:
            if line.startswith("#"): 
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "exon": 
                continue

            # Parse attributes dynamically
            attrs = dict(x.split("=") for x in parts[8].split(";") if "=" in x)
            
            # Extract ID and strip "exon:" prefix if present
            exon_id = attrs.get("ID", "")
            if exon_id.startswith("exon:"):
                exon_id = exon_id[5:]
            if not exon_id: 
                continue

            enst = attrs.get("transcript_id", attrs.get("Parent", "").replace("transcript:", ""))
            ensg = attrs.get("gene_id", "")
            ense = attrs.get("exon_id", "")
            exon_num = attrs.get("exon_number", "0")

            if exon_id not in unique_exons:
                unique_exons[exon_id] = {
                    "chrom": parts[0], "start": parts[3], "end": parts[4],
                    "strand": parts[6], "ENSE": ense, "ENSG": ensg, "ENSTs": set()
                }
            unique_exons[exon_id]["ENSTs"].add(enst)

            if enst not in tx_exons:
                tx_exons[enst] = []
                tx_gene[enst] = ensg
            
            tx_exons[enst].append((exon_id, int(exon_num) if exon_num.isdigit() else 0))

    print(f"Writing feature GFF to {gff_out}...")
    with open(gff_out, "w") as f_gff:
        for ex_id, data in unique_exons.items():
            ensts_str = ",".join(data["ENSTs"])
            # Format requested: chr . ID start end . strand ENSE ENST ENSG
            row = [data["chrom"], ".", ex_id, data["start"], data["end"], ".", data["strand"], data["ENSE"], ensts_str, data["ENSG"]]
            f_gff.write("\t".join(row) + "\n")

    print(f"Writing transcript-exon map to {map_out}...")
    with open(map_out, "w") as f_map:
        for enst, ex_list in tx_exons.items():
            ex_list.sort(key=lambda x: x[1])
            map_str = ";".join([f"{ex_id}:{ex_num}" for ex_id, ex_num in ex_list])
            f_map.write(f"{enst}\t{map_str}\t{tx_gene.get(enst, '')}\n")

def format_fasta(fasta_in, gff_in, fasta_out):
    print(f"Appending strand information to {fasta_out}...")
    strands = {}
    with open(gff_in) as f:
        for line in f:
            if line.strip():
                parts = line.strip().split("\t")
                if len(parts) >= 7:
                    strands[parts[2]] = parts[6]  # ID is col 3, strand is col 7

    with open(fasta_in) as fin, open(fasta_out, "w") as fout:
        header = None
        seq = []
        for line in fin:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    clean_id = header[1:]
                    strand = strands.get(clean_id, "")
                    fout.write(f">{clean_id}{strand}\n")
                    fout.write("".join(seq) + "\n")
                header = line
                seq = []
            else:
                seq.append(line)
        # Process the final record
        if header is not None:
            clean_id = header[1:]
            strand = strands.get(clean_id, "")
            fout.write(f">{clean_id}{strand}\n")
            fout.write("".join(seq) + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python prep_reference.py [parse_gff|format_fasta] ...")
        sys.exit(1)

    mode = sys.argv[1]
    if mode == "parse_gff":
        parse_gff_and_map(sys.argv[2], sys.argv[3], sys.argv[4])
    elif mode == "format_fasta":
        format_fasta(sys.argv[2], sys.argv[3], sys.argv[4])