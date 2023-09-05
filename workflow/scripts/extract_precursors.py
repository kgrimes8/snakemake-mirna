import gzip
import sys

from Bio import SeqIO


def read_line(line: str) -> dict:
    columns = line.split("\t")
    return {
        "chr": columns[0],
        "start": columns[3],
        "end": columns[4],
        "strand": columns[6],
        "attributes": columns[8],
    }


def parse_gff(path: str) -> list:
    dict_array = []
    # read the gz file as a string and split into lines
    with gzip.open(path, "rt") as file:
        lines = file.read().split("\n")
        # read through the file one line at a time
        for line in lines:
            if "pre_miRNA" in line:
                dict_array.append(read_line(line))
    return dict_array


def parse_fasta(path: str) -> dict:
    with gzip.open(path, "rt") as fa_file:
        seq_dict = SeqIO.to_dict(SeqIO.parse(fa_file, "fasta"))
    return seq_dict


def get_seqs(premirna_array: list, fasta_dict: dict) -> list:
    data_list = []
    for premirna_dict in premirna_array:
        chr_seq = fasta_dict[premirna_dict["chr"]]
        start_loc = int(premirna_dict["start"])
        end_loc = int(premirna_dict["end"])
        premirna_seq = chr_seq.seq[start_loc - 1 : end_loc]
        if premirna_dict["strand"] == "-":
            premirna_seq = premirna_seq.reverse_complement()

        # parse attributes to get transcript ID
        attributes_split = premirna_dict["attributes"].split(";")
        transcript_id = attributes_split[(len(attributes_split) - 1)].split("=")[1]

        data_list.append(f">{transcript_id}")
        data_list.append(str(premirna_seq))

    return data_list


def write_file(data: list, file: str) -> None:
    with gzip.open(file, "wt") as ofh:
        ofh.write("\n".join(data))


if __name__ == "__main__":
    with open(snakemake.log[0], "w") as f:
        sys.stderr = sys.stdout = f

        seq = snakemake.input[0]
        annot = snakemake.input[1]
        output_file = snakemake.output[0]

        print(f"File inputs: Fasta file - {seq}, Annotation file - {annot}")

        # read in the gff for pre_miRNA, return list of dicts
        print("Parse gff ...")
        premirna_array = parse_gff(annot)
        print("Complete!")

        # parse fasta using seqIO to dict
        print("Parse fasta ...")
        fasta_dict = parse_fasta(seq)
        print("Complete!")

        # get the sequence from the FASTA file for each premirna
        print("Get seqs ...")
        output_data = get_seqs(premirna_array, fasta_dict)
        print("Complete!")

        # write output file
        print("Write to file ...")
        write_file(output_data, output_file)
        print("Complete!")
