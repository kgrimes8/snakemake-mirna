import gzip


def convert_to_tsv(list: list) -> str:
    return '\t'.join(list)


def remove_brackets(string: str) -> str:
    return string.strip('()')


def parse_file(path: str) -> list:
    data_list = []
    lines = None

    # read the gz file as a string and split into lines
    with gzip.open(path, "rt") as file:
        lines = file.read().split("\n")

    # read through the file one line at a time
    for count, line in enumerate(lines):
        if line.startswith(">"):

            if count != 0:
                # new premirna so append existing data
                line_data[3] = remove_brackets(line_data[3])
                data_list.append(convert_to_tsv(line_data))

            # create new line_data, remove >
            line_data = [line.replace('>','')]

        else:
            for item in line.split(' '):
                line_data.append(item)
                
    # add last premirna to datalist
    line_data[3] = remove_brackets(line_data[3])
    data_list.append(convert_to_tsv(line_data))

    return data_list


def write_file(data: list, file: str) -> None:
    with gzip.open(file, "wt") as ofh:
        ofh.write("\n".join(data))


if __name__ == "__main__":

    structure_file = snakemake.input[0]
    output_file = snakemake.output[0]

    output_data = parse_file(structure_file)

    write_file(output_data, output_file)






