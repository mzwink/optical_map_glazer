from Bio import SeqIO

path_to_files = "/Users/madisonzwink/Desktop/Glazer_project/"

scaff_location_file = open(path_to_files + "NewScaffoldOrder.tsv")
fasta_file = path_to_files + "revisedAssemblyUnmasked.fa"
#scaffold_fasta = open(path_to_files + "glazer_scaffolds.fa", 'w')

def find_scaffold_pos(location_file, chr_num):
    location_file = open(location_file)
    header = location_file.readline()
    scaffold_file = location_file.readlines()
    scaffold_positions = []

    for scaffold in scaffold_file:

        scaffold_strip = scaffold.rstrip()
        scaffold_split = scaffold.split("\t")

        chromosome = scaffold_split[2]

        if chromosome == chr_num:

            scaffold_num = scaffold_split[0]
            start_pos = scaffold_split[3]
            end_pos = scaffold_split[4]

            scaffold_positions.append([scaffold_num, start_pos, end_pos])
    return scaffold_positions

def parse_fasta(fasta_file):

    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')

    for fasta in fasta_sequences:
        chromosome = "NA"
        name, sequence = fasta.id, str(fasta.seq)
        if name == "chrI": chromosome = "1"
        if name == "chrII": chromosome = "2"
        if name == "chrIII": chromosome = "3"
        if name == "chrIV": chromosome = "4"
        if name == "chrV": chromosome = "5"
        if name == "chrVI": chromosome = "6"
        if name == "chrVII": chromosome = "7"
        if name == "chrVIII": chromosome = "8"
        if name == "chrIX": chromosome = "9"
        if name == "chrX": chromosome = "10"
        if name == "chrXI": chromosome = "11"
        if name == "chrXII": chromosome = "12"
        if name == "chrXIII": chromosome = "13"
        if name == "chrXIV": chromosome = "14"
        if name == "chrXV": chromosome = "15"
        if name == "chrXVI": chromosome = "16"
        if name == "chrXVII": chromosome = "17"
        if name == "chrXVIII": chromosome = "18"
        if name == "chrXIX": chromosome = "19"
        if name == "chrXX": chromosome = "20"
        if name == "chrXXI": chromosome = "21"
        if name == "chrM": chromosome = "M"
        if name == "chrUn": chromosome = "Un"

        output = open("chr_" + str(chromosome) + "_scaffolds.fa", 'w')
        scaffold_positions = find_scaffold_pos("/Users/madisonzwink/Desktop/Glazer_project/NewScaffoldOrder.tsv", chromosome)

        for pos in scaffold_positions:
            #i = 0
            scaffold_id = pos[0]
            start_pos = int(pos[1])
            end_pos = int(pos[2])

            scaffold = sequence[start_pos-1:end_pos-1]

            output.write(">" + str(name) + "_scaffold_" + str(scaffold_id) + "\n")
            output.write(str(scaffold) + "\n\n")

parse_fasta(fasta_file)
