
#Have the gap info for each fragment- want to summarize for each chromosome

#Go through stats, append gap lists for all fragments in same chromosome
#Make sure there is no duplicate sites, sort lists
#Then find the longest continuous gap - output this information

xmap_stats = "/Users/madisonzwink/Desktop/Glazer_project/Refaligner_output/glazer_chr_scaffolds_wrap_BspQI_to_GAST_ACUL_2014_057_DEFAULT_T_REFINEFINAL1_stats.txt"
cmap = "/Users/madisonzwink/Desktop/Glazer_project/glazer_chr_scaffolds_BspQI.cmap"

def parse_xmap_stats(xmap_stats, cmap):
    output = xmap_stats.replace(".txt", "_chr.txt")
    output = open(output, 'w')
    xmap_stats = open(xmap_stats).readlines()
    cmap_key = cmap.replace(".cmap", "_key.txt")

    output.write("CHR\tRESTR_SITES\tGAP_RESTR_SITES(TOTAL)\tCONTINUOUS_GAP_START\tCONTINUOUS_GAP_END\tGAP_LEN\tGAP_POS\n")


    for chromosome in range(1,22):

        total_restr_sites = 0
        total_gaps = []
        ref_id_list = []
        largest_gap_start = 0.0
        largest_gap_end = 0.0
        largest_gap_len = 0.0

        gap_sites = {}
        total_gap_positions = []

        if chromosome == 1: chromosome = "chrI_"
        if chromosome == 2: chromosome = "chrII_"
        if chromosome == 3: chromosome = "chrIII_"
        if chromosome == 4: chromosome = "chrIV_"
        if chromosome == 5: chromosome = "chrV_"
        if chromosome == 6: chromosome = "chrVI_"
        if chromosome == 7: chromosome = "chrVII_"
        if chromosome == 8: chromosome = "chrVIII_"
        if chromosome == 9: chromosome = "chrIX_"
        if chromosome == 10: chromosome = "chrX_"
        if chromosome == 11: chromosome = "chrXI_"
        if chromosome == 12: chromosome = "chrXII_"
        if chromosome == 13: chromosome = "chrXIII_"
        if chromosome == 14: chromosome = "chrXIV_"
        if chromosome == 15: chromosome = "chrXV_"
        if chromosome == 16: chromosome = "chrXVI_"
        if chromosome == 17: chromosome = "chrXVII_"
        if chromosome == 18: chromosome = "chrXVIII_"
        if chromosome == 19: chromosome = "chrXIX_"
        if chromosome == 20: chromosome = "chrXX_"
        if chromosome == 21: chromosome = "chrXXI_"

        for stat in xmap_stats:
            if stat.startswith("REF"):
                continue
            else:
                stat_strip = stat.rstrip()
                stat_split = stat_strip.split("\t")
                chr_id = stat_split[0]
                restr_sites = int(stat_split[2])
                gap_positions = stat_split[3]
                num_gap_sites = stat_split[5]
                continuous_gap_start = float(stat_split[6])
                continuous_gap_end = float(stat_split[7])
                continuous_gap_len = float(stat_split[8])

            if chr_id.startswith(chromosome):
                total_restr_sites += restr_sites
                total_gaps.append(gap_positions)
                ref_id_list.append(chr_id)

                #print("largest_gap_len: " + str(largest_gap_len))

                if continuous_gap_len > largest_gap_start:
                    largest_gap_start = continuous_gap_start
                    largest_gap_end = continuous_gap_end
                    largest_gap_len = continuous_gap_len

        #print(chromosome)
        #print(largest_gap_len)
        #print(ref_id_list)
        #print(total_restr_sites)

        for gap_list in total_gaps:
            #print(gap_list)
            gap_list_split = gap_list.split(",")
            for pos in gap_list_split:
                pos = pos.replace("]","").replace("[", "")

                if pos in gap_sites.keys():
                    continue
                else:
                    gap_sites[pos] = [0]

        for key in gap_sites.keys():
            #print(key)
            total_gap_positions.append(key)


        total_gap_positions = sorted(total_gap_positions)

        continuous_gap_list = []
        for pos in total_gap_positions:
            #print(pos)
            if float(pos) >= largest_gap_start and float(pos) <= largest_gap_end:
                continuous_gap_list.append(pos)


        output_string = ''
        counter = 0

        for pos in continuous_gap_list:
            pos = pos.replace("\t", "")
            pos = pos.replace(" ", "")
            if counter == 0:
                output_string += "(" + str(pos) + ", "

                counter += 1
            if counter == len(continuous_gap_list):
                output_string += str(pos) + ")"

            else:
                output_string += str(pos) + ", "
                counter += 1


        gap_counter = len(total_gap_positions)
        chromosome = chromosome.replace("_", "")

        #print(chromosome)
        #print(len(total_gap_positions))
        output.write(str(chromosome) + "\t" + str(total_restr_sites)+ "\t" + str(gap_counter) + "\t" + str(largest_gap_start) + "\t" + str(largest_gap_end) + "\t" + str(largest_gap_len) + "\t" + str(output_string) + "\n")


parse_xmap_stats(xmap_stats, cmap)
