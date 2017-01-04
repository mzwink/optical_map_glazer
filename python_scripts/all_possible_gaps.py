
##FIX GAP LENGTH INFO (START/END)


######################## Parse_Xmap ###########################
# Input: xmap file, cmap file and extracts info from files,
# returns statistics on the reference and query fragments as well as
# the restriction sites on the fragments and gaps in alignment
################################################################

def parse_xmap(xmap, cmap, cmap_key):

    #xmap format:
    #XmapEntryID QryContigID RefContigID QryStartPos QryEndPos RefStartPos
    # RefEndPos Orientation Confidence  HitEnum QryLen  RefLen  LabelChannel Alignment

    #Output that will summarize xmap and cmap statistics
    output_file = xmap.replace(".xmap", "gaps.txt")
    output = open(output_file, 'w')
    output.write("REF_ID\tREF_START\tREF_END\tREF_LEN\tRESTR_SITES\tQRY_ID\tALL_GAPS\tLARGEST_GAP_START\tLARGEST_GAP_END\tGAP_LEN\tGAP_RESTR_SITES\n")

    #Initialize variables that will hold xmap stats
    xmap_values = open(xmap).readlines()
    xmap_stats = []
    #Call function to make a dictionary of chromosome to reference id conversions
    ref_to_chromosome = chr_scaffold_info(cmap_key)
    chr_keys = ref_to_chromosome.keys()
    #print(chr_keys)
    ref_id_list = {}

    #Go through xmap for each chromosome, append all info to lists for analysis
    for key in chr_keys:
        for ref_id_chr in ref_to_chromosome[key]:
            chromosome = key.replace("_", "")
            print(chromosome)
            #print(ref_id_chr)

            if ref_id_chr in ref_id_list.keys():
                ref_id_list[ref_id_chr].append(0)
            else:
                ref_id_list[ref_id_chr] = [0]

    

    #######################################################







            aligment_list = []
            gap_list = []
            ref_coordinate_list = {}
            ref_min = 0
            ref_max = 0



            for value in xmap_values:

                if value.startswith('#'):
                    continue

                else:
                    value_strip = value.rstrip()
                    value_split = value_strip.split("\t")


                    ref_contig_ID = int(value_split[2])
                    if ref_contig_ID in ref_id_list.keys():

                        ref_contig_ID = value_split[2]
                        ref_start_pos = float(value_split[5])
                        ref_end_pos = float(value_split[6])
                        ref_len = value_split[11]
                        alignment = value_split[13]

                        aligment_list.append(alignment)


            #Separate alignment coordinates

                    #print(chromosome)
                    #print(ref_min)
                    #print(ref_max)

                    for align in aligment_list:
                        #print(align)
                        align_split = align.split(")(")

                        #Get range of reference alignment (min -> max)
                        #Use this range to find possible alignment gaps
                        #ref_min = 0
                        #ref_max = 0
                        for coordinate in align_split:
                            coordinate = coordinate.replace("(", "")
                            coordinate = coordinate.replace(")", "")

                            ref_coordinate = coordinate.split(",")
                            ref_coordinate = int(ref_coordinate[0])
                            if ref_min == 0:
                                ref_min = ref_coordinate
                            else:
                                if ref_coordinate < ref_min:
                                    ref_min = ref_coordinate

                            if ref_coordinate > ref_max:
                                ref_max = ref_coordinate
                        #Add all reference coordinates to a dictionary for easy look up
                            if int(ref_coordinate) in ref_coordinate_list.keys():
                                ref_coordinate_list[int(ref_coordinate)].append(0)
                            else:
                                ref_coordinate_list[int(ref_coordinate)] = [0]

                        restr_site = 0
                        gap_counter = 0
                        continuous_gap_list = []
                        if chromosome == "chrI":
                            print(chromosome)
                            print(ref_min)
                            print(ref_max)

                        #Go through the reference range to find gaps
                        #if value in range not in ref coordinate, add value to gap list
                        for i in range(ref_min, ref_max+1):

                            if i in ref_coordinate_list.keys():
                                restr_site +=1
                            else:
                                gap_counter += 1
                                gap_list.append(i)



def chr_scaffold_info(chr_key):#, chr_num):

    chr_id_key = open(chr_key).readlines()
    chr_dict = {}

    for chr_num in chr_id_key:
        if chr_num.startswith("#") or chr_num.startswith("CompntId"):
            continue
        else:
            chr_num_strip = chr_num.rstrip()
            chr_num_split = chr_num_strip.split("\t")
            ref_id = int(chr_num_split[0])
            chr_id = chr_num_split[1]

            for num in range(1,22):

                if num == 1: chromosome = "chrI_"
                if num == 2: chromosome = "chrII_"
                if num == 3: chromosome = "chrIII_"
                if num == 4: chromosome = "chrIV_"
                if num == 5: chromosome = "chrV_"
                if num == 6: chromosome = "chrVI_"
                if num == 7: chromosome = "chrVII_"
                if num == 8: chromosome = "chrVIII_"
                if num == 9: chromosome = "chrIX_"
                if num == 10: chromosome = "chrX_"
                if num == 11: chromosome = "chrXI_"
                if num == 12: chromosome = "chrXII_"
                if num == 13: chromosome = "chrXIII_"
                if num == 14: chromosome = "chrXIV_"
                if num == 15: chromosome = "chrXV"
                if num == 16: chromosome = "chrXVI_"
                if num == 17: chromosome = "chrXVII_"
                if num == 18: chromosome = "chrXVIII_"
                if num == 19: chromosome = "chrXIX_"
                if num == 20: chromosome = "chrXX_"
                if num == 21: chromosome = "chrXXI_"

                if chr_id.startswith(chromosome):


                    if chromosome in chr_dict.keys():
                        chr_dict[chromosome].append(ref_id)

                    else:
                        chr_dict[chromosome] = [ref_id]

    return chr_dict

######################################################################################
############################### Path to files ########################################

cmap_key = "/Volumes/MW_18TB/Madison_Zwink/optical_map/glazer_chr_scaffolds_BspQI_key.txt"
xmap = "/Volumes/MW_18TB/Madison_Zwink/optical_map/ref_Nov_17/glazer_chr_scaffolds_wrap_BspQI_to_GAST_ACUL_2014_057_DEFAULT_T_REFINEFINAL1.xmap"
cmap = "/Volumes/MW_18TB/Madison_Zwink/optical_map/glazer_chr_scaffolds_wrap_BspQI.cmap"

parse_xmap(xmap, cmap, cmap_key)
#chr_scaffold_info(cmap_key)


            #Use gap list in function to find largest continuous gap in fragment
            #continuous_gap_list = find_gap_sites(gap_list, ref_min, ref_max)

            #print(continuous_gap_list)

            #position_1 = 0
            #position_2 = 0

            #if len(continuous_gap_list) > 1:

            #    cmap_info = open(cmap).readlines()

            #    for value in cmap_info:
            #        if value.startswith('#'):
            #            continue
            #        else:

            #            value_strip = value.rstrip()
            #            value_split = value_strip.split("\t")
            #            cmap_id = value_split[0]
            #            numsite_id = value_split[3]
            #            position = value_split[5]

            #            #print(ref_contig_ID)
            #            if cmap_id == ref_contig_ID:#
            #                #print(numsite_id)
            #                #print(continuous_gap_list[0])
            #                if int(numsite_id) == int(continuous_gap_list[0]):
            #                #if int(numsite_id) == int(continuous_gap_list[0]):
            #                    position_1 = float(position)
            #                elif int(numsite_id) == int(continuous_gap_list[len(continuous_gap_list)-1]):
            #                   position_2 = float(position)

            #    gap_length = position_2 - position_1
            #    gap_counter = len(continuous_gap_list)
            #    scaffold_id = ref_to_scaff[int(ref_contig_ID)]
            #    for value in scaffold_id:
            #        scaffold_id = value#.replace("[", "").replace("]", "").replace("'", "")
            #        print(scaffold_id)
            #    #output.write("REF_ID\tREF_START\tREF_END\tREF_LEN\tRESTR_SITES\tGAP_START\tGAP_END\tGAP_LEN\tGAP_RESTR_SITES\n")
            #    output.write(str(scaffold_id) + "\t" + str(ref_start_pos) + "\t" + str(ref_end_pos) + "\t" + str(ref_len) + "\t" + str(qry_contig_ID) + "\t"  + str(qry_start_pos) + "\t" + str(qry_end_pos) + "\t" + str(qry_len)
            #    + "\t" + str(restr_site) + "\t" + str(position_1) + "\t" + str(position_2) + "\t"
            #    + str(gap_length) + "\t" + str(gap_counter) + "\n")



            ######################################
            #### Now have info for xmap stats ####
            #### Use xmap id to find gap length ##
            #### from the cmap values ############
            ######################################

#def find_gap_sites(gap_list, ref_min, ref_max):

#        continuous_gap_list = []
#        past_list = []
#        current_list = []


#        if len(gap_list) > 1:
#            #print(gap_list)
#            for site in range(0,len(gap_list)-1):
#                if site == 0:
#                    if int(gap_list[site+1]) - int(gap_list[site]) == 1:
#                        current_list.append(gap_list[site])
#                        #print(current_list)

#                elif site > 0 and site < len(gap_list)-1:
#                    if int(gap_list[site+1]) - int(gap_list[site]) == 1:
#                        current_list.append(gap_list[site])

#                    else:
#                        if len(past_list) == 0:
#                            past_list = current_list
#                            current_list = []

#                        else:
#                            if len(current_list) >= len(past_list):
#                                past_list = current_list
#                                current_list = []
#                            else:
#                                 current_list = []

#                elif site == len(gap_list) -1:
#                    if int(gap_list[site]) - int(gap_list[site-1]) == 1:
#                        current_list.append(gap_list[site])


#            if len(past_list) >= len(current_list):
#                continuous_gap_list = past_list
#            else:
#                continuous_gap_list = current_list

#        return continuous_gap_list
