
##FIX GAP LENGTH INFO (START/END)


######################## Parse_Xmap ###########################
# Input: xmap file, cmap file and extracts info from files,
# returns statistics on the reference and query fragments as well as
# the restriction sites on the fragments and gaps in alignment
################################################################

def parse_xmap(xmap, cmap):

    #xmap format:
    #XmapEntryID QryContigID RefContigID QryStartPos QryEndPos RefStartPos
    # RefEndPos Orientation Confidence  HitEnum QryLen  RefLen  LabelChannel Alignment

    #Output that will summarize xmap and cmap statistics
    output_file = xmap.replace(".xmap", "_stats.txt")
    output = open(output_file, 'w')
    output.write("REF_ID\tREF_START\tREF_END\tREF_LEN\tQRY_ID\tQRY_START\tQRY_END\tQRY_LEN\tRESTR_SITES\tGAP_START\tGAP_END\tGAP_LEN\tGAP_RESTR_SITES\n")

    xmap_values = open(xmap).readlines()
    xmap_stats = []
    ref_to_scaff = chr_scaffold_info("/Users/madisonzwink/Desktop/Glazer_project/glazer_chr_scaffolds_BspQI_key.txt")

    for value in xmap_values:

        #Skip comments in header
        if value.startswith('#'):
            continue

        else:
            value_strip = value.rstrip()
            value_split = value_strip.split("\t")

            #Extract all info from xmap
            xmap_ID = int(value_split[0])
            qry_contig_ID = value_split[1]
            ref_contig_ID = value_split[2]
            qry_start_pos = float(value_split[3])
            qry_end_pos = float(value_split[4])
            ref_start_pos = float(value_split[5])
            ref_end_pos = float(value_split[6])
            qry_len = 0 #value_split[10]
            ref_len = 0 #value_split[11]
            alignment = value_split[13]
            ref_coordinate_list = {}


            if qry_start_pos > qry_end_pos:
                qry_len = qry_start_pos - qry_end_pos

            else:
                qry_len = qry_end_pos - qry_start_pos

            ref_len = ref_end_pos - ref_start_pos


            alignment_split = alignment.split(")(")

            ref_min = 0
            ref_max = 0
            for coordinate in alignment_split:
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

                if int(ref_coordinate) in ref_coordinate_list.keys():
                    ref_coordinate_list[int(ref_coordinate)].append(0)
                else:
                    ref_coordinate_list[int(ref_coordinate)] = [0]

                restr_site = 0
                gap_counter = 0
                gap_list = []
                continuous_gap_list = []


                for i in range(ref_min, ref_max+1):

                    if i in ref_coordinate_list.keys():
                        restr_site +=1
                    else:
                        gap_counter += 1
                        gap_list.append(i)



            continuous_gap_list = find_gap_sites(gap_list, ref_min, ref_max)

            #print(continuous_gap_list)

            position_1 = 0
            position_2 = 0

            if len(continuous_gap_list) > 1:

                cmap_info = open(cmap).readlines()

                for value in cmap_info:
                    if value.startswith('#'):
                        continue
                    else:

                        value_strip = value.rstrip()
                        value_split = value_strip.split("\t")
                        cmap_id = value_split[0]
                        numsite_id = value_split[3]
                        position = value_split[5]

                        #print(ref_contig_ID)
                        if cmap_id == ref_contig_ID:#
                            #print(numsite_id)
                            #print(continuous_gap_list[0])
                            if int(numsite_id) == int(continuous_gap_list[0]):
                            #if int(numsite_id) == int(continuous_gap_list[0]):
                                position_1 = float(position)
                            elif int(numsite_id) == int(continuous_gap_list[len(continuous_gap_list)-1]):
                               position_2 = float(position)

                gap_length = position_2 - position_1
                gap_counter = len(continuous_gap_list)
                scaffold_id = ref_to_scaff[int(ref_contig_ID)]
                for value in scaffold_id:
                    scaffold_id = value#.replace("[", "").replace("]", "").replace("'", "")
                    print(scaffold_id)
                #output.write("REF_ID\tREF_START\tREF_END\tREF_LEN\tRESTR_SITES\tGAP_START\tGAP_END\tGAP_LEN\tGAP_RESTR_SITES\n")
                output.write(str(scaffold_id) + "\t" + str(ref_start_pos) + "\t" + str(ref_end_pos) + "\t" + str(ref_len) + "\t" + str(qry_contig_ID) + "\t"  + str(qry_start_pos) + "\t" + str(qry_end_pos) + "\t" + str(qry_len)
                + "\t" + str(restr_site) + "\t" + str(position_1) + "\t" + str(position_2) + "\t"
                + str(gap_length) + "\t" + str(gap_counter) + "\n")



            ######################################
            #### Now have info for xmap stats ####
            #### Use xmap id to find gap length ##
            #### from the cmap values ############
            ######################################

def find_gap_sites(gap_list, ref_min, ref_max):

        continuous_gap_list = []
        past_list = []
        current_list = []


        if len(gap_list) > 1:
            #print(gap_list)
            for site in range(0,len(gap_list)-1):
                if site == 0:
                    if int(gap_list[site+1]) - int(gap_list[site]) == 1:
                        current_list.append(gap_list[site])
                        #print(current_list)

                elif site > 0 and site < len(gap_list)-1:
                    if int(gap_list[site+1]) - int(gap_list[site]) == 1:
                        current_list.append(gap_list[site])

                    else:
                        if len(past_list) == 0:
                            past_list = current_list
                            current_list = []

                        else:
                            if len(current_list) >= len(past_list):
                                past_list = current_list
                                current_list = []
                            else:
                                 current_list = []

                elif site == len(gap_list) -1:
                    if int(gap_list[site]) - int(gap_list[site-1]) == 1:
                        current_list.append(gap_list[site])


            if len(past_list) >= len(current_list):
                continuous_gap_list = past_list
            else:
                continuous_gap_list = current_list

        return continuous_gap_list

def chr_scaffold_info(chr_key):

    chr_id_key = open(chr_key).readlines()
    scaffold_id_dict = {}

    for scaff in chr_id_key:
        if scaff.startswith("#") or scaff.startswith("CompntId"):
            continue
        else:
            scaff_id = scaff.rstrip()
            scaff_split = scaff_id.split("\t")
            ref_id = int(scaff_split[0])
            chr_id = scaff_split[1]

            if scaff_id in scaffold_id_dict.keys():
                scaffold_id_dict[ref_id].append(chr_id)

            else:
                scaffold_id_dict[ref_id] = [chr_id]

    return scaffold_id_dict




parse_xmap("Refaligner_output/glazer_chr_scaffolds_wrap_BspQI_to_GAST_ACUL_2014_057_DEFAULT_T_REFINEFINAL1.xmap", "glazer_chr_scaffolds_BspQI.cmap")
