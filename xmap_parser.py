

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
    output.write("REF_START\tREF_END\tREF_LEN\tQUERY_START\tQUERY_END\tQUERY_LEN\tGAP_LEN\tRESTR_SITES\tGAP_RESTR_SITES\n")


    xmap_values = open(xmap).readlines()
    xmap_stats = []

    for value in xmap_values:
        gap_length = 0.0

        #Skip comments in header
        if value.startswith('#'):
            continue

        else:
            value_strip = value.rstrip()
            value_split = value_strip.split("\t")

            #Extract all info from xmap
            xmap_ID = value_split[0]
            qry_contig_ID = value_split[1]
            ref_contig_ID = value_split[2]
            qry_start_pos = value_split[3]
            qry_end_pos = value_split[4]
            ref_start_pos = value_split[5]
            ref_end_pos = value_split[6]
            orientation = value_split[7]
            confidence = value_split[8]
            hitenum = value_split[9]
            qry_len = value_split[10]
            ref_len = value_split[11]
            label_channel = value_split[12]
            alignment = value_split[13]

            #find length of aligned fragments
            ref_len = float(ref_end_pos) - float(ref_start_pos)
            qry_len = float(qry_end_pos) - float(qry_start_pos)

            if ref_len > qry_len:
                gap_length = ref_len - qry_len

            else:
                gap_length = qry_len - ref_len

        # call function to find restr sites in region and gaps from cmap
        # Restriction site info is appended to xmap stats, go through and write to output
            xmap_stats = find_restr_sites(cmap, ref_contig_ID, ref_start_pos, ref_end_pos, ref_len, qry_start_pos, qry_end_pos, qry_len, gap_length)
            for stats in xmap_stats:
                counter = 0
                for i in stats:
                    if counter == len(stats)-1:
                        output.write(str(i) + "\n")
                    else:
                        output.write(str(i) + "\t")
                    counter += 1


###################### Find_Restr_Sites ###########################
# Input: cmap file, reference contig id, reference start, reference end,
# reference length, query start, query end, query end, and gap length
# returns statistics on the reference and query fragments as well as
# the restriction sites on the fragments and gaps in alignment
###################################################################

def find_restr_sites(cmap, ref_id, ref_start, ref_end, ref_len, qry_start, qry_end, qry_len, gap_length):

    counter = 0
    gap_counter = 0
    num_restr_sites = 0

    output_list = []

    cmap_values = open(cmap).readlines()
    for cmap in cmap_values:

        if cmap.startswith('#'):
            continue

        else:

            #Extract info from cmap
            cmap_strip = cmap.rstrip()
            cmap_split = cmap_strip.split("\t")

            cmap_ID = cmap_split[0]
            contig_len = cmap_split[1]
            num_sites = cmap_split[2]
            site_ID = cmap_split[3]
            label_channel = cmap_split[4]
            position = cmap_split[5]

            #Find restriciton sites for the correct contig
            if ref_id  == cmap_ID:
                new_contig_len = contig_len
                num_restr_sites += 1

                if float(position) >= float(ref_start) and float(position) <= float(ref_end):
                    counter += 1
                else:
                    gap_counter += 1
    #print("cmap ID: " + ref_id + "\ncontig length: " + new_contig_len + "\n" + str(ref_start) + "\t" + str(ref_end) + "\t" + str(ref_len) + "\t" + str(qry_start) + "\t" + str(qry_end) + "\t" + str(qry_len) + "\t" + str(counter) + "\t" + str(gap_counter) + "\n")
    output_list.append([ref_start, ref_end, ref_len, qry_start, qry_end, qry_len, gap_length, counter, gap_counter])
    return output_list

###################### Commands #########################

parse_xmap("/Users/madisonzwink/xmap_test/NC_010473_mock_scaffolds_to_ESCH_COLI_1_2015_000_STRICT_T_150_REFINEFINAL1.xmap", "/Users/madisonzwink/xmap_test/NC_010473_mock_scaffolds_BspQI.cmap")
