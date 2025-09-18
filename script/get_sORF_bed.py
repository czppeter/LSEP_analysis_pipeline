
f = open('final_nocoding_trans.gtf', 'r')
f1 = open('all_sORF.txt', 'r')
f2 = open('predict_sORF.bed', 'w')

ts_id_pos = {}
for line in f:
    line_list = line.strip().split('\t')
    feature = line_list[2]
    start = line_list[3]
    end = line_list[4]
    chr_name = line_list[0]
    strand = line_list[6]
    if strand == ".":
        strand = "+"

    message = line_list[8]
    ts_id = message.split(';')[0].split('"')[1]
    new_ts_id = ts_id + ' ' + chr_name + ' ' + strand

    pos = str(start) + '_' + str(end)

    if feature == "exon":
        if ts_id not in ts_id_pos:
            ts_id_pos[ts_id] = []
            ts_id_pos[ts_id].append(new_ts_id)
            ts_id_pos[ts_id].append(pos)
        else:
            ts_id_pos[ts_id].append(pos)


def get_gff_pos(sORF_pos, add_exon_length_list, trans_pos):
    count = 0
    for m in add_exon_length_list:
        if int(sORF_pos) > m:
            count += 1
        else:
            count = count

    sORF_exon = trans_pos[count]
    sORF_exon_start = int(sORF_exon.split('_')[0])

    if count == 0:
        sORF_gff_start_pos = int(sORF_pos) + sORF_exon_start - 1

    if count >0:

        sORF_gff_start_pos = int(sORF_pos) - add_exon_length_list[count - 1] + sORF_exon_start - 1

    return sORF_gff_start_pos,count



def get_gff_pos2(sORF_pos, add_exon_length_list, trans_pos, strand):
    if strand == "+":
        trans_length = add_exon_length_list[-1]
        if int(sORF_pos) + 3 > trans_length:
            sORF_pos = sORF_pos
        else:
            sORF_pos = int(sORF_pos) +3

    if strand == "-":
        if int(sORF_pos) - 3 < 0:
            sORF_pos = sORF_pos
        else:
            sORF_pos = sORF_pos - 3

    count = 0
    for m in add_exon_length_list:
        if int(sORF_pos) > m:
            count += 1
        else:
            count = count

    sORF_exon = trans_pos[count]
    sORF_exon_start = int(sORF_exon.split('_')[0])

    if count == 0:
        sORF_gff_start_pos = int(sORF_pos) + sORF_exon_start - 1

    if count >0:

        sORF_gff_start_pos = int(sORF_pos) - add_exon_length_list[count - 1] + sORF_exon_start - 1

    return sORF_gff_start_pos,count




def get_bed12(block_num, sORF_start_count, sORF_end_count,trans_pos,sORF_gff_start_pos,sORF_gff_end_pos):

    if block_num == 1:
        block_length_str = str(int(sORF_gff_end_pos) - int(sORF_gff_start_pos) + 1)
        block_distance_str = str(0)
    if block_num > 1:

        block_start_list = [int(sORF_gff_start_pos)-1]
        for k in range(sORF_start_count + 1, sORF_end_count + 1):
            block_start = trans_pos[k].split('_')[0]
            block_start_list.append(int(block_start)-1)

        block_end_list = []

        for j in range(sORF_start_count, sORF_end_count):
            block_end = trans_pos[j].split('_')[1]
            block_end_list.append(int(block_end)-1)
        block_end_list.append(int(sORF_gff_end_pos)-1)

        block_length_list = []
        for a in range(0, len(block_end_list)):
            block_start1 = block_start_list[a]
            block_end1 = block_end_list[a]
            block_length = int(block_end1) - int(block_start1) + 1
            block_length_list.append(str(block_length))


        new_count = 0
        for b in block_length_list:
            new_count += int(b)
        #sORF_length = int(sORF_end) - int(sORF_start) + 1

        block_distance_list = [str(0)]
        for d in block_start_list[1:]:
            distance = int(d) - int(block_start_list[0])
            block_distance_list.append(str(distance))

        block_length_str = ','.join(block_length_list)
        block_distance_str = ','.join(block_distance_list)

    return block_length_str, block_distance_str






for line in f1:
    line_list = line.strip().split('>')[1].split(' ')
    sORF_id = line_list[0]
    sORF_start = line_list[1][1:]
    sORF_end = line_list[3][:-1]

    #print(line_list,sORF_start,sORF_end)

    new_trans_id = sORF_id.split('_')[0]
    trans_pos = ts_id_pos[new_trans_id][1:]
    message = ts_id_pos[new_trans_id][0]

    strand = message.split(' ')[2]
    message_list = message.split(' ')

    ts = message_list[0]
    chr_name = message_list[1]


   # print(trans_pos, message,sORF_start,sORF_end,strand)

    length_list = []
    for i in trans_pos:
        i_start = int(i.split('_')[0])
        i_end = int(i.split('_')[1])
        #print(i_start,i_end)
        length = (i_end - i_start) + 1
       # print(length)
        length_list.append(length)

    exon_num = len(length_list)
    add_exon_length_list = []
    count = 0
    #print(length_list)
    for n in range(0, exon_num):
        count += length_list[n]
        add_exon_length_list.append(count)
  #  print(add_exon_length_list)

    # count = 0
    # for m in add_exon_length_list:
    #     if int(sORF_start) > m:
    #         count +=1
    #     else:
    #         count = count
    #
    # sORF_exon = trans_pos[count]
    # sORF_exon_start = int(sORF_exon.split('_')[0])
    #
    # sORF_gff_start_pos = int(sORF_start) - add_exon_length_list[count-1] + sORF_exon_start -1
    # print(sORF_gff_start_pos)
    if strand == "+":
        sORF_gff_start_pos = get_gff_pos(sORF_start, add_exon_length_list, trans_pos)[0]
        sORF_gff_end_pos = get_gff_pos2(sORF_end, add_exon_length_list, trans_pos, strand)[0]
        sORF_start_count = get_gff_pos(sORF_start, add_exon_length_list, trans_pos)[1]
        sORF_end_count = get_gff_pos2(sORF_end, add_exon_length_list, trans_pos, strand)[1]

        block_num = sORF_end_count - sORF_start_count + 1

        # if block_num > 1:
        #
        #     block_start_list = [sORF_gff_start_pos]
        #     for k in range(sORF_start_count+1, sORF_end_count+1):
        #         block_start = trans_pos[k].split('_')[0]
        #         block_start_list.append(block_start)
        #
        #     block_end_list = []
        #
        #     for j in range(sORF_start_count,sORF_end_count):
        #
        #         block_end = trans_pos[j].split('_')[1]
        #         block_end_list.append(block_end)
        #     block_end_list.append(sORF_gff_end_pos)
        #
        #     block_length_list = []
        #     for a in range(0,len(block_end_list)):
        #         block_start1 = block_start_list[a]
        #         block_end1 = block_end_list[a]
        #         block_length = int(block_end1) - int(block_start1) + 1
        #         block_length_list.append(block_length)
        #
        #     new_count = 0
        #     for b in block_length_list:
        #         new_count += b
        #     sORF_length = int(sORF_end) - int(sORF_start) + 1
        #
        #
        #     block_distance_list = [0]
        #     for d in block_start_list[1:]:
        #         distance = d - block_distance_list[0] + 1
        #         block_distance_list.append(distance)
        #
        #     print(block_start_list, block_end_list,block_length_list,new_count,sORF_length,)

        block_length_str,block_distance_str = get_bed12(block_num, sORF_start_count, sORF_end_count,trans_pos,sORF_gff_start_pos,sORF_gff_end_pos)




        #print(sORF_gff_start_pos, sORF_gff_end_pos, block_num,sORF_start_count,sORF_end_count)



    if strand == "-":

        new_sORF_start = add_exon_length_list[-1] - int(sORF_end) + 1
        new_sORF_end = add_exon_length_list[-1] - int(sORF_start) + 1

        print(new_sORF_start,new_sORF_end)

        sORF_gff_start_pos = get_gff_pos2(new_sORF_start, add_exon_length_list, trans_pos, strand)[0]
        sORF_gff_end_pos = get_gff_pos(new_sORF_end, add_exon_length_list, trans_pos)[0]

        sORF_start_count = get_gff_pos2(new_sORF_start, add_exon_length_list, trans_pos, strand)[1]
        sORF_end_count = get_gff_pos(new_sORF_end, add_exon_length_list, trans_pos)[1]

        block_num = sORF_end_count - sORF_start_count + 1

        block_length_str, block_distance_str = get_bed12(block_num, sORF_start_count, sORF_end_count, trans_pos,
                                                         sORF_gff_start_pos, sORF_gff_end_pos)


    #print(sORF_gff_start_pos,sORF_gff_end_pos,block_num)
    #print(add_exon_length_list[count])

    all_message = chr_name + '\t' + str(sORF_gff_start_pos-1) + '\t' + str(sORF_gff_end_pos) + '\t' + sORF_id + '\t' + "sORF" + '\t' + strand + '\t' + str(sORF_gff_start_pos-1) + '\t'  + \
                            str(sORF_gff_end_pos) + '\t' + "0" + '\t' + str(block_num) + "\t" + block_length_str + "\t" + block_distance_str

    f2.write(all_message + "\n")

    print(chr_name,sORF_gff_start_pos,sORF_gff_end_pos,ts,'sORF',strand, sORF_gff_start_pos,sORF_gff_end_pos,'0', block_num,block_length_str,block_distance_str)



















