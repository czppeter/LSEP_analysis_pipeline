
import numpy as np
#f = open('final_nocoding_trans_rmChrm.gtf', 'r')
#f1 = open('predict_rm_overlap_sORF.bed', 'r')
#f2 = open('all_sORF_related_pos.txt', 'w')


#f = open('final_nocoding_trans_rmChrm.gtf', 'r')
#f1 = open('PD_sORF_filter.bed12', 'r')
#f2 = open('PD_sORF_related_pos_filter.txt', 'w')

f = open('final_nocoding_trans_rmChrm.gtf', 'r')
f1 = open('Maxquant_sORF_filter.bed12', 'r')
f2 = open('maxquant_sORF_related_pos_filter.txt', 'w')





lncrna_pos_dic = {}
for line in f:
    line_list = line.strip().split('\t')
    chr_name, start, end, feature, strand = line_list[0], line_list[3], line_list[4], line_list[2], line_list[6]
    message = line_list[8]
    lncrna_id = message.split(';')[0].split('"')[1]
    if feature == "exon":
        pos_message = start + '_' + end
        lncrna_id_strand = lncrna_id + '_' + strand
        if lncrna_id_strand not in lncrna_pos_dic:
            lncrna_pos_dic[lncrna_id_strand] = []
            lncrna_pos_dic[lncrna_id_strand].append(pos_message)

        else:
            lncrna_pos_dic[lncrna_id_strand].append(pos_message)


##获取每一个lnRNA 的长度
lnc_length_dic = {}
for key, value in lncrna_pos_dic.items():
    lnc_length = 0

    for i in value:
        start = i.split('_')[0]
        end = i.split('_')[1]
        length = int(end) - int(start) + 1
        lnc_length += length
    lnc_length_dic[key] = lnc_length


##获取每一个sORF 的起始位点在基因组上
sORF_start_dic = {}
for line in f1:
    line_list = line.strip().split('\t')
    strand = line_list[5]
    sROF_name = line_list[3]
    if strand == "+":
        start = line_list[1]
    if strand == "-":
        start = line_list[2]

    message = sROF_name + ' ' + strand
    sORF_start_dic[message] = start


##获取sORF 相对在lncRNA 上的起始位置
sORF_start_pos_dic = {}
for key, value in sORF_start_dic.items():
    message_list = key.split(' ')
    print(message_list)
    lncrna_id = message_list[0].split('_')[0]
    strand = message_list[1]

    new_key = lncrna_id + '_' + strand
    if new_key not in lncrna_pos_dic:
        continue
    else:
        lncrna_pos_list = lncrna_pos_dic[new_key]
    count = 0
    length = ''
    region_length_list = []
    all_length = lnc_length_dic[new_key]
    new_count = ''
    print(key,value,lncrna_pos_list)
    for i in lncrna_pos_list:
        start = int(i.split('_')[0])
        end = int(i.split('_')[1])
        count += 1

        region_length = end - start + 1
        region_length_list.append(region_length)
        if int(value) + 1 >= int(start) and int(value) <= int(end):
            length = int(value) - int(start) + 1
            new_count = count
            break

    start_pos = ''
    if new_count == 1 and strand == "+":
        start_pos = length + 1
        #print(key, start_pos)

    if new_count >1 and strand == "+":
        length1 = sum(region_length_list[:new_count-1])
        start_pos = length1 + length + 1
       # print(key,start_pos,new_count,region_length_list)

    if new_count == 1 and strand == "-":
        start_pos = int(all_length) - length + 1
        #print(key, start_pos)

    if new_count>1 and strand == "-":
        length1 = sum(region_length_list[:new_count - 1]) + length
        start_pos = int(all_length) -length1 + 1
        #print(key, start_pos)


    sORF_start_pos_dic[key] = start_pos



#将lncRNA 归一化到50，并计算sORF所在的相对位置


scale_pos_value_list = []
for key,value in sORF_start_pos_dic.items():
    message_list = key.split(' ')
    ln_name = message_list[0].split('_')[0]
    strand = message_list[1]
    new_key = ln_name + '_' + strand
    lnc_length = int(lnc_length_dic[new_key])
   # print(lnc_length, value)

    scale_factor = 50/lnc_length


    pos = round(scale_factor*int(value))

    pos1 = scale_factor*int(value)
    print(lnc_length, value, pos, pos1)

    value_list = [0]*50
    print(value_list)
    value_list[pos-1] = 1
    print(value_list)

    scale_pos_value_list.append(value_list)

arr1 = np.array(scale_pos_value_list)

print(arr1)
sum_list = list(arr1.sum(axis=0))
sum_list2 = sum_list/sum(sum_list)

print(sum_list,sum_list2)

for m,n in zip(sum_list, sum_list2):
    f2.write(str(m) + '\t' + str(n) + '\t' "all" + '\n')



















