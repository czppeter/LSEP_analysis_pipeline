

import pandas as pd
import matplotlib.pyplot as plt



def get_all_count_dic(count_file):
    f = open(count_file, 'r')
    all_count_dic = {}
    for line in f:
        line_list = line.strip().split('\t')
        chr_name = line_list[0]
        pos_count_list = line_list[1].split(' ')
        pos = pos_count_list[0]
        count = pos_count_list[1]
        pos_message = chr_name + '_' + pos
        all_count_dic[pos_message] = int(count)
    return all_count_dic


def get_sORF_flank_pos(pos_file,all_count_dic):
    f2 = open(pos_file, 'r')
    up_down_count_dic = {}
    for line2 in f2:
        line2_list = line2.strip().split('\t')
        chr_name2 = line2_list[0]
        sORF_name = line2_list[1]
        pos_list = line2_list[5].split(',')
        up_down_count_list = []
        for i in pos_list:
            pos_message2 = chr_name2 + '_' + str(i)
            try:
                count2 = all_count_dic[pos_message2]
            except KeyError:
                count2 = 0
            up_down_count_list.append(count2)

        up_down_count_dic[sORF_name] = up_down_count_list

    return up_down_count_dic


def get_filter_sORF(pick_sORF,up_down_count_dic):
    f3 = open(pick_sORF,'r')
    pick_sORF = []
    for line in f3:
        sORF_name3 = line.strip()
        pick_sORF.append(sORF_name3)

    filtered_dic = {k: up_down_count_dic[k] for k in pick_sORF if k in up_down_count_dic}
    return filtered_dic


def get_pos_count(count_dic):
    data = pd.DataFrame(count_dic)
    data['Col_sum'] = data.apply(lambda x: x.sum(), axis=1)

    all_count = data['Col_sum']
    count_list = []
    ribosome_phase_list = []
    for a in range(0, len(all_count)):
        i_number = all_count[a]
        count_list.append(all_count[a])
        ribosome_phase_list.append(i_number)
    return ribosome_phase_list


count_27_dic = get_all_count_dic('Mock-48h.count.27.txt')
count_28_dic = get_all_count_dic('Mock-48h.count.28.txt')
count_29_dic = get_all_count_dic('Mock-48h.count.29.txt')
count_30_dic = get_all_count_dic('Mock-48h.count.30.txt')
count_31_dic = get_all_count_dic('Mock-48h.count.31.txt')
count_32_dic = get_all_count_dic('Mock-48h.count.32.txt')
count_33_dic = get_all_count_dic('Mock-48h.count.33.txt')
count_34_dic = get_all_count_dic('Mock-48h.count.34.txt')
count_35_dic = get_all_count_dic('Mock-48h.count.35.txt')
count_36_dic = get_all_count_dic('Mock-48h.count.36.txt')



up_down_count_dic_27 = get_sORF_flank_pos('sORF_flank_pos.txt',count_27_dic)
up_down_count_dic_28 = get_sORF_flank_pos('sORF_flank_pos.txt',count_28_dic)
up_down_count_dic_29 = get_sORF_flank_pos('sORF_flank_pos.txt',count_29_dic)
up_down_count_dic_30 = get_sORF_flank_pos('sORF_flank_pos.txt',count_30_dic)
up_down_count_dic_31 = get_sORF_flank_pos('sORF_flank_pos.txt',count_31_dic)
up_down_count_dic_32 = get_sORF_flank_pos('sORF_flank_pos.txt',count_32_dic)
up_down_count_dic_33 = get_sORF_flank_pos('sORF_flank_pos.txt',count_33_dic)
up_down_count_dic_34 = get_sORF_flank_pos('sORF_flank_pos.txt',count_34_dic)
up_down_count_dic_35 = get_sORF_flank_pos('sORF_flank_pos.txt',count_35_dic)
up_down_count_dic_36 = get_sORF_flank_pos('sORF_flank_pos.txt',count_36_dic)

filtered_dic_27 = get_filter_sORF('Identity_non_gene_pos_sORF.txt', up_down_count_dic_27)
filtered_dic_28 = get_filter_sORF('Identity_non_gene_pos_sORF.txt', up_down_count_dic_28)
filtered_dic_29 = get_filter_sORF('Identity_non_gene_pos_sORF.txt', up_down_count_dic_29)
filtered_dic_30 = get_filter_sORF('Identity_non_gene_pos_sORF.txt', up_down_count_dic_30)
filtered_dic_31 = get_filter_sORF('Identity_non_gene_pos_sORF.txt', up_down_count_dic_31)
filtered_dic_32 = get_filter_sORF('Identity_non_gene_pos_sORF.txt', up_down_count_dic_32)
filtered_dic_33 = get_filter_sORF('Identity_non_gene_pos_sORF.txt', up_down_count_dic_33)
filtered_dic_34 = get_filter_sORF('Identity_non_gene_pos_sORF.txt', up_down_count_dic_34)
filtered_dic_35 = get_filter_sORF('Identity_non_gene_pos_sORF.txt', up_down_count_dic_35)
filtered_dic_36 = get_filter_sORF('Identity_non_gene_pos_sORF.txt', up_down_count_dic_36)



pos_count_list_27 = get_pos_count(filtered_dic_27)
pos_count_list_28 = get_pos_count(filtered_dic_28)
pos_count_list_29 = get_pos_count(filtered_dic_29)
pos_count_list_30 = get_pos_count(filtered_dic_30)
pos_count_list_31 = get_pos_count(filtered_dic_31)
pos_count_list_32 = get_pos_count(filtered_dic_32)
pos_count_list_33 = get_pos_count(filtered_dic_33)
pos_count_list_34 = get_pos_count(filtered_dic_34)
pos_count_list_35 = get_pos_count(filtered_dic_35)
pos_count_list_36 = get_pos_count(filtered_dic_36)



all_lists = [pos_count_list_27, pos_count_list_28, pos_count_list_29, pos_count_list_30, pos_count_list_31, pos_count_list_32,
             pos_count_list_33,pos_count_list_34,pos_count_list_35,pos_count_list_36]

ribosome_phase_list = [sum(values) for values in zip(*all_lists)]

print(ribosome_phase_list)

plt.figure(figsize=(4, 4))
x_list = []
for b in range(-31, 30):
    x_list.append(b)
plt.bar(x_list, ribosome_phase_list, color=['r', 'g', 'b'])
plt.title('sORF_start')

plt.savefig('mock_start_flank.pdf')

#plt.savefig('mock_end_flank.pdf')

plt.show()