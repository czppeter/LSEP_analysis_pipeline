
#f = open('PD_and_maxquant_sORF_count.txt', 'r')
#f1 = open('centromere_pos.txt', 'r')
#f3 = open('relate_distance_count.txt', 'w')

f = open('PD_and_maxquant_sORF_count_filter.txt', 'r')
f1 = open('centromere_pos.txt', 'r')
f3 = open('relate_distance_count_filter.txt', 'w')



chr_length = {'Chr1':43270923, 'Chr2':35937250, 'Chr3':36413819, 'Chr4':35502694, 'Chr5':29958434, 'Chr6':31248787,
              'Chr7':29697621, 'Chr8':28443022, 'Chr9':23012720, 'Chr10':23207287, 'Chr11':29021106, 'Chr12':27531856}


centro_dic = {}
for line in f1:
    line_list = line.strip().split('\t')
    chr_name = line_list[0]
    ce_start = line_list[1]
    ce_end = line_list[2]
    ce = ce_start + '_' +ce_end
    centro_dic[chr_name] = ce



all_count_dic = {}
for i in f:
    i_list = i.strip().split('\t')
    chr_name = i_list[0]

    if chr_name not in all_count_dic:
        all_count_dic[chr_name] = []
        all_count_dic[chr_name].append(i_list)
    else:
        all_count_dic[chr_name].append(i_list)



distance_dic = {}
for key,vaue in all_count_dic.items():
    for line_list2 in vaue:
       # print(line_list2)
        chr_name = line_list2[0]
        start = int(line_list2[1])
        end = int(line_list2[2])
        trans_count = line_list2[3]
        notrans_count = line_list2[4]
        all_count = trans_count + '_' + notrans_count

        ce_start = int(centro_dic[chr_name].split('_')[0])
        ce_end = int(centro_dic[chr_name].split('_')[1])

       # print(ce_start,start)
        if chr_name not in distance_dic:
            distance_dic[chr_name] = {}

     #   print(start)
        if ce_start >= start:
            distance = int(ce_start) - int(start)
            relate_distance = distance / int(ce_start)


            if relate_distance not in distance_dic[chr_name]:
                distance_dic[chr_name][relate_distance] = []
                distance_dic[chr_name][relate_distance].append(all_count)

            else:
                distance_dic[chr_name][relate_distance].append(all_count)

        if start > ce_start and start <= ce_end:

            distance = 0
            relate_distance = 0


            if relate_distance not in distance_dic[chr_name]:
                distance_dic[chr_name][relate_distance] = []
                distance_dic[chr_name][relate_distance].append(all_count)

            else:
                distance_dic[chr_name][relate_distance].append(all_count)

        if start > ce_end:

            distance = int(start) - int(ce_end)
            relate_distance = distance / (chr_length[chr_name] - ce_end)

         #   print(ce_start,start,distance)

            if relate_distance not in distance_dic[chr_name]:
                distance_dic[chr_name][relate_distance] = []
                distance_dic[chr_name][relate_distance].append(all_count)

            else:
                distance_dic[chr_name][relate_distance].append(all_count)

for key,value in distance_dic.items():
    for m,n in value.items():

        count1_sum = 0
        count2_sum = 0

        for i in n:
            count1= int(i.split('_')[0])
            count2 = int(i.split('_')[1])

            count1_sum += count1
            count2_sum += count2
        print(key, m, count1_sum, count2_sum, n)
        all_str = key + '\t' + str(m) + '\t' +str(count1_sum)+ '\t' + str(count2_sum)
        f3.write(all_str + '\n')

