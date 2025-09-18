
f = open('all_predict_sORF.bed', 'r')

f2 = open('predict_rm_overlap_sORF.bed','w')
all_dic = {}
for line in f:
    line_list = line.strip().split('\t')
    pre_list = line_list[0:3]
    af_list = line_list[4:]
    print(af_list)
    ts_id = line_list[3]
    key_id = '_'.join(pre_list+af_list)
    if key_id not in all_dic:
        all_dic[key_id] = []
        all_dic[key_id].append(ts_id)
    else:
        all_dic[key_id].append(ts_id)

for key,value in all_dic.items():
    key_list = key.split('_')
    af_list1 = key_list[:3]
    af_list2 = key_list[3:]

    ts_id1 = value[0]

    new_message = '\t'.join(af_list1) + '\t' + ts_id1 + '\t' + '\t'.join(af_list2)
    f2.write(new_message + '\n')




















