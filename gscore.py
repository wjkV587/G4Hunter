import re
import copy
import json
from flask import Flask
from flask import request
import requests
app = Flask(__name__)


#转化为分值数组
def translate(str):
    str_s = str
    score = []
    score_per = []
    score_per_c = []

    #处理c
    str_c = str
    str_c = re.sub(r'[GTA]',' ',str_c);
    str_is_c = re.sub(' +',' ',str_c).split(' ');
    #print(str_is_c)

    #处理G
    str = re.sub(r'[C]','1',str);
    str = re.sub(r'[TA]','0',str);
    str_temp = str
    #str = ['GGG0101GG0GGGGG01010GGG00G1GGG100GG']  
    str = re.sub(r'[01]',' ',str);
    str_not_s = re.sub(' +',' ',str).split(' ');
    #print(str_not_s)

    for x in str_not_s:
        if len(x)>0:
            score_per.append(len(x)) # score_per = [3, 2, 5, 3, 1, 3, 2]
    for x in str_is_c:
        if len(x)>0:
            score_per_c.append(len(x)) 

    #print(score_per_c)
    i = 0
    j = 0
    k = 0
    while i<len(str):
        w = str_s[i]
        #G:评分以G个数为准；
        # T：0分；
        # C：-1分；
        # A：0分
        if w == 'G':
            #判断有多少个G,append多少个
            n = score_per[j]
            i = i + score_per[j]          
            while n>0:
                score.append(4 if score_per[j]>3 else score_per[j])
                n = n-1
            j = j + 1
        elif w == 'C':
            #判断有多少个C。append多少个
            m = score_per_c[k]
            i = i + score_per_c[k]
            while m > 0 :
                score.append(-4 if score_per_c[k]>3 else -1*score_per_c[k])
                m = m - 1
            k = k + 1
        elif w == 'T' or w == 'A':
            score.append(0)
            i = i+1

    #str_sp = re.sub(' +',' ',str);
    #print(score)
    return score

#递归去重叠
def deduplication(data_view_overlaps):
    #print(data_view_overlaps)
    d = data_view_overlaps #列表有序
    #每次合并两个序列
    i = 0
    while i<len(d)-1:
        l_end = d[i]['end']
        r_start = d[i+1]['start']

        
        if l_end > r_start:
           
            d[i]['end']= d[i+1]['end']
            #合并序列seq
            
            seq_2 = d[i+1]['sequence'][l_end-r_start+1 : d[i+1]['end']-d[i+1]['start']+1]
            
            seq_1 = d[i]['sequence']
            d[i]['sequence'] = seq_1+seq_2
            arr = translate(d[i]['sequence'])
            d[i]['G4Hscore'] = sum(arr)/len(arr) #更新GHscore值
            d.remove(d[i+1])
            i = i - 1
        i = i + 1
    return d 
#计算序列串的GHscore
def calculate_score(str,windowed_value,threshold):
    #str:基因序列
    #windowed_value:可选值，滑动窗值
    #threshold:可选值，筛选

    G4Hscore = []
    #滑动窗口,l = 0,r = windowed_value
    l = 0
    array = translate(str)
    #print(array)

    #计算（2）
    while l < (len(array)-windowed_value+1):
        i = l
        sum_window = 0
        while i<l+windowed_value:
            #计算窗口平均值
            sum_window = sum_window + array[i]
            i = i+1
        l = l+1
        G4Hscore.append(sum_window/windowed_value)
    #print(G4Hscore)

    #计算（3）（4）  
    data_json = {}
    data = [] #定义返回字典
    dic = {}
    #筛选窗口段
    G4Hscore_list = []
    score_per = []
    flag = 1
    for i,x in enumerate(G4Hscore):
        
        if x > threshold:
            if flag == 1:            
                start = i+1
                flag = 0
            score_per.append(x)
        else:
            flag = 1
            if len(score_per)>0:
                G4Hscore_list.append(score_per)
                dic['start']=start
                dic['end'] = i+windowed_value-1
                dic['sequence']=str[start-1:i+windowed_value-1]
                #dic['G4Hscore(windows=25)']=score_per
                array = translate(dic['sequence'])
                dic['G4Hscore'] = sum(array)/len(array)
                data.append(dic)
                dic = {}
            score_per = []
    data_json['data_view']=data
    data_copy = copy.deepcopy(data)

    

    #print(G4Hscore_list)
    #去重叠
    data_copy = deduplication(data_copy)

    #序列首尾规范限制
    for x in data_copy:
        seq = x['sequence']
        l = 0
        r = len(seq)-1
        while l<len(seq):
            if seq[l]!= 'G':
                l = l + 1
            else:
                break
        while r>0:
            if seq[r]!= 'G':
                r = r - 1
            else:
                break
        x["sequence"] = seq[l:r+1]
        x["start"] = x['start']+l
        x["end"] = x['end'] - (len(seq)-1-r)
        arr = translate(x['sequence'])
        x["G4Hscore"] = sum(arr)/len(arr)

    #overlaps 去重叠

    data_json["data_view_overlaps"] = data_copy
    return data_json;

#定义Restful接口，返回格式JSON
@app.route('/getG4hunter',methods = ['POST','GET'])
def index():
    sequence = None
    if request.method == 'POST':
        data = request.get_data()
        data = data.decode("utf-8")
        json_data = json.loads(data)
        sequence = json_data['sequence']
        windowed_value = int(json_data['windowed_value'])
        threshold = float(json_data['threshold'])
        res = calculate_score(sequence,windowed_value,threshold)
        res = json.dumps(res)
        return str(res)

#def work():
if __name__ == '__main__':
    #str = 'GGGTCTCGGTGGGGGTCACTGGGTTGCGGGCTTGG';
    #str = 'GGGTTAGGGTTAGGGTTAGGG'
    #str = 'TGAGGGTGGGTAGGGTGGGTAA'
    #str = 'GAGGACAAGGAGGTGCGAGGAAAGGGGTTGGGGGATGGTCCCACAGGCAGCCACACCTGAGGCGTGGGCGGCCGGTAGGAGCTGGGGGAGGGCGGGGAGAAGAGGGGTTTCTGTGTGTAGA'
    #calculate_score(str,25,1.5)
    #array = translate(str)
    #print(sum(array)/len(array))
    #d = [{'sequence': 'GGAGGTGCGAGGAAAGGGGTTGGGGGATGG', 'start': 9, 'end': 38, 'G4Hscore': 1.7666666666666666}, {'sequence': 'GGGCGGCCGGTAGGAGCTGGGGGAGGG', 'start': 66, 'end': 92, 'G4Hscore': 1.6666666666666667}, {'sequence': 'GGCCGGTAGGAGCTGGGGGAGGGCGGGGAGAAGAGGGGTTTCTGTG', 'start': 70, 'end': 115, 'G4Hscore': 1.5434782608695652}]
    #print(deduplication(d));
    app.run(host='0.0.0.0',port=8080)
    