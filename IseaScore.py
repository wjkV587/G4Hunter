"""
@Author：wangjunkai

"""
import re
import copy
import json
from flask import Flask
from flask import request,make_response
from data_clean import *
import requests
from flask_cors import *  # 导入模块

app = Flask(__name__)

"""

1. I_motif() 函数为自定义算分方法
2. translate() 函数为G4结构算分方法
3. calculate_score(): I_motif算分方法的on a long sequence

"""

# function : I motif
def I_motif(str):
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
        """
        G:每个G取值=(-1)*G个数；
        T：0分；
        C：每个C取值=1*C个数
        A：0分；
        N：这里多了一个N，记为0分
        """
        if w == 'C':
            #判断有多少个C,append多少个
            n = score_per_c[j]
            i = i + score_per_c[j]          
            while n>0:
                score.append(4 if score_per_c[j]>3 else 1*score_per_c[j])
                n = n-1
            j = j + 1
        elif w == 'G':
            #判断有多少个G。append多少个
            m = score_per[k]
            i = i + score_per[k]
            while m > 0 :
                score.append(-4 if score_per[k]>3 else -1*score_per[k])
                m = m - 1
            k = k + 1
        elif w == 'T' or w == 'A' or w == 'N':
            score.append(0)
            i = i+1

    #str_sp = re.sub(' +',' ',str);
    #print('第1步：\n',score)
    return score

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
    print('第1步：\n',score)
    return score

#迭代,重叠融合
def deduplication(data_view_overlaps):
    #print(data_view_overlaps)
    d = data_view_overlaps #列表有序
    #每次合并两个序列
    i = 0
    while i<len(d)-1:
        l_end = d[i]['end']
        r_start = d[i+1]['start']

        
        if l_end > r_start:
            # 相邻序列存在可以重叠融合的部分
            d[i]['end']= d[i+1]['end']
            #合并序列seq
            
            seq_2 = d[i+1]['sequence'][l_end-r_start+1 : d[i+1]['end']-d[i+1]['start']+1]
            
            seq_1 = d[i]['sequence']
            d[i]['sequence'] = seq_1+seq_2
            arr = I_motif(d[i]['sequence'])
            d[i]['ISEAscore'] = sum(arr)/len(arr) #更新GHscore值
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
    #array = translate(str)
    array = I_motif(str)
    print('第1步：\n',array)

    #计算（2）
    while l < (len(array)-windowed_value+1):
        i = l
        sum_window = 0
        while i<l+windowed_value:
            #计算窗口平均值
            sum_window = sum_window + array[i]
            i = i+1
        l = l+1
        # 取g4score的绝对值
        average = sum_window/windowed_value if sum_window/windowed_value >=0 else -1*sum_window/windowed_value 
        G4Hscore.append(average)
    print('第2步：\n',G4Hscore)

    #计算（3）（4）  
    data_json = {}
    data = [] #定义返回字典
    dic = {}
    #筛选窗口段
    G4Hscore_list = []
    score_per = []
    flag = 1
    #print('调试=',G4Hscore)

    # 进行序列的划分，每一段序列的值都需要大于threshold值
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
                dic['end'] = i+windowed_value-1 #当前分数代表(当前字母+windowed_value-1)组成的序列取值
                dic['sequence']=str[start-1:i+windowed_value-1]
                dic['sequence_size'] = len(dic['sequence'])
                #dic['G4Hscore(windows=25)']=score_per
                array = I_motif(dic['sequence'])
                dic['ISEAscore'] = sum(array)/len(array)
                data.append(dic)
                dic = {}
            score_per = []
    data_json['data_view']=data #包含重叠
    data_copy = copy.deepcopy(data)

    #print(G4Hscore_list)

    #迭代，进行重叠融合
    data_copy = deduplication(data_copy)

    #序列首尾规范限制
    for x in data_copy:
        seq = x['sequence']
        l = 0
        r = len(seq)-1
        #从第一个为G的开始（‘首’C开头）
        while l<len(seq):
            if seq[l]!= 'C':
                l = l + 1
            else:
                break
        #‘尾’C结束
        while r>0:
            if seq[r]!= 'C':
                r = r - 1
            else:
                break
        
        x["sequence"] = seq[l:r+1]
        #x["sequence_size"] = len(x['sequence'])
        x["start"] = x['start']+l
        x["end"] = x['end'] - (len(seq)-1-r)
        arr = I_motif(x['sequence'])
        x["ISEAscore"] = sum(arr)/len(arr)

    #overlaps 去重叠

    data_json["data_view_overlaps"] = data_copy
    return data_json;

@app.route('/getISEAscore',methods = ['POST'])
def getISEAscore():
    sequence = None
    if request.method == 'POST':
        data = request.get_data()
        data = data.decode('utf-8')
        json_data = json.loads(data)

        try:
            sequence = json_data['sequence'].upper()
            if sequence.startswith('>'):
                resp = index_bash(request) #批量查询，处理翻页
            else:
                resp = index(request) #单次查询,不处理翻页
            return resp
        except:
            pass
        
        #except:
        #    resp = make_response({"msg":"服务器错误"})
        #    resp.status = '403'
        #    return resp
    else:
        resp = make_response({"msg":"请求方式不支持"})
        resp.statu = '403'
        return resp
#定义Restful接口，返回格式JSON
#@app.route('/getG4hunter',methods = ['POST'])
def index(request):
    sequence = None
    if request.method == 'POST':
        data = request.get_data()
        data = data.decode("utf-8")
        json_data = json.loads(data)

        #threshold_range = [1,1.2,1.5,1.75,2]

        try:
            #接受序列，不区分大小写，一律转换为大写
            sequence = json_data['sequence'].upper()
            #限制序列串长度不超过1w
            print (len(sequence))
            seq_size = len(sequence)
            if seq_size > 10000:
                resp = make_response({"msg":"序列长度过长"})
                resp.statu = '403'
                return resp
            else:
                windowed_value = int(json_data['windowed_value']) # default=25,范围：15-70
                threshold = float(json_data['threshold']) # value = [1.1,1.2,......,2.9,3.0] default = 1.5

        except:
            resp = make_response({"msg":"请求参数错误"})
            resp.statu = '403'
            return resp

        if not re.match(r'[ATGC]+$',sequence):
            resp = make_response({"msg":"序列内容包含除ATGC外的其它字符"})
            resp.statu = '403'
            return resp
        
        res = calculate_score(sequence,windowed_value,threshold)
        #res = json.dumps(res)
        
        temp = {}
        temp['gene_name'] = ''
        temp['data'] = res
        res = make_response(json.dumps(temp))
        res.status = '200'
        return res
    else:
        resp = make_response({"msg":"请求方式不支持"})
        resp.statu = '403'
        return resp

#批量计算I motif分数
#@app.route('/getG4hunter_bash',methods = ['POST'])
def index_bash(request):
    #批量查询，每段基因序列以>XXXX开头，然后ATGC字母再单独作为一行
    sequence_list = []
    pagsiz = 10 #每页记录数
    if request.method == 'POST':
        data = request.get_data()
        data = data.decode('utf-8')
        json_data = json.loads(data)
        result = []
        
        try:
            #以>划分序列为数组
            sequence_list = json_data['sequence'].split('>')
            windowed_value = int(json_data['windowed_value']) # default=25,范围：15-70
            threshold = float(json_data['threshold']) # value = [1.1,1.2,......,2.9,3.0] default = 1.5
            pagnum = int(json_data['pagnum'])
            if pagnum < 0 or pagnum >1000:
                resp = make_response({"msg":"翻页参数值不符合规范"})
                resp.statu = '403'
                return resp
            sequence_list = clean_nulldata(sequence_list)
            k_min = pagnum * pagsiz
            k_max = (pagnum+1) * pagsiz

            for k,i in enumerate(sequence_list):

                if k >= k_min and k < k_max:#筛选翻页记录
                    i = i.split('\n')
                    gene_name = i[0]
                    gene_seq = i[1]
                    seq_size = len(gene_seq)
                    if seq_size>10000:
                        resp = make_response({"msg":"序列长度过长"})
                        resp.statu = '403'
                        return resp
                    gene_data = {}
                    gene_data['gene_name'] = gene_name
                    
                    res = calculate_score(gene_seq,windowed_value,threshold)
                    gene_data['data'] = res
                    result.append(gene_data)
                else:
                    continue
            resp = make_response(json.dumps(result))
            resp.status = '200'
            return resp

        except:
            resp = make_response({"msg":"内部服务器错误"})
            resp.statu = '500'
            return resp

#文件上传，批量查询G4hunter分数
@app.route('/getG4hunter_uploadfile',methods = ['POST'])
def index_uploadfile():
    if request.method == 'POST':
        if upload_file(request):
            return "SUCCESS"
        else:
            return "ERROR"
    else:
        return "Message: The File format does not meet the requirements, please try again"


if __name__ == '__main__':
    
    #str = 'GGGTCTCGGTGGGGGTCACTGGGTTGCGGGCTTGG';
    #str = 'ACCCCCTGCATCTGCATGCCCCCTCCCACCCCCT'
    #str = 'TGAGGGTGGGTAGGGTGGGTAA'
    #str = 'GAGGACAAGGAGGTGCGAGGAAAGGGGTTGGGGGATGGTCCCACAGGCAGCCACACCTGAGGCGTGGGCGGCCGGTAGGAGCTGGGGGAGGGCGGGGAGAAGAGGGGTTTCTGTGTGTAGA'
    #print(calculate_score(str,25,-1.5))
    '''
    print('\n*** G4 ***')
    array = translate(str)
    print('length=',len(array))
    print('score=',sum(array)/len(array))
    print('\n*** I motif ***')
    array = I_motif(str)
    print('length=',len(array))
    print('score=',sum(array)/len(array))
    #d = [{'sequence': 'GGAGGTGCGAGGAAAGGGGTTGGGGGATGG', 'start': 9, 'end': 38, 'G4Hscore': 1.7666666666666666}, {'sequence': 'GGGCGGCCGGTAGGAGCTGGGGGAGGG', 'start': 66, 'end': 92, 'G4Hscore': 1.6666666666666667}, {'sequence': 'GGCCGGTAGGAGCTGGGGGAGGGCGGGGAGAAGAGGGGTTTCTGTG', 'start': 70, 'end': 115, 'G4Hscore': 1.5434782608695652}]
    #print(deduplication(d));
    #app.run(host='0.0.0.0',port=8080)
    '''
    CORS(app, supports_credentials=True)
    #app.run(host='47.100.189.241',port=5000,debug=True)
    app.run(host='0.0.0.0',port=5000,debug=True)
    