'''
1 清洗数据
2 读取用户上传文件
'''
import re
import json,os
from flask import request
from werkzeug.utils import secure_filename
allow_format = ['fasta','txt']

def clean_nulldata(arraylist):
    newlist = []
    for i in arraylist:
        if i != '': #去空
            newlist.append(i)
    for i in newlist:
        print('基因：')
        print(i)
        i = i.split('\n')[1] #以换行切割，得到想要的基因序列    
    return newlist

def upload_file(request):

    if (request.method == 'POST'):
        print('文件上传')
        #获取上传文件对象
        file = request.files['file']
        #判断文件后缀是否为fasta/txt格式，否则报错
        filename = secure_filename(file.filename)
        print(filename)
        if not allowed_file(filename):#文件后缀格式错误
            return False
        else:
            #读取文件内容
            fo = file.read()
            #切割，分析每段基因，返回数组
            file_data = str(fo)
            print('文件内容=',file_data)

            
            return True
    else:
        return False


def allowed_file(filename):
    # 获取文件扩展名，以'.'为右分割然后取第二个值
    return '.' in filename and filename.rsplit('.',1)[1].lower() in allow_format

if (__name__ == '__main__'):
    f = allowed_file('aaa.fasta')
    print(f)
