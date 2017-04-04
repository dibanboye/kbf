from flask import Flask
from flask import Markup
from flask import Flask, request, flash, redirect, url_for
from flask import render_template

from parse import parse_log
import sys, os
app = Flask(__name__)
 
def parse(dataset, K, query):
    d = ''
    k = ''
    q = ''

    if dataset == 'test':
        d = 'test.fasta'
    elif dataset == 'yeast':
        d = 'yeast.fasta'
    elif dataset == 'ecoli':
        d = 'sra_data.fasta'

    if query == '0.1M':
        q = '500000'
    elif query == '1M':
        q = '1000000'
    elif query == '10M':
        q = '10000000'
    
    if K == 'all':
        k = [5, 10, 20, 25]
    else:
        k = int(K)
    
    return d, k, q

@app.route("/", methods=['GET', 'POST'])
@app.route("/main", methods=['GET', 'POST'])
def index():
    # if request.method == 'POST' or request.method == 'GET':
    dataset = request.form.get("dataset")
    K = request.form.get("K")
    query = request.form.get("query")

    if dataset == None or K == None or query == None:
        return render_template('index.html', start=0)
    else:    
        print ('###', request.form['submit'])
        dataset, K, query = parse(dataset, K, query)

        print (dataset)
        print (K)
        print (query)
        path = '/home/bioinf/project_zizhao_test/KBF/cpp-src/'
        demo_name = 'log/log_demo_{}_{}_{}.txt'.format(dataset,K,query)
        

        if not os.path.isfile(demo_name) or request.form['submit'] == 'go':
            cmd = "{}main data/{} {} {} {} | tee {}".format(path, dataset, K, 'test_demo', query, demo_name)
            print (cmd)
            os.system(cmd)
        else:
            print ('{} file exist!'.format(demo_name))

        if os.path.isfile(demo_name):
            res = parse_log(demo_name)
            return render_template('index.html', 
                        start=1,
                        methods=res['methods'],
                        populate=res['populate'],
                        accuracy=res['accuracy'],
                        query=res['query'],
                        )    
        else:
            return render_template('index.html', start=0)


if __name__ == "__main__":
    app.run(host='0.0.0.0', port=9999, debug=True)