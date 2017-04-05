import copy


def parse_log(file_name):
    f = open(file_name)

    scores = dict()
    data = f.readlines()
    
    methods = []
    populate = []
    accuracy = []
    query = []
    num_of_query = []
    kmer = []
    BF_size = 0
    for c, lin in enumerate(data):
        if lin[0:3] == '###':
            method = lin.split(' ')
            method = ' '.join(method[1:-1])
            method = '_'.join(method.split(' '))
            methods.append(method)
        if lin[0:2] == '@@':
            if 'populate' in lin:
                temp  = lin.split(' ')[-1].rstrip()
                populate.append(temp)
            elif 'query time' in lin:
                temp  = lin.split(' ')[-1].rstrip()
                query.append(temp)
            elif 'Accuracy' in lin:
                temp  = lin.split(' ')[-1].rstrip()
                accuracy.append(temp)
            elif 'queries' in lin:
                temp  = lin.split(' ')[-1].rstrip()
                num_of_query.append(temp)
            elif 'Getting' in lin:
                getting_time = lin.split(' ')[-2].rstrip()
                kmer.append('%.3f'%(float(lin.split(' ')[1])))
        if 'Input finished,' in lin:
            kmer.append(int(lin.split()[6]))
            # print(lin.split())
           
        if 'BF_size' in lin:
        
            temp  = lin.split(' ')[-1].rstrip()
            BF_size = '%.3f'%(float(temp)/(1024*1024))

    res = {
        'methods': methods,
        'populate': populate,
        'accuracy': accuracy,
        'query': query,
        # 'num_of_query': num_of_query,
        'kmer': kmer,
        'kmer_reading_time': getting_time,
        'BF_size': BF_size
    }
    print (res)
    return res

# parse multi files
def parse_logs(file_names):
    
    res = []
    for i, file_name in enumerate(file_names):
        print ('parse ' + file_name)
        res.append(parse_log(file_name))

    num_methods = len(res[0]['methods'])
    print ('-'*10)
    print (res)
    print ('number of methods', num_methods)
    out = []
    for i in range(num_methods):
        local_res = {
            'methods': [],
            'populate': [],
            'accuracy': [],
            'query': [],
            # 'num_of_query': [],
            'kmer': [],
            'kmer_reading_time': [],
            'BF_size': []
        }
        for k in local_res.keys():
            print (k)
            if type(res[0][k]) != list:
                local_res[k] = [res[c][k] for c in range(len(res))]
            else:
                tk = len(res[c][k])
                local_res[k] += [res[c][k][i%tk] for c in range(len(res))]

       
        out.append(copy.deepcopy(local_res))
    return out

if __name__ == '__main__':
    file_names = ['temp/log_demo_yeast.fasta_20_500000.txt',  'temp/log_demo_yeast.fasta_25_500000.txt']

    out = parse_logs(file_names)
