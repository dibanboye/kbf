
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
        if 'Input finished,' in lin:
            kmer.append(int(lin.split()[6]))
            # print(lin.split())
    res = {
        'methods': methods,
        'populate': populate,
        'accuracy': accuracy,
        'query': query,
        'num_of_query': num_of_query,
        'kmer': kmer,
        'kmer_reading_time': getting_time
    }
    print (res)
    return res