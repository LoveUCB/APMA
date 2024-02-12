def other_feature(path):
    '''
    This function is used to bind the feature provided by users to the final paras file
    :return: lists
    '''

    with open(path, 'r') as file:
        current_data = []
        current_cata = []
        for line in file:
            columns = line.strip().split('\t')
            current_data.append(columns[1:])
