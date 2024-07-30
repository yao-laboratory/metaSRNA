import pandas as pd
import os
from os import listdir, path, makedirs



def main(sample_id, start, end, filename):
    folder_path = path.abspath(path.join(path.dirname(__file__), f"../../../mytuber/results/{sample_id}/blast_result/"))
    f = path.join(folder_path, f'{filename}.txt')
    f_new = path.join(folder_path, f'{filename}_after_filter.csv')
    # Read the .txt file
    df = pd.read_csv(f, header=None)
    df.columns = ['qseqid','sacc','sstart','send','evalue','bitscore','qcovhsp','pident']
    # print(df)
    df.to_csv(f_new, sep='\t', index=False)
    # Read the CSV file into a DataFrame

    df_new = pd.read_csv(f_new, sep='\t')
    # print(df_new['send'] )
    # Calculate the absolute difference between the first two columns
    if df_new['send'].isna().any() or df_new['sstart'].isna().any():
        print("There are NaN values in the columns.")
    else:
        df_new['difference'] = abs(df_new['send'] - df_new['sstart']) + 1
        #print(df_new['difference'])
        # Filter the DataFrame to keep only the rows in difference conditions
        df_filtered = df_new[(df_new['difference'] >=int(start)) & (df_new['difference'] <=int(end))]
        #print(df_filtered)
        df_filtered.to_csv(f_new, index=False)



if __name__ == "__main__":
    sample_id = os.getenv('PARAM1', 'default_value1')
    start = os.getenv('PARAM2', 'default_value2')
    end = os.getenv('PARAM3', 'default_value3')
    filename = os.getenv('PARAM4', 'default_value4')
    print(sample_id, start, end, filename)
    main(sample_id, start, end, filename)
    print("finished filtering.")