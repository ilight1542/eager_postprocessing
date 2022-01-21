import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='merged output file')
parser.add_argument('-i', '--input', metavar="summary.txt files", required=True, nargs='+', help="multiple summary files to be merged")
parser.add_argument('-o', '--output', required=False, nargs=1, help="output file name (RunSummary.txt or TotalCount.txt)")
args = parser.parse_args()

def parse_summary_files(inputs):
    first=True
    for i in inputs:
        if first:
            df_out=pd.read_csv(i, sep='\t')
            first=False
        else:
            next=pd.read_csv(i, sep='\t')
            df_out=pd.merge(df_out,next,"outer")
    for cols in df_out.columns:
        if cols != "Node":
            df_out[cols]=df_out[cols].fillna(0)
            df_out[cols]=df_out[cols].astype(int)
    return df_out

def write_output(df_out, output_name):
    if ".txt" not in output_name:
        output_name=output_name+".txt"
    df_out.to_csv(output_name, index=False, sep='\t')

if __name__ == '__main__':
    write_output(parse_summary_files(args.input),args.output[0])