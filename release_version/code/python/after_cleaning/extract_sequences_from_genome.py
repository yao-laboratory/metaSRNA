import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

EXTRA_LEFT_GAP = 40
EXTRA_RIGHT_GAP = 40

def extract_sequences_from_genome(input_blast_result, input_genome, output_fasta):
    # input_blast_result = '../../../mytuber/results/SRR20723831/blast_result/blastn_mapping_SRR20723831_test.txt'
    df_blast_score = pd.read_csv(input_blast_result, sep=',')
    df_blast_score.columns = ["qseqid", "sacc", "sstart", "send", "evalue", "bitscore", "qcovhsp", "pident"]
    # Initialize an empty DataFrame to store the results
    result_df = pd.DataFrame(columns=df_blast_score.columns)

    # input_genome = '../../input_data/database/blast_mytuber_fna/test.fasta'
    seq = Seq("")
    #seq start index is 0
    #this following fasta is whole genome fasta file
    for record in SeqIO.parse(input_genome, "fasta"):
        # Create a Seq object from the sequence data
        seq = Seq(record.seq)
        # print(input_genome, len(seq))
    start_border = 1
    end_border = len(seq)
    print("part1 done")
    records = []
    with open(output_fasta, "w") as output_file:
        for index, row in df_blast_score.iterrows():
            seq_id = row['qseqid']
            start = int(row['sstart'])
            end = int(row['send'])
            
            if start < end:
                #number border protection
                df_start_pos = (start - EXTRA_LEFT_GAP) if (start - EXTRA_LEFT_GAP) >= start_border else start_border
                df_end_pos = (end + EXTRA_RIGHT_GAP) if (end + EXTRA_RIGHT_GAP) <= end_border else end_border
                # records.append({'seq_id': seq_id, 'seq': str(seq[df_start_pos: df_end_pos])})
                SeqIO.write(SeqRecord(Seq(str(seq[df_start_pos: df_end_pos + 1])), id=str(seq_id), description="index_"+str(index)), output_file, "fasta")
                # print(records[-1],df_start_pos,df_end_pos)
            else:
                df_start_pos = (end - EXTRA_LEFT_GAP) if (end - EXTRA_LEFT_GAP) >= start_border else start_border
                df_end_pos = (start + EXTRA_RIGHT_GAP) if (start + EXTRA_RIGHT_GAP) <= end_border else end_border
                SeqIO.write(SeqRecord(Seq(str(seq[df_start_pos: df_end_pos + 1].reverse_complement())), id=str(seq_id), description="index_"+str(index)), output_file, "fasta")
                #records.append({'seq_id': seq_id, 'seq': str(seq[df_start_pos: df_end_pos].reverse_complement())})
                # print(records[-1],df_end_pos,df_start_pos)
    # Write to a FASTA file
    #'../../../mytuber/results/SRR20723831/blast_result/test_predict_use.fasta'
    # with open(output_fasta, "w") as output_file:
    #     SeqIO.write(records, output_file, "fasta")
    # print(records)
    print("finish!")