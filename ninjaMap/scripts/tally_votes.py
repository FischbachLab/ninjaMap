#!/usr/bin/env python3

import sys
import pandas as pd

def check_for_fraud(votes_file, discard_pile=1):
    fraud = True
    df = pd.read_csv(votes_file, index_col='Read_Name', usecols=['Read_Name', 'cSingular_Vote', 'cEscrow_Vote'])
    df['Total_Votes'] = round((df['cSingular_Vote'] + df['cEscrow_Vote']), 5)
    df = df.drop(['cSingular_Vote',  'cEscrow_Vote'], axis = 1)
    votes_df = df.groupby('Read_Name').sum()
    gt_row, gt_col = votes_df[(votes_df.Total_Votes > 1.001)].shape
    mid_row, mid_col = votes_df[(votes_df.Total_Votes > 0.001) & (votes_df.Total_Votes < 0.999)].shape
    if gt_row == 0 and mid_row == 0:
        return False

fraud_committed = check_for_fraud(votes_file)

if __name__ == '__main__':
    votes_file = sys.argv[1]
    fraud_committed = check_for_fraud(votes_file)
# if fraud_committed:
#     print('Evidence of voter fraud was found! :(')
# else:
#     print('No evidence of voter fraud detected! :)')
    fraud_committed.to_csv()