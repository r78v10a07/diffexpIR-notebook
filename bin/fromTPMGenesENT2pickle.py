import os
import pandas
import argparse
import numpy as np

from six.moves import cPickle as pickle


def extract_data(ctrl_dir, ctl_samples, treat_dir, treat_samples, log_file):
    samples = {}
    for file in ctl_samples:
        log_file.write("Parsing file %s\n" % os.path.join(ctrl_dir, file))
        f_log.flush()
        data = pandas.read_csv(os.path.join(ctrl_dir, file), sep='\t')
        for chr in data.Chr.unique():
            if '_' not in chr:
                if chr not in samples:
                    samples[chr] = {}
                for g in data.loc[data['Chr'] == chr, 'Gene_Id'].unique():
                    if g not in samples[chr]:
                        samples[chr][g] = {}
        for index, row in data.iterrows():
            chr = row['Chr']
            if '_' not in chr:
                g = row['Gene_Id']
                c = int(row['Type_Number'])
                start = int(row['start']) + 1
                end = int(row['end']) + 1
                if c not in samples[chr][g]:
                    samples[chr][g][c] = {'type': row['Type'], 'start': start, 'end': end,
                                          'TPM_list_ctrl': [], 'Reads_list_ctrl': [], 'TPM_ctrl': 0,
                                          'Reads_ctrl': 0,
                                          'TPM_list_treat': [], 'Reads_list_treat': [], 'TPM_treat': 0,
                                          'Reads_treat': 0}
                tpm = float(row['TPM'])
                if tpm < 10E-5:
                    tpm = 10E-5
                samples[chr][g][c]['TPM_list_ctrl'].append(tpm)
                samples[chr][g][c]['Reads_list_ctrl'].append(int(row['Count_Reads']))
                samples[chr][g][c]['TPM_ctrl'] = np.mean(samples[chr][g][c]['TPM_list_ctrl'])
                samples[chr][g][c]['Reads_ctrl'] = np.mean(samples[chr][g][c]['Reads_list_ctrl'])
    for file in treat_samples:
        log_file.write("Parsing file %s\n" % os.path.join(treat_dir, file))
        f_log.flush()
        data = pandas.read_csv(os.path.join(treat_dir, file), sep='\t')
        for chr in data.Chr.unique():
            if '_' not in chr:
                if chr not in samples:
                    samples[chr] = {}
                for g in data.loc[data['Chr'] == chr, 'Gene_Id'].unique():
                    if g not in samples[chr]:
                        samples[chr][g] = {}
        for index, row in data.iterrows():
            chr = row['Chr']
            if '_' not in chr:
                g = row['Gene_Id']
                c = int(row['Type_Number'])
                start = int(row['start']) + 1
                end = int(row['end']) + 1
                if c not in samples[chr][g]:
                    samples[chr][g][c] = {'type': row['Type'], 'start': start, 'end': end,
                                          'TPM_list_ctrl': [], 'Reads_list_ctrl': [], 'TPM_ctrl': 0,
                                          'Reads_ctrl': 0,
                                          'TPM_list_treat': [], 'Reads_list_treat': [], 'TPM_treat': 0,
                                          'Reads_treat': 0}
                tpm = float(row['TPM'])
                if tpm < 10E-5:
                    tpm = 10E-5
                samples[chr][g][c]['TPM_list_treat'].append(tpm)
                samples[chr][g][c]['Reads_list_treat'].append(int(row['Count_Reads']))
                samples[chr][g][c]['TPM_treat'] = np.mean(samples[chr][g][c]['TPM_list_treat'])
                samples[chr][g][c]['Reads_treat'] = np.mean(samples[chr][g][c]['Reads_list_treat'])
    return samples


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Clean SNPdata table')
    parser.add_argument('-s', help='cll_samples.tsv', required=True)
    parser.add_argument('-l', help='Log file', required=True)
    parser.add_argument('-p', help='Pickle file', required=True)
    args = parser.parse_args()
    sample_data = args.s
    log_file = args.l
    pickle_file = args.p

    with open(log_file, "w") as f_log:
        f_log.write('Reading CLL sample files\n')
        cll_samples = pandas.read_csv(sample_data, sep='\t')

        f_log.write('Creating control and CLL sample list\n')
        samples = {'ctrl': [], 'cll': [], 'mut': [], 'unm': []}

        files = [name for root, dirs, files in os.walk('./') for name in files if name.endswith("_genes.ent")]
        for f in files:
            n = f.split('_')[0]
            row = cll_samples.loc[cll_samples['ECGA_ID'] == n]
            if row.values[0][1] == 'NBC':
                samples['ctrl'].append(f)
            else:
                samples['cll'].append(f)
                if row.values[0][1] == 'MUT':
                    samples['mut'].append(f)
                elif row.values[0][1] == 'UNM':
                    samples['unm'].append(f)
        f_log.write('There are ' + str(len(samples['ctrl'])) + ' control samples\n')
        f_log.write('There are ' + str(len(samples['cll'])) + ' CLL samples\n')
        f_log.flush()
        if not os.path.exists(pickle_file):
            try:
                f_log.write('Extracting sample data\n')
                f_log.flush()
                samples_data = extract_data('./', samples['ctrl'], './', samples['cll'], f_log)
                samples_data_1 = extract_data('./', samples['ctrl'], './', samples['mut'], f_log)
                samples_data_2 = extract_data('./', samples['ctrl'], './', samples['unm'], f_log)
                samples_data_3 = extract_data('./', samples['mut'], './', samples['unm'], f_log)
                f_log.write('Pickling data\n')
                f_log.flush()
                with open(pickle_file, 'wb') as f:
                    pickle.dump(samples_data, f, pickle.HIGHEST_PROTOCOL)
                    pickle.dump(samples_data_1, f, pickle.HIGHEST_PROTOCOL)
                    pickle.dump(samples_data_2, f, pickle.HIGHEST_PROTOCOL)
                    pickle.dump(samples_data_3, f, pickle.HIGHEST_PROTOCOL)
            except Exception as e:
                print('Unable to save data to', pickle_file, ':', e)
f_log.write('Done\n')
