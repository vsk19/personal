# Translates Italian data to English and splits up data using standard delimiters


import pandas as pd
import numpy as np
import re
from googletrans import Translator


# returns true if the text is English, false if not
def is_english(text):
    lang = Translator().detect(text).lang
    return lang == 'en'


# returns english translation
def translate(to_translate):
    return Translator().translate(to_translate).text


def main():

    data = pd.read_excel('covid_data.xlsx', dtype = str)
    data = data[['Linked ID', 'Sample ID', 'co-morbidità']]
    data.rename(columns = {'co-morbidità': 'Comorbidity',
                              'Sample ID': 'Sample ID Raw'}, inplace = True)

    sample_key = pd.read_excel('sample_key.xlsx')
    sample_key = sample_key[['Database ID', 'External Sample ID']]
    sample_key.drop_duplicates(inplace = True)
    sample_key.rename(columns = {'External Sample ID': 'Sample ID'}, inplace = True)
    sample_key['Sample ID'] = sample_key['Sample ID'].str.strip()

    # clean up data sample ids in separate column
    data['Matched Sample ID'] = None
    data['Database ID'] = None
    data['ID Notes'] = None

    for i, row in data.iterrows():
        ids = row['Sample ID Raw']

        # split string into list of ids
        if '=' in ids:
            ids = ids.split("=")
        elif ',' in ids:
            ids = ids.split(',')
        else:
            ids = re.split(' {2,}', ids)

        # Check every ID. Stop when we've found an ID
        for id in ids:

            # don't care about blanks
            if id.isspace() or id == '':
                continue
            test_id = id.upper()
            test_id = test_id.strip()

            # Look for exact match
            matches = sample_key[sample_key['Sample ID'] == test_id]
            if matches.shape[0] > 0:
                database_ids = matches['Database ID'].unique()

                # if multiple exact matches, we have a problem!
                if len(database_ids) > 1:
                    print('EXACT TESTING: More than one Database ID match!')
                    print(database_ids)
                    print(test_id)
                    print('*****************************')
                else:
                    data.at[i, 'Matched Sample ID'] = id.strip()
                    data.at[i, 'Database ID'] = database_ids[0]
                    data.at[i, 'ID Notes'] = 'Exact match'
                    break

            # Look for timestamp at the end like " T0" and replace space with "-"
            if len(re.findall(r"\sT\d+$", test_id)) == 1:
                print('This one has a timestamp: ' + test_id)
                timestamp_id = test_id.replace(' ', '-')
                matches = sample_key[sample_key['Sample ID'].str.contains(timestamp_id, regex=False)].reset_index(drop=True)
                if matches.shape[0] > 0:
                    database_ids = matches['Database ID'].unique()

                    # if multiple exact matches, we have a problem!
                    if len(database_ids) > 1:
                        print('Timestamp already there: More than one Database ID match!')
                        print(database_ids)
                        print(test_id)
                        print(timestamp_id)
                        print('------------------------------')
                    elif len(database_ids) == 1:
                        data.at[i, 'Matched Sample ID'] = id.strip()
                        data.at[i, 'Database ID'] = database_ids[0]
                        data.at[i, 'ID Notes'] = 'Matched time point (' + str(matches.at[0, 'Sample ID']) + ')!!!!'
                        break
            # Timestamp isn't already in the ID. Add it and see if match
            else:
                # See if listed sample IDs contain this test_id with time point info (Ex: ID-T0)
                time_id = test_id + '-T'
                matches = sample_key[sample_key['Sample ID'].str.contains(time_id, regex = False)].reset_index(drop = True)
                if matches.shape[0] > 0:
                    database_ids = matches['Database ID'].unique()

                    # if multiple exact matches, we have a problem!
                    if len(database_ids) > 1:
                        print('Added timestamp: More than one Database ID match!')
                        print(database_ids)
                        print(test_id)
                        print(time_id)
                        print('------------------------------')
                    elif len(database_ids) == 1:
                        data.at[i, 'Matched Sample ID'] = id.strip()
                        data.at[i, 'Database ID'] = database_ids[0]
                        data.at[i, 'ID Notes'] = 'Matched time point (' + str(matches.at[0, 'Sample ID']) + ')'
                        break

    data = data[['Database ID', 'Matched Sample ID', 'ID Notes', 'Linked ID', 'Sample ID Raw', 'Comorbidity']]

    processed = pd.DataFrame(columns = ['Database ID', 'Matched Sample ID', 'ID Notes', 'Linked ID', 'Sample ID Raw', 'Original Comorbidity', 'Split', 'Final Comorbidity'])

    # Split up the comorbidities line by line
    for i, row in data.iterrows():
        comorbs_raw = row['Comorbidity']
        comorbs = str(comorbs_raw)

        # if text has periods, split by that
        if comorbs.count('.') >= 3:
            comorbs = comorbs.split('.')
        else:
            comorbs = [comorbs]

        for comorb in comorbs:
            translated_comorb = ''
            if ~is_english(comorb):
                translated_comorb = translate(comorb)
            else:
                translated_comorb = comorb
            new_row = pd.Series([row['Database ID'], row['Matched Sample ID'], row['ID Notes'], row['Linked ID'], row['Sample ID Raw'], comorbs_raw, comorb, translated_comorb],
                                index = ['Database ID', 'Matched Sample ID', 'ID Notes', 'Linked ID', 'Sample ID Raw', 'Original Comorbidity', 'Split', 'Final Comorbidity'])
            processed = processed.append(new_row, ignore_index=True)

    processed.to_csv('comorb_processed.csv', index = False)


if __name__ == '__main__':
    main()
