import sys
import os
import re
import datetime
import glob
import pandas as pd
import numpy as np
import argparse
from docx2python import docx2python
import argparse
from argparse import RawTextHelpFormatter
from ncbr_huse import send_update, err_out, pause_for_input
from ncbr_bsi import read_conf, send_curl, get_bsi_session, get_bsi_name, bsi_query


# database_df passes in all variants that are shared across families.
# Returns subject level variant table. Each row is a Patient ID and a variant that the person has
def build_subject_variant_tab(database_df, id_dict, fam_dict, sample_id_dict):
    subj_variants = pd.DataFrame(
        columns=['Sample_ID', 'Patient_Name', 'Patient_ID', 'database_ID', 'Patient_Fam_ID', 'database_Fam_ID',
                 'hgvsc', 'hgvsp', 'Gene', 'Chrom', 'Pos', 'Ref', 'Alt', 'Genotype',
                 'in_HGMD', 'Zygosity', 'Tags', 'Notes', 'Effect', 'sift', 'polyphen'])

    # loop thru all rows of database which contains variants shared across families
    for ind, row in database_df.iterrows():
        # go thru all samples in database row
        for i in range(0, 10):
            colname = 'GT_genotype_' + str(i)
            samplename = 'sample_id_' + str(i)
            genotype = str(row[colname])

            database_id = str(row[samplename])
            database_id = database_id[:database_id.find('_')]
            phen_id = id_dict.get(database_id)
            sample_id = sample_id_dict.get(database_id)

            # Don't allow blanks
            if '/' in genotype and genotype != './.':

                # don't allow any hom ref subjects in the variants output
                if '(hom_ref)' not in genotype:
                    database_fam = str(row['family'])
                    database_fam = database_fam[:database_fam.find('_')]
                    phen_fam = fam_dict.get(database_fam)

                    hom_alt = row['alt'] + '/' + row['alt']
                    htz = row['ref'] + '/' + row['alt']

                    zygosity = 'Unknown'

                    if genotype == hom_alt:
                        zygosity = 'Homozygous Alt'
                    elif genotype == htz:
                        zygosity = 'Heterozygous'
                    elif genotype[:genotype.find('/')] != row['ref'] and genotype[genotype.find('/') + 1:] != row['ref']:
                        zygosity = 'Compound Het'

                    # populate individual level row from database data into subj_variants
                    subj_variants = subj_variants.append({
                        'Sample_ID': sample_id,
                        'Patient_ID': phen_id,
                        'database_ID': database_id,
                        'Patient_Fam_ID': phen_fam,
                        'database_Fam_ID': database_fam,
                        'hgvsc': row['hgvsc'],
                        'hgvsp': row['hgvsp'],
                        'Gene': row['gene'],
                        'Chrom': row['chrom'],
                        'Pos': row['pos'],
                        'Ref': row['ref'],
                        'Alt': row['alt'],
                        'Genotype': genotype,
                        'in_HGMD': row['in_hgmd'],
                        'Zygosity': zygosity,
                        'Tags': row['tags'],
                        'Notes': row['notes'],
                        'Effect': row['effect'],
                        'sift': row['sift'],
                        'polyphen': row['polyphen']}, ignore_index=True, sort = False)

    subj_variants.sort_values(['Patient_ID'], inplace=True)

    return subj_variants


# Returns all rows from database that are tagged with Inconclusive Negative, but are missing actual Negative Patient reports
def check_inconclusive_neg_variants(variants, database, sample_id_df):

    # makes dicts out of our Sample_ID dictionary
    id_dict = dict(zip(sample_id_df.database_ID, sample_id_df.Patient_ID))
    fam_dict = dict(zip(sample_id_df.database_Family_ID, sample_id_df.Patient_Family_ID))
    sample_id_dict = dict(zip(sample_id_df.database_ID, sample_id_df.Sample_ID))

    # keep only Inconclusive negative report tagged variants
    database = database[database['tags'].str.contains('Inconclusive negative report', na=False)]

    # Get subject level table for easier comparison
    database_individual = build_subject_variant_tab(database, id_dict, fam_dict, sample_id_dict)

    no_reports = pd.DataFrame(columns=['Sample_ID', 'Patient_ID', 'database_ID', 'Patient_Fam_ID', 'database_Fam_ID',
                                       'hgvsc', 'hgvsp', 'Gene', 'Chrom', 'Pos', 'Ref', 'Alt', 'Genotype',
                                       'Zygosity', 'Tags', 'Notes', 'Effect', 'sift', 'polyphen', 'in_HGMD'])

    # Checks each database row against list of Negative reports
    for ind, row in database_individual.iterrows():
        sample_id = row['Sample_ID']

        # cut down variants to those that match the Sample_ID
        have_neg_reports = variants[variants['Sample_ID'] == sample_id]

        # if no reports to be found, save the database row
        if have_neg_reports.empty:
            no_reports = no_reports.append(row, sort=False)

    return no_reports


# Returns all rows from database that are tagged with Patient Report, but are missing actual Positive Patient reports
def check_patient_report_variants(variants, database, sample_id_df):

    # makes dicts out of our sample_id dictionary
    id_dict = dict(zip(sample_id_df.database_ID, sample_id_df.Patient_ID))
    fam_dict = dict(zip(sample_id_df.database_Family_ID, sample_id_df.Patient_Family_ID))
    sample_id_dict = dict(zip(sample_id_df.database_ID, sample_id_df.Sample_ID))

    # Limit to variants with Patient Report Tag
    database = database[database['tags'].str.contains('Patient Report', na=False)]

    # Get subject level table for easier comparison
    database_individual = build_subject_variant_tab(database, id_dict, fam_dict, sample_id_dict)

    no_reports = pd.DataFrame(columns=['Sample_ID', 'Patient_ID', 'database_ID', 'Patient_Fam_ID', 'database_Fam_ID',
                                       'hgvsc', 'hgvsp', 'Gene', 'Chrom', 'Pos', 'Ref', 'Alt', 'Genotype',
                                       'Zygosity', 'Tags', 'Notes', 'Effect', 'sift', 'polyphen', 'in_HGMD'])

    # Checks each database row against list of Patient Reports
    for ind, row in database_individual.iterrows():
        sample_id = row['Sample_ID']
        dna_change = row['hgvsc']

        # cut down variants to those that match the Sample_ID and the DNA Change
        have_patient_reports = variants[variants['Sample_ID'] == sample_id]
        have_patient_reports = have_patient_reports[have_patient_reports['hgvsc'] == dna_change]

        # if no reports to be found, save the database row
        if have_patient_reports.empty:
            no_reports = no_reports.append(row, sort=False)

    return no_reports


# Copies info from database into the variants table:
#   variants = df with all scraped variants
#   variant_ind = row in variants table to look for ID / assign variant data
#   database_row = database row that has matched dna/protein change
def copy_info(variants, variant_ind, database_row):
    matched_id = False

    # Check all IDs in database row for a match with the variant row ID
    for col in range(0, 10):
        sample_name = 'sample_id_' + str(col)
        database_id = str(database_row[sample_name])

        # Found a match! Adds database variant info to df
        if variants.at[variant_ind, 'database_ID'] in database_id:
            matched_id = True
            variants.at[variant_ind, 'Chrom'] = database_row['chrom']
            variants.at[variant_ind, 'Pos'] = database_row['pos']
            variants.at[variant_ind, 'Ref'] = database_row['ref']
            variants.at[variant_ind, 'Alt'] = database_row['alt']
            variants.at[variant_ind, 'hgvsc'] = database_row['hgvsc']
            variants.at[variant_ind, 'hgvsp'] = database_row['hgvsp']
            variants.at[variant_ind, 'in_HGMD'] = database_row['in_hgmd']
            variants.at[variant_ind, 'Tags'] = database_row['tags']
            variants.at[variant_ind, 'Notes'] = database_row['notes']
            break

    return variants, matched_id


# Adds variants that match DNA Change from database
def add_database_variants(variants, database):
    # Add database fields
    variants['Chrom'] = None
    variants['Pos'] = None
    variants['Ref'] = None
    variants['Alt'] = None
    variants['hgvsc'] = None
    variants['hgvsp'] = None
    variants['in_HGMD'] = None
    variants['Tags'] = None
    variants['Notes'] = None
    variants['Comparison Notes'] = None

    # Loop through each row in scraped variants
    for i in range(0, variants.shape[0]):

        # stores if family id exists in the database download (for diagnostic purposes)
        family_exists = database['family'].str.contains(str(variants.at[i, 'database_Family_ID']), na=False).any()

        # Removes any spaces from rows that look like c.2164 G>A for example
        dna_change = str(variants.at[i, 'DNA Change'])
        if dna_change.startswith('c.'):
            dna_change = dna_change.replace(' ', '')
            variants.at[i, 'DNA Change'] = dna_change

        # Removes spaces from protein change field
        protein_change = str(variants.at[i, 'Protein Change'])
        if protein_change.startswith('p.'):
            protein_change = protein_change.replace(' ', '')
            variants.at[i, 'Protein Change'] = protein_change

        # Checks for Gene field match with database
        database_rows = database[database['gene'] == variants.at[i, 'Gene']]
        if database_rows.empty:
            variants.at[i, 'Comparison Notes'] = 'No Gene match'
            continue

        matched_id = False
        matched_dna_change = False
        matched_notes = False
        matched_protein = False
        fam_id = variants.at[i, 'database_Family_ID']

        # At this point, database_rows contains all database rows with a matching gene
        # Loop through each gene match database_row and check dna_change, notes, protein_change
        for ind, row in database_rows.iterrows():
            # checks IDs if the DNA Change is actually there
            if dna_change in row['hgvsc']:
                matched_dna_change = True
                variants, matched_id = copy_info(variants, i, row)

            else:
                # Try to see if the variant nomenclature shows up in the notes
                if 'nomenclature' in str(row['notes']) or 'Nomenclature' in str(row['notes']):
                    if dna_change in str(row['notes']):
                        matched_notes = True
                        variants, matched_id = copy_info(variants, i, row)

                # if not in the notes, try looking at protein change
                else:
                    # if there's a protein change match, check if the DNA match field has same number of base deletion
                    if protein_change in str(row['hgvsp']):

                        num_bases = dna_change.split('del')

                        # get number of bases from DNA Change field and compare to database
                        if len(num_bases) > 1 and num_bases[1].isdigit():
                            num_bases = int(num_bases[1])
                            num_bases_database = row['hgvsc'].split('del')

                            # compares dna change base deletion length to database base deletion length
                            if len(num_bases_database) > 1:
                                num_bases_database = len(num_bases_database[1])
                                if num_bases == num_bases_database:
                                    matched_protein = True
                                    variants, matched_id = copy_info(variants, i, row)

            if matched_id:
                if matched_notes:
                    variants.at[i, 'Comparison Notes'] = 'Variant nomenclature was found in notes. Not an exact match.'
                elif matched_protein:
                    variants.at[
                        i, 'Comparison Notes'] = 'Protein Change match found in database. DNA Change was not an exact match.'
                break

        if not matched_id:
            if matched_dna_change:
                variants.at[i, 'Comparison Notes'] = 'Other subjects in database have this variant, but not this subject.'
            elif matched_notes:
                variants.at[
                    i, 'Comparison Notes'] = 'Variant nomenclature was found in notes in other subjects, but not this subject.'
            elif matched_protein:
                variants.at[
                    i, 'Comparison Notes'] = 'Protein Change match found in database for other subjects, but not this subject.'
            else:
                variants.at[i, 'Comparison Notes'] = 'No DNA Change match'

        if not family_exists:
            variants.at[i, 'Comparison Notes'] = 'database Family ID not found in database download. ' + str(
                variants.at[i, 'Comparison Notes'])

    return variants


# Returns Sample_ID, database, Patient ID dictionary
def make_id_dict():
    # Set up the variables, bsi info
    cnf = os.path.expanduser('~/.my.cnf.bsi')
    url_session = 'https://rest.bsisystems.com/api/rest/EBMS/common/logon'
    url_reports = 'https://rest.bsisystems.com/api/rest/EBMS/reports/list'
    curl_get = "curl -s -X GET --header 'Accept: application/json' --header 'BSI-SESSION-ID: "

    print('Making Sample_ID to Patient dictionary from BSI...')
    # Establish a BSI Connection with user's credentials
    user, pw = read_conf(cnf)
    session = get_bsi_session(url_session, user, pw)

    # Build sample_id to Patient dictionary
    fields = ['Sample_ID', 'Patient Name', 'Patient ID', 'Patient Family ID', 'Batch Received', 'Patient Order Status', 'Exome ID']
    sample_id_dict = bsi_query(curl_get, url_reports, session, fields, ['BATCH*', 'Gene*'], 'Batch Received', islike=True)
    sample_id_dict.drop_duplicates(inplace=True)
    sample_id_dict = sample_id_dict[~sample_id_dict['Patient_Order_Status'].str.contains('Cancel', na = False, case = False)]
    sample_id_dict = sample_id_dict[~sample_id_dict['Patient_Order_Status'].str.contains('Auto Complete', na = False)]
    sample_id_dict = sample_id_dict[~sample_id_dict['Exome ID'].isna()]
    sample_id_dict = sample_id_dict[['Sample_ID', 'Patient_Name', 'Patient_ID', 'Patient_Family_ID', 'Batch_Received']]
    sample_id_dict.columns = ['Sample_ID', 'Patient_Name', 'Patient_ID', 'Patient_Family_ID', 'Batch_Received']
    sample_id_dict.drop_duplicates(inplace=True)

    # Import the database to Patient dictionary
    database_dict = pd.read_csv('samplekey.txt', dtype=str, sep='\t', names=['Patient_ID', 'database_ID'])
    database_fam_dict = pd.read_csv('samplekey_fam.txt', dtype=str, sep='\t', names=['Patient_Family_ID', 'database_Family_ID'])

    sample_id_dict = sample_id_dict.merge(database_dict, on='Patient_ID', how='left')
    sample_id_dict = sample_id_dict.merge(database_fam_dict, on='Patient_Family_ID', how='left')

    return sample_id_dict


def get_secondary_findings(variants):

    secondary = variants[variants['Tags'].str.contains('Secondary finding', na = False)]
    secondary = secondary[secondary['Tags'].str.contains('Patient Report', na = False)]
    secondary['to_remove'] = ''

    # Loop through all secondary findings reports
    # and check if the words "Secondary findings: Not detected." are in the report
    # if so, remove from the secondary findings list, and put a note in Comparison Notes that
    # the report has been incorrectly tagged in database for this person
    for i, row in secondary.iterrows():

        received = secondary.at[i, 'Batch_Received']
        sample_id = secondary.at[i, 'sample_id']
        sample_id = sample_id.replace("-", "")
        batch = int(re.search(r'\d+', received).group())
        batch_name = 'Batch ' + str(batch)

        report_types = pd.Series(glob.glob('/Path/to/Patient/Reports/' + batch_name + '/*'))
        report_dir = report_types[report_types.str.contains('Positive', case=False)]
        if report_dir.empty:
            continue
        report_dir = report_dir.reset_index(drop=True)[0]

        report_path = report_dir + "/" + sample_id + "*.docx"
        print(report_path)

        fnames = glob.glob(report_path)
        fnames = [x for x in fnames if "~$" not in x]
        print(fnames)

        doc = docx2python(fnames[0])
        text = doc.text.lower()

        if "secondary findings: not detected" in text:
            if variants.at[i, 'Comparison Notes'] is None:
                variants.at[i, 'Comparison Notes'] = "Tagged in database with secondary finding, but not detected in actual patient report."
            else:
                variants.at[i, 'Comparison Notes'] = variants.at[i, 'Comparison Notes'] + ", Tagged in database with secondary finding, but not detected in actual patient report."
            secondary.at[i, 'to_remove'] = 'Remove'

    secondary = secondary[~secondary['to_remove'].str.match('Remove', na = False)]
    secondary.drop(columns = ['to_remove'], inplace = True)

    return secondary, variants


def main():

    #### Usage statement ####
    parsestr = 'Scrapes all positive reports for variant information.\n\
            Includes Negative reports if flagged. Outputs a CSV.\n\n\
            Usage:\n\
                scrape_report_docs.py -o outdir -s database_files_dir \n\n'

    parser = argparse.ArgumentParser(description=parsestr, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-o', '--outdir', required=True, type=str, help='Output file directory')
    parser.add_argument('-n', '--negative_include', required=False, action='store_true', default=False,
                        help='Also include negative reports in output.')
    parser.add_argument('-s', '--database_dir', required=True, type=str, help='Directory where database downloads are stored')

    args = parser.parse_args()
    outdir = args.outdir
    include_neg = args.negative_include
    database_dir = args.database_dir

    if not database_dir.endswith("/"):
        database_dir = database_dir + "/"
    if not outdir.endswith("/"):
        outdir = outdir + "/"

    ####### Make a combined database file from all the files in database_dir #######
    print('Combining database files')
    database_files = glob.glob(database_dir + '*.csv')

    database_combined = pd.DataFrame()
    for file in database_files:
        curr_table = pd.read_csv(file, dtype = str)
        database_combined = database_combined.append(curr_table)

    database_combined.drop_duplicates(inplace = True)
    database_combined.to_csv('database_combined.csv', index = False)
    flip_reports = 2
    if include_neg:
        flip_reports = 3
        neg_patients = pd.DataFrame(columns=['sample_id', 'Batch'])
        print('Including negative reports.')

    # sample_id_dict = make_id_dict()
    # sample_id_dict.to_csv(outdir + 'dict.csv', index=False)
    sample_id_dict = pd.read_csv(outdir + 'dict.csv', dtype=str)

    variants = pd.DataFrame()
    cma = pd.DataFrame()

    # Make sure we're scraping all available batches in the FINAL directory
    all_batches = glob.glob('/Path/to/Patient/Reports/Batch*')
    max_batch = 0
    for b in all_batches:
        batch_num = int(re.search(r'\d+', b).group())
        if batch_num > max_batch:
            max_batch = batch_num

    batch_nums = range(1, max_batch + 3)

    print('Scraping Reports...')

    # Loop through all report directories
    for batch in batch_nums:

        if batch == 8:
            continue

        batch_name = 'Batch ' + str(batch)
        if batch == max_batch + 1:
            batch_name = 'GeneDX'
        elif batch == max_batch + 2:
            batch_name = 'Sanger only'
        print('Scraping ' + batch_name)

        report_types = pd.Series(glob.glob('/Path/to/Patient/Reports/' + batch_name + '/*'))

        # Flips between positive and negative reports
        for i in range(1, flip_reports):

            # Switch between positive and negative reports
            if i == 1:
                report_dir = report_types[report_types.str.contains('Positive', case = False)]
                if report_dir.empty:
                    continue
                report_dir = report_dir.reset_index(drop = True)[0]

                dir_name_doc = report_dir + '/*.docx'
                fnames = glob.glob(dir_name_doc)
                fnames = [x for x in fnames if "~$" not in x]
            else:
                report_dir = report_types[report_types.str.contains('Negative', case = False)]
                if report_dir.empty:
                    continue
                report_dir = report_dir.reset_index(drop = True)[0]
                dir_name_pdf = report_dir + '/*.pdf'
                fnames = glob.glob(dir_name_pdf)

            # Loop through all files in the directory
            for file in fnames:
                # print(file)

                # Get sample_id from the file name
                sample_id = file.replace(report_dir + '/', '')
                sample_id = sample_id[0:7]
                sample_id = "-".join([sample_id[a:a + 2] for a in range(0, len(sample_id), 2)])

                # If we're doing Positive Reports
                if i == 1:
                    doc = docx2python(file)

                    # Split up doc text by periods or new lines
                    text_list = re.split('\. |\n', doc.text)
                    text_list = pd.Series(text_list)

                    # Check to see if this is a CMA Positive report
                    if not text_list[text_list.str.contains('COPY NUMBER VARIANT')].empty:
                        # CMA Positive Reports don't have tables, need to scrape the text in paragraphs

                        # first find 'DNA Change.' It'll be the sentence with the copy number info
                        dna_change = text_list[text_list.str.contains('CMA', case = False)]
                        dna_change = dna_change[dna_change.str.contains('showed')
                                                | dna_change.str.contains('identified')
                                                | dna_change.str.contains('detected')]
                        dna_change = dna_change.str.strip()
                        dna_change = dna_change[dna_change != 'CMA detected copy number variant:']
                        dna_change.rename('DNA_Change')
                        dna_change = dna_change.to_frame()
                        dna_change['sample_id'] = sample_id

                        # find sentence with associated disease. look for 'OMIM'
                        disease = text_list[text_list.str.contains('OMIM', case = False)]
                        disease = disease[~disease.str.contains('OMIM \*')]
                        disease = disease.str.strip()
                        disease.rename('Associated_Disease')
                        disease = disease.to_frame()
                        disease['sample_id'] = sample_id

                        # find sentence with pathogenicity classification
                        pathogenicity = ''
                        classification = text_list[text_list.str.contains('pathogenic | likely pathogenic | uncertain significance | likely benign | benign',
                                                                          case = False)]
                        classification = classification.str.strip()
                        classification.reset_index(drop = True, inplace = True)
                        # classification = classification[classification != 'Likely Benign'].reset_index(drop = True)
                        if not classification.empty:

                            # get actual classification from the sentence
                            line = classification[0].lower()
                            if 'likely pathogenic' in line:
                                pathogenicity = 'Likely Pathogenic'
                            elif 'pathogenic' in line:
                                pathogenicity = 'Pathogenic'
                            elif 'likely benign' in line:
                                pathogenicity = 'Likely Benign'
                            elif 'benign' in line:
                                pathogenicity = 'Benign'
                            elif 'uncertain significance' in line:
                                pathogenicity = 'Uncertain Significance'

                        df = dna_change.merge(disease, on = 'sample_id')
                        df['Classification'] = pathogenicity
                        df.rename(columns = {'0_x': 'CMA Result', '0_y': 'Associated Disease'}, inplace = True)
                        df = df[['sample_id', 'CMA Result', 'Associated Disease', 'Classification']]
                        cma = cma.append(df, sort = False)

                        continue

                    # Deals with all non-CMA positive reports
                    else:
                        # Loop through all tables in the file
                        counter = 1
                        has_table = False

                        for table in doc.body:
                            has_table = True
                            df = pd.DataFrame(table)

                            # check if Gene is in the first cell
                            if 'Gene' in df.iloc[0, 0]:
                                # df.to_csv(outdir + sample_id + '_raw' + str(counter) + '.csv', index=False)
                                counter = counter + 1

                                num_rows = df.shape[0]
                                num_cols = df.shape[1]

                                # Get all scraped text out of the lists into strings
                                for row in range(0, num_rows):
                                    for col in range(0, num_cols):
                                        df.at[row, col] = ' '.join(df.at[row, col])

                                # Clean up the text
                                df = df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
                                df = df.apply(lambda x: x.str.replace('- ', '-') if x.dtype == "object" else x)
                                df = df.apply(lambda x: x.str.replace('  ', ' ') if x.dtype == "object" else x)

                                # Assign first row to column names
                                df.columns = df.loc[0]
                                df = df.loc[1:]
                                df.reset_index(drop=True, inplace=True)

                                # Some reports don't have detailed info
                                # if don't, then just write it out
                                if 'DNA Change' not in df.columns:
                                    # Sometimes the Benign/Likely Benign charts get picked up... ignore them
                                    if 'Benign' in df.columns:
                                        continue
                                    else:
                                        df.insert(0, 'sample_id', sample_id)
                                        df.rename(columns = {'Variant': 'DNA Change', 'Inheritance': 'Disease Inheritance'}, inplace = True)
                                        df['Protein Change'] = None

                                        for j in range(0, df.shape[0]):
                                            # Pull actual dna change out of the text
                                            dna_change = df.at[j, 'DNA Change'].split(' ')
                                            if df.at[j, 'DNA Change'].startswith('chr'):
                                                df.at[j, 'DNA Change'] = dna_change[3]
                                                df.at[j, 'Protein Change'] = dna_change[4].strip('()')

                                        variants = variants.append(df, sort=False)

                                # "Normal" format
                                else:
                                    if 'Gene' in df.columns:

                                        # Clean up column name inconsistencies
                                        if 'OMIM #' in df.columns:
                                            df.rename(columns={'OMIM #': 'OMIM'}, inplace=True)
                                        if 'Associated disease' in df.columns:
                                            df.rename(columns={'Associated disease': 'Associated Disease'}, inplace=True)
                                        if 'Disease inheritance' in df.columns:
                                            df.rename(columns={'Disease inheritance': 'Disease Inheritance'}, inplace=True)
                                        if 'Inheritance' in df.columns:
                                            df.rename(columns={'Inheritance': 'Disease Inheritance'}, inplace=True)

                                        if df.shape[0] > 1:
                                            # Loops through all rows and matches missing gene info or missing disease info
                                            for row_ind in range(1, df.shape[0]):
                                                # fills in empty Gene info from row above
                                                if df.at[row_ind, 'Gene'] == '':
                                                    df.at[row_ind, 'Gene'] = df.at[row_ind - 1, 'Gene']
                                                    df.at[row_ind, 'DNA Change'] = df.at[row_ind - 1, 'DNA Change']
                                                    df.at[row_ind, 'Protein Change'] = df.at[row_ind - 1, 'Protein Change']
                                                    df.at[row_ind, 'Zygosity'] = df.at[row_ind - 1, 'Zygosity']
                                                    df.at[row_ind, 'Classification'] = df.at[row_ind - 1, 'Classification']

                                                # fills in empty disease info from row above
                                                if df.at[row_ind, 'Associated Disease'] == '':
                                                    # Deals with "C9orf66, DOCK8" vs. "DOCK8" issue
                                                    # Checks if the gene above matches the current gene
                                                    prev_genes = df.at[row_ind - 1, 'Gene'].split(', ')
                                                    curr_genes = df.at[row_ind, 'Gene'].split(', ')
                                                    shared_gene = list(set(prev_genes) & set(curr_genes))

                                                    # if there is a shared gene between the rows, link the missing disease info and pathogenicity from row above
                                                    if len(shared_gene) > 0:
                                                        df.at[row_ind, 'Associated Disease'] = df.at[
                                                            row_ind - 1, 'Associated Disease']
                                                        df.at[row_ind, 'OMIM'] = df.at[row_ind - 1, 'OMIM']
                                                        df.at[row_ind, 'Disease Inheritance'] = df.at[
                                                            row_ind - 1, 'Disease Inheritance']
                                                        df.at[row_ind, 'Classification'] = df.at[row_ind - 1, 'Classification']

                                    df.insert(0, 'sample_id', sample_id)
                                    variants = variants.append(df, sort=False)

                        if not has_table:
                            variants.append({
                                'sample_id': sample_id, 'Gene': 'Unable to scrape report', 'DNA Change': 'Unable to scrape report',
                                'Protein Change': 'Unable to scrape report', 'Zygosity': 'Unable to scrape report',
                                'Classification': 'Unable to scrape report', 'Associated Disease': 'Unable to scrape report',
                                'OMIM': 'Unable to scrape report', 'Disease Inheritance': 'Unable to scrape report'},
                                ignore_index=True, sort = False)

                # Handle negative reports
                else:
                    neg_rpt = pd.DataFrame(columns = ['sample_id', 'Batch'])
                    neg_rpt['sample_id'] = [sample_id]
                    neg_rpt['Batch'] = [batch_name]
                    neg_patients = neg_patients.append(neg_rpt, sort = False)

    print('Scraped all reports.')

    # Add Patient, database IDs to variants
    variants.drop_duplicates(inplace = True)
    variants = variants.merge(sample_id_dict, on='sample_id', how='left')
    variants = variants[['Patient_ID', 'database_ID', 'Patient_Family_ID', 'database_Family_ID', 'Batch_Received', 'sample_id', 'Patient_Name', 'Gene', 'DNA Change',
                         'Protein Change', 'Zygosity', 'Classification', 'Associated Disease', 'OMIM', 'Disease Inheritance']]
    neg_patients = neg_patients.merge(sample_id_dict, on = 'sample_id', how = 'left')

    # Isolate all positive reports that weren't able to be scraped
    unscraped = variants[variants['Gene'].str.match('Unable to scrape report', na = False)]
    variants = variants[~variants['Gene'].str.match('Unable to scrape report', na = False)]

    cma = cma.merge(sample_id_dict, on = 'sample_id', how = 'left')
    cma.drop_duplicates(inplace = True)
    cma = cma[['Patient_ID', 'database_ID', 'Patient_Family_ID', 'database_Family_ID',
               'Sample_ID', 'Patient_Name', 'CMA Result', 'Associated Disease', 'Classification']]
    cma.to_csv(outdir + 'all_cma.csv', index = False)

    print('Matching database variants with Patient Report scraped variants...')
    variants_merged = add_database_variants(variants, database_combined)
    variants_merged.to_csv(outdir + 'all_variants_database.csv', index = False)

    print('Finding "Patient Report" tagged variants that are missing Positive Reports...')
    missing_positive_rpt = check_patient_report_variants(variants_merged, database_combined, sample_id_dict)

    print('Finding "Inconclusive negative report" tagged variants that are missing Negative Reports...')
    missing_negative_rpt = check_inconclusive_neg_variants(neg_patients, database_combined, sample_id_dict)

    print('Double checking secondary findings...')
    secondary_findings, variants_merged = get_secondary_findings(variants_merged)

    missing_in_database = variants_merged[variants_merged['Chrom'].isnull()]
    matched_in_database = variants_merged[~variants_merged['Chrom'].isnull()]

    # Write out all dfs into one Excel spreadsheet
    writer = pd.ExcelWriter(outdir + 'scraping_results.xlsx', engine = 'xlsxwriter')
    variants_merged.to_excel(writer, sheet_name = 'All Scraped Pos Rpt Variants', index = False)
    matched_in_database.to_excel(writer, sheet_name = 'Pos Rpt variants in database', index = False)
    missing_in_database.to_excel(writer, sheet_name = 'Pos Rpt variants not in database', index = False)
    missing_positive_rpt.to_excel(writer, sheet_name = 'Patient Rpt variants no Pos Rpt', index = False)
    unscraped.to_excel(writer, sheet_name = 'Cannot scrape Pos Rpt', index = False)
    if include_neg:
        neg_patients.to_excel(writer, sheet_name = 'All patients w Negative Reports', index = False)
        missing_negative_rpt.to_excel(writer, sheet_name = 'Inconc Neg variants no Neg Rpt', index = False)
    cma.to_excel(writer, sheet_name = 'CMA Pos Rpts', index = False)
    secondary_findings.to_excel(writer, sheet_name = 'Secondary Findings', index = False)

    writer.save()


if __name__ == '__main__':
    main()


