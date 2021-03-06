#!/usr/bin/env python

import logging, sys, os, shutil, subprocess, argparse
from pathlib import Path
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None  # default='warn'

def create_log(filename, filepath = os.getcwd()):
    if filepath == os.getcwd():
	    logname = os.getcwd() + '/' + filename
    else:
        try:
            logname = os.path.abspath(filepath) + '/' + filename
        except:
            print(f"Unable to resolve {filepath}. Creating log in cwd.")
            logname = os.getcwd() + '/' + filename 
    LOG_FORMAT = "%(levelname)s %(asctime)s - %(message)s"
    logging.basicConfig(filename = logname, format = LOG_FORMAT, level = logging.DEBUG)
    logger = logging.getLogger()
    return(logger)


def file_exists(filename):
    """Checks for existence of a file

    Params
    ------
    file: String
        Name of the file

    Returns
    ------
    Path to file

    """
    if Path(filename).is_file():
        print(f"{filename} found.")
        return(Path(filename).resolve())
    else: 
        print(f"{filename} could not be found.")
        sys.exit()


def cmd_exists(command):
    """Checks for existence of a command on user's path

    Params
    ------
    command: String
        Name of the command

    Returns
    ------
    Path to executable command
    
    """
    if shutil.which(command) is None:
        print(f"{command} not found on PATH")
        sys.exit()
    else:
        return(shutil.which(command))

def make_blast_db(reference, dbtype, makeblastdb):
    """Makes a Blast database from a set of reference amplicons

    Params
    ------
    reference: String
        Path to reference amplicon set

    dbtype: String
	Specify nucleotide or protein database

    Returns
    ------

    """
    # Would like to add db_prefix = 'amr2020' so that it's easy to clean up files
    db_name = os.path.splitext(reference)[0]
    p = subprocess.run([makeblastdb, '-in', reference, '-out', db_name, '-dbtype', dbtype])
    if p.returncode != 0:
        print("makeblastdb could not be executed.  Error encountered.")
        print(p.returncode)
    else:
        print("makeblastdb run successfully.")
        return(db_name)


def run_blast(db, fasta, outfile, blastn, maxhits = 10):
    """Runs Blast on a fasta against a reference database and saves output

    Params
    ------
    ref: String
        Path to reference blast database

    fasta: String
        Specify query fasta sequence

    outfile: String
        Specify output file name

    Returns
    ------
    blast_out: String
        Name of the file containing the blast results
    """
    blast_out = outfile + '_blast.out'
    if maxhits != 10:
        try:
            maxhits + 1 - 1
        except TypeError: 
            print(f"You must supply an integer value for max. number of hits. The default is 10.")
            maxhits = 10
        finally:
            maxhits = maxhits
    print(f"Running blast with {maxhits} maximum hits to generate")
    p = subprocess.run(
                [blastn, '-db', db, '-query', fasta, \
                         '-outfmt', \
                         '6 qseqid sseqid qlen slen evalue qcovs pident mismatch', \
                         '-max_target_seqs', str(maxhits), '-out', blast_out])
    if p.returncode != 0:
        print("blastn could not be executed.  Error encountered.")
        print(p.returncode)
    else:
        print("blastn run successfully.")
        print("format: query ID, subject ID, query length, subj length, e val, query cov/subj, percent identical matches, # mismatches")
        return(blast_out) 

    # Log errors currently printed to stdout

def get_data(filepath, sep, colnames = None):
    if colnames is not None:
        return(pd.read_csv(filepath, sep = sep, names = colnames))
    else:
        return(pd.read_csv(filepath, sep = sep))

def main_analysis(sample_keyword, data_df, data_type, primer_df, group_df, logger):

    # get only the Twist control samples (containing all post-QC seqs)
    final = data_df
    final_twist = data_df.loc[data_df['sample'].str.contains(sample_keyword) ]

    conditions = [
        (final_twist.isnull().any(axis=1)), # 0
        (final_twist['pident'] == 100) & (final_twist['cov'] >= 98), 
        (final_twist['pident'] < 100) | (final_twist['cov'] < 98)
    ]

    values = [0,2,1]
    final_twist['match'] = np.select(conditi

    logger.info(f"""[OUT] 
        {final_twist.groupby('match')['nseqs'].sum()[2]:,} total sequences match their expected target with 100% identity and 98% coverage or better.")
        These comprise {final_twist['match'].value_counts()[2]:,} unique sequences.
        {final_twist.groupby('match')['nseqs'].sum()[1]:,} total sequences match their expected target with < 100% identity and/or < 98% coverage.
        These comprise {final_twist['match'].value_counts()[1]:,} unique sequences.
        {final_twist.groupby('match')['nseqs'].sum()[0]:,} total sequences did not hit their expected target at all.
        These comprise {final_twist['match'].value_counts()[0]:,} unique sequences.""")

    final_twist['match'].value_counts() # unique seqs in each group
    final_twist.groupby('match')['nseqs'].sum() # total seqs in each group
    final_twist.groupby('match')['primer'].nunique() # Not that helpful because of overlap between groups

    # if primer has a perfect match, remove it from match=0 
    ids = final_twist.loc[final_twist['match'] == 2, 'primer'].drop_duplicates() # list of ids with perfect match
    index_ids = final_twist[ (final_twist['match'] == 0) & (final_twist['primer'].isin(ids)) ].index
    #final_twist.drop(index_ids, inplace = True)
    newf = final_twist.drop(index_ids)

    logger.info(f"""[OUT] 
        After removing primers that do have a perfect match to their target, there are
        {newf.groupby('match')['nseqs'].sum()[0]:,} total sequences remain that did not hit their target. 
        These comprise {newf['match'].value_counts()[0]:,} unique sequences. """)


    # examine the primers with imperfect hits to their target (match = 1)
    perfect = final_twist.loc[ (final_twist['match'] == 1) & (final_twist['primer'].isin(ids)) & \
                            (final_twist['pident'] == 100) & (final_twist['cov'] >= 98), 'primer'].drop_duplicates()
    middle = final_twist.loc[ (final_twist['match'] == 1) & (final_twist['primer'].isin(ids)) & \
                            (final_twist['pident'] > 99) & (final_twist['cov'] >= 98), 'primer'].drop_duplicates()
    low = final_twist.loc[ (final_twist['match'] == 1) & (final_twist['primer'].isin(ids)) & \
                            (final_twist['pident'] < 100) & (final_twist['cov'] < 98), 'primer'].drop_duplicates()                                             
    logger.info(f"""[OUT]
        Examining the {final_twist['match'].value_counts()[1]:,} unique sequences that hit their intended target with <100% identity and/or <98% coverage, 
        {perfect.shape[0]:,} primers had at least 1 match with 100% ID and >= 98% coverage.
        {middle.shape[0]:,} primers had at least 1 match with >99% ID and >=98% coverage.
        {low.shape[0]:,} primers had a match with <100% ID and <98% coverage.""")

    # 'perfect' is the df we would be looking more closely at. Everything else gets removed. For this data, perfect is empty.

    # if primer has a perfect match, remove from match=1 
    index_ids = final_twist[ (final_twist['match'] == 1) & (final_twist['primer'].isin(ids)) ].index
    newf2 = newf.drop(index_ids)

    logger.info(f"""[OUT] After removing primer pairs with perfect hits from the rest of the data, there are 
        {newf2['match'].value_counts()[2]:,} unique sequences that "perfectly" hit their target (abundance = {newf2.groupby('match')['nseqs'].sum()[2]:,})
        {newf2['match'].value_counts()[1]:,} unique sequences that hit their target with <100% &/or <98% coverage (abundance = {newf2.groupby('match')['nseqs'].sum()[1]:,})
        {newf2['match'].value_counts()[0]:,} unique sequences that did not hit their intended target (abundance = {newf2.groupby('match')['nseqs'].sum()[0]:,}) """)
     
    # Remove any seqs with match = 0 if the primer has an imperfect match
    ids = newf2.loc[(newf2['match'] == 1), 'primer'].drop_duplicates()
    index_ids = newf2[ (newf2['match'] == 0) & (newf2['primer'].isin(ids)) ].index
    newf3 = newf2.drop(index_ids)
    logger.info(f"""[OUT] 
        After removing {len(index_ids)} primer pairs that imperfectly hit their target from the group that did not hit their target,
        {newf3.groupby('match')['primer'].nunique()[2]} primers hit their target with 100% identity and >= 98% coverage
        {newf3.groupby('match')['primer'].nunique()[1]} primers hit their target with <100% and/or <98% coverage
        {newf3.groupby('match')['primer'].nunique()[0]} primers did not hit their target""")

    primers_in_pools = group_df['sample_primer'].str.split('.').str[1].drop_duplicates()

    logger.info(f"[OUT] {len(primers_in_pools)} primers made it past QC filters with at least 10 reads. ")
    logger.info(f"[OUT {newf3['primer'].nunique()} of these primers amplified sequences in the positive control samples.")
    logger.info(f"[OUT] There are {primer_df.shape[0]} total primers on the primer panel used in this assay" )


    t = primers_in_pools.isin(newf3['primer'].drop_duplicates())
    stray = primers_in_pools[~t] # The 11 primers that passed QC but not in Twist samples
    stray_df = final.loc[(final['primer'].isin(stray))]
    stray_df.to_csv('primers_'+sample_keyword+'_passQC_nomatch.txt',sep='\t', index = False, float_format = "%.2E")

    nm = newf3.loc[(newf['match'] == 0)]
    ps = nm['primer'].drop_duplicates()
    newp = primer_df[primer_df['name'].isin(ps)]
    newp.to_csv('primers_'+sample_keyword+'_nomatch.txt',sep='\t', index = False, float_format = "%.2E")

def main():

    parser = argparse.ArgumentParser(description = 'Run HMAS Internal Validation pipeline after the HMAS QC pipeline.')
    parser.add_argument('-f', '--fasta', metavar = '', required = True, help = 'Specify fasta file output from HMAS QC pipeline (should have "final" in the filename).')
    parser.add_argument('-g', '--group', metavar = '', required = True, help = 'Specify group file output from HMAS QC pipeline (should have "final" in the filename).')
    #parser.add_argument('-n', '--name', metavar = '', required = True, help = 'Specify name file output from HMAS QC pipeline (should have "final" in the filename).')
    parser.add_argument('-c', '--count', metavar = '', required = True, help = 'Specify count file output from HMAS QC pipeline (should have "final" in the filename).')
    parser.add_argument('-r', '--reference', metavar = '', required = True, help = 'Specify fasta file containing the positive control targets.')
    parser.add_argument('-p', '--primers', metavar = '', required = True, help = 'Specify primer design file.')
    args = parser.parse_args()

    # primer file: enterics.amr.fixed.genotypes_twistAll.txt (twist pos ctrls) 10932_ORP_19.O1_Design_File_06242019.txt (amr data)
    # reference file: twist_pc_design_amplicons.fasta (twist pos ctrls) amr_design_primers.fasta (amr data)

    # Generate log file
    logger = create_log('amr2020.log')

    # Check for files
    query_fasta = file_exists(args.fasta)
    reference_fasta = file_exists(args.reference)  
    group_file = file_exists(args.groups)
    name_file = file_exists(args.names)
    primer_file = file_exists(args.primers)

    path_to_makeblastdb = cmd_exists('makeblastdb')
    path_to_blastn = cmd_exists('blastn')

    db = make_blast_db(reference_fasta, 'nucl', path_to_makeblastdb)
    blast_file = run_blast(db, query_fasta, 'amr2020', path_to_blastn, 20)

    # Actions on count file
    c = get_data(count_file, '\t')
    c = c[['Representative_Sequence','total']]
    c.columns = ['seq','nseqs']

    # Actions on name file (get unique sequences and num of consensus seqs)
    #ncolnames = ["seq", "sameSeqs"]
    #n = get_data(name_file, '\t', ncolnames)
    #n['nseqs'] = n.sameSeqs.str.split(',', expand = False).agg(len)
    #n = n.drop(columns = ['sameSeqs'])

    # Note that # of unique seqs after QC steps may not be correct?  (Inconsistent with number in fasta file)
    logger.info(f"[OUT] {len(c.index):,} unique sequences remain after QC steps performed.")
    logger.info(f"[OUT] These unique sequences comprise {c.nseqs.sum():,} total sequences.")
    logger.info(f"[OUT] The minimum number of identical sequences per unique sequence is {c.nseqs.min()}")
    logger.info(f"""[OUT] The max is {c.nseqs.max():,}, the median is {c.nseqs.median()}, 
        and the mean is {int(c.nseqs.mean()):,} (stdev {int(c.nseqs.std()):,})""")

    # Actions on group file (get primers for all unique sequences)
    gcolnames = ["seq", "sample_primer"]
    g = get_data(group_file, '\t', gcolnames)

    # Actions on primer file  (for new panel)
    #p['name'] = p['assay_id'] + p['assay_name']

    p = get_data(primer_file, '\t')
    p = p.rename(columns = str.lower)

    # Create a dataframe with all seqs that passed QC, the # of consensus seqs, and associated primer 
    full = pd.merge(c, g, on = 'seq', how = 'left')
    full['primer'] = full['sample_primer'].str.split('.').str[1]
    full['sample'] = full['sample_primer'].str.split('.').str[0]
    full = full.drop(columns = ['sample_primer'])

    # Get the blast output
    bcolnames = ["seq", "primer", "query_len", "subj_len", "eval", "cov", "pident", "mismatch"]
    b = get_data(blast_file, '\t', bcolnames)

    # Merge the blast hits with the full dataset (containing all post-QC seqs)
    final = pd.merge(full, b, on = ['seq', 'primer'], how = 'left') # All post-QC seqs matched to blast hits on exact primer match

    # sample_keyword, data_df, data_type, primer file, group file, logger name
    main_analysis('Twist', final, 'positive controls', p, g, logger)

    logger.info("Validation completed.")

    print("[OUT] Validation completed.  Please see output and log files for details.")

if __name__ == "__main__":
    main()
