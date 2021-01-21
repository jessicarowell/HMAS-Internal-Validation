#!/usr/bin/env python

import logging, sys, os, shutil, subprocess, argparse, sqlite3
from pathlib import Path
import pandas as pd

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
    db_prefix = os.path.splitext(reference)[0]
    db_name = db_prefix + "_db"

    p = subprocess.run([makeblastdb, '-in', reference, '-out', db_name, '-dbtype', dbtype])

    if p.returncode != 0:
        print("makeblastdb could not be executed.  Error encountered.")
        print(p.returncode)
    else:
        print("makeblastdb run successfully.")
        return(db_name)

    # Put the date/time in log name
    # Log the errors currently printed to sdout


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

    p = subprocess.run([blastn, '-db', db, '-query', fasta, \
                                '-outfmt', '6 qseqid sseqid qlen slen evalue ndent mismatch', \
                                '-max_target_seqs', str(maxhits), '-out', blast_out])

    if p.returncode != 0:
        print("blastn could not be executed.  Error encountered.")
        print(p.returncode)
    else:
        print("blastn run successfully.")
        print("format: query ID, subject ID, query length, subj length, e val, # mismatches")
        return(blast_out)

    # Log errors currently printed to stdout

def create_connection(db_file):
    #if Path(db_file).is_file():
    #    db_init = db_file + '_sqlDB'
    #else:
    #    db_init = db_file

    conn = None
    try:
        conn = sqlite3.connect(db_file)
        print(sqlite3.version)
    except Error as e:
        print(e)
    finally:
        if conn:
            conn.close()
        return(db_file)

def get_data(filepath, sep, colnames):
        return(pd.read_csv(filepath, sep = sep, names = colnames))


def main():

    parser = argparse.ArgumentParser(description = 'Run Blast on some data.')
    parser.add_argument('-f', '--query_fasta', metavar = '', required = True, \
                        help = 'Specify query fasta')
    parser.add_argument('-r', '--reference_fasta', metavar = '', required = True, \
                        help = 'Specify data to blast against')
    parser.add_argument('-g', '--group_file', metavar = '', required = True, \
                        help = 'Specify name of group file')
    parser.add_argument('-n', '--name_file', metavar = '', required = True, \
                         help = 'Specify name of name file')
    parser.add_argument('-o', '--outfile', metavar = '', required = True, \
                        help = 'Specify identifying prefix for output files')
    args = parser.parse_args()

    # Generate log file
    logger = create_log(args.outfile + ".log")

    # Check for files
    query_fasta = file_exists(args.query_fasta)
    reference_fasta = file_exists(args.reference_fasta)
    group_file = file_exists(args.group_file)
    name_file = file_exists(args.name_file)

    path_to_makeblastdb = cmd_exists('makeblastdb')
    path_to_blastn = cmd_exists('blastn')

    print(f"All files exist. All commands are executable.")

    db = make_blast_db(reference_fasta, 'nucl', path_to_makeblastdb)
    blast_file = run_blast(db, query_fasta, args.outfile, path_to_blastn, 20)

    print("Building database")
    db_init = create_connection(args.outfile + '.db')

    # Actions on name file (get unique sequences and num of consensus seqs)
    ncolnames = ["seq", "consensusSeqs"]
    n = get_data(name_file, '\t', ncolnames)
    n['nseqs'] = n.consensusSeqs.str.split(',', expand = False).agg(len)
    n = n.drop(columns = ['consensusSeqs'])

    logger.info(f"[OUT] {len(n.index)} unique sequences remain after QC steps performed.")
    logger.info(f"[OUT] These consensus sequences comprise {n.nseqs.sum()} total sequences.")
    logger.info(f"[OUT] The minimum number of identical sequences per consensus sequence is {n.nseqs.min()}")
    #logger.info(f"""[OUT] The max is {n.nseqs.max()}, the median is {n.nseqs.median()}, 
    #           and the mean is {n.nseqs.mean()} (stdev {n.nseqs.std()})""")

    logger.info(f"[OUT] Summary statistics on the number of identical sequences per consensus sequence: \n {n.nseqs.describe()}")
    # Make a plot showing distribution of num_consensus_sequences

    # Actions on group file (get primers for all unique sequences)
    gcolnames = ["seq", "sample_primer"]
    g = get_data(group_file, '\t', gcolnames)
    conn = sqlite3.connect(db_init) 
    c = conn.cursor()
    try:
        g.to_sql('groups', conn, if_exists = 'fail')
    except ValueError:
        print("A table with this name already exists. Using that one.")
         
    ids = n['seq'].tolist()
    select = f"SELECT seq, sample_primer FROM groups WHERE seq IN ({','.join('?' * len(ids))})"
    rows = c.execute(select, ids).fetchall()    
    rows_df = pd.DataFrame(rows, columns = ['seq', 'sample_primer'])
    full = pd.merge(n, rows_df, on = "seq")
    full['primer'] = full['sample_primer'].str.split('.').str[1]
    full = full.drop(columns = ['sample_primer']) 

    # Get the blast output and compare it to the expected
    bcolnames = ["seq", "primer", "query_len", "subj_len", "eval", "mismatch"]
    b = get_data(blast_file, '\t', bcolnames)
    
    # Want to add max num hits to the log file but it's initialized inside function
    #logger.info("[OUT] Blastn run with max target number of sequences set to {maxhits}.")
    logger.info(f"[OUT] {len(b.index)} total Blast hits against the high-quality, unique sequences")
    logger.info(f"[OUT] {b[b.mismatch < 2].shape[0]} total Blast hits with 0 or 1 mismatches.") 


    final = pd.merge(full, b, on = ["seq", "primer"], how = "left")
    nomatch = final[final.isnull().any(axis=1)]
    hits = pd.merge(full, b, on = ["seq", "primer"], how = "inner")

    logger.info(f"[OUT] There are {len(hits.index)} sequences with a blast hit to the primer from the design file.")
    logger.info(f"[OUT] Of these, {hits[hits.mismatch < 2].shape[0]} have 0 or 1 mismatches.")

    hits.to_csv(args.outfile + '_perfect_matches.txt', sep='\t', index = False, float_format = "%.2E")
    nomatch.to_csv(args.outfile + '_no_perfect_matches.txt', columns = ['seq','nseqs','primer'], sep='\t', index = False)
    logger.info(f"[OUT] The perfect matches can be found here: {args.outfile}_perfect_matches.txt ")
    logger.info(f"[OUT] The sequences without a perfect match can be found here: {args.outfile}_no_perfect_matches.txt ")    
    # Wording of the file name ("no matches") is confusing
    m = b.seq.isin(hits.seq)
    m2 = b.loc[m]
    m3 = b.loc[~m]
    m2.to_csv(args.outfile + '_all_matches.txt', sep='\t', index = False, float_format = "%.2E")
    m3.to_csv(args.outfile + '_mismatches.txt', sep='\t', index = False, float_format = "%.2E")

    logger.info("[OUT] The final reads that matched perfectly to their intended primers also matched to other primers in the design file")
    logger.info(f"[OUT] A list of all such matches can be found here: {args.outfile}_all_matches.txt")
    logger.info("[OUT] Some reads that did not match their intended primers did match other primers")
    logger.info(f"[OUT] A list of those can be found here: {args.outfile}_mismatches.txt")

    n = pd.merge(full, b, on = ["seq"], how = "left")
    n2 = n[n.isnull().any(axis=1)]
    
    n2 = n[n.isnull().any(axis=1)].drop(labels = ['primer_y','query_len','subj_len','eval','mismatch'], \
           axis = 1).rename(columns ={'primer_x':'primer'})
    n2.to_csv(args.outfile + '_orphans.txt', columns = ['seq','nseqs','primer'], sep='\t', index = False)
    logger.info("[OUT] Some reads did not hit any of the amplicons in the design file.")
    logger.info(f"[OUT] A list of those can be found here: {args.outfile}_orphans.txt")
    
    

    print("Validation completed.  Please see output and log files for details.")

if __name__ == "__main__":
    main()

