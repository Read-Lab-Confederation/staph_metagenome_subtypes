#! /usr/bin/env python
"""Determine the MLST for a set of genomes."""


def run_command(cmd, stdout=False, stderr=False, verbose=False):
    """
    Execute a single command and return STDOUT and STDERR.

    If stdout or stderr are given, output will be written to given file name.
    """
    import subprocess

    cmd = filter(None, cmd)
    if verbose:
        print(' '.join(cmd))
    stdout = open(stdout, 'w') if stdout else subprocess.PIPE
    stderr = open(stderr, 'w') if stderr else subprocess.PIPE
    p = subprocess.Popen(cmd, stdout=stdout, stderr=stderr)

    return p.communicate()


def blast_alleles(input_file, blastn_results, blastdb_path, blastn_path,
                  num_cpu):
    """Blast assembled contigs against MLST blast database."""
    from collections import OrderedDict
    import json
    outfmt = "6 sseqid bitscore slen length gaps mismatch pident evalue"
    alleles = ['arcC', 'aroE', 'glpF', 'gmk', 'pta', 'tpi', 'yqiL']
    results = OrderedDict()
    top_hits = []

    for allele in alleles:
        blastdb = '{0}/{1}'.format(blastdb_path, allele)
        blastn = run_command(
            [blastn_path, '-db', blastdb, '-query', input_file,
             '-outfmt', outfmt, '-max_target_seqs', '1', '-num_threads',
             num_cpu, '-evalue', '10000']
        )
        top_hit = blastn[0].split('\n')[0].split('\t')

        # Did not return a hit
        if not top_hit[0]:
            top_hit = ['0'] * 9
            top_hit[0] = '{0}_0'.format(allele)

        results[allele] = OrderedDict((
            ('sseqid', top_hit[0]),
            ('bitscore', top_hit[1]),
            ('slen', top_hit[2]),
            ('length', top_hit[3]),
            ('gaps', top_hit[4]),
            ('mismatch', top_hit[5]),
            ('pident', top_hit[6]),
            ('evalue', top_hit[7])
        ))

        if float(top_hit[6]) == 100 and top_hit[2] == top_hit[3]:
            top_hits.append(top_hit[0].split('_')[1])
        else:
            top_hits.append('?')

    with open(blastn_results, 'w') as fh:
        json.dump(results, fh, indent=4, separators=(',', ': '))

    return '-'.join(top_hits)

if __name__ == '__main__':
    from os.path import basename
    import glob
    import argparse as ap

    parser = ap.ArgumentParser(
        prog='determine-mlst.py',
        conflict_handler='resolve',
        description=(''))
    group1 = parser.add_argument_group('Options', '')
    group1.add_argument('input', metavar="INPUT_DIR", type=str,
                        help=('Directory of completed genomes to blast.'))
    group1.add_argument('output', metavar="BLAST_OUT", type=str,
                        help=('Directory to output blast results to.'))
    group1.add_argument('blastn', metavar="BLASTN_PATH", type=str,
                        help=('Path to blastn.'))
    group1.add_argument('blastdb', metavar="BLASTDB", type=str,
                        help=('Path to mlst blastdb.'))
    group1.add_argument('mlst', metavar="BLASTDB", type=str,
                        help=('MLST mappings.'))

    args = parser.parse_args()

    STs = {}
    with open(args.mlst, 'r') as fh:
        for line in fh:
            if line.startswith('ST'):
                pass
            else:
                line = line.rstrip()
                cols = line.split('\t')
                STs['-'.join(cols[1:8])] = cols[0]

    for genome in glob.glob('{0}/*.fna'.format(args.input)):
        blastn_result = '{0}/{1}'.format(
            args.output,
            basename(genome).replace('.fna', '.json')
        )
        st = blast_alleles(
            genome, blastn_result, args.blastdb, args.blastn, '23'
        )
        if st in STs:
            print('{0}\t{1}'.format(genome, STs[st]))
        else:
            print('{0}\tUndetermined'.format(genome))
  
