#!/usr/bin/env python3
"""
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/DASCRUBBER_wrapper

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not,
see <http://www.gnu.org/licenses/>.
"""

import argparse
import os
import sys
import shutil
import gzip
import itertools
import random
import subprocess


def main():
    args = get_arguments()
    check_tools_exist()
    check_files_and_directories(args)
    starting_dir = os.getcwd()
    create_temp_dir(args.tempdir)

    try:
        read_name_dict, read_comment_dict, depth = \
            rename_reads_with_fake_pacbio_names(args.input_reads, args.genome_size)
    except AssertionError:
        sys.exit('Error parsing input read file')

    create_dazzler_db(args.dbsplit_options)
    align_reads(args.daligner_options)
    mask_repeats_1(args.repmask_options, depth, args.repeat_depth)
    mask_repeats_2(args.datander_options)
    mask_repeats_3(args.tanmask_options)
    align_reads_with_masking(args.daligner_options)
    estimate_coverage(args.dascover_options)
    intrinsic_quality(args.dasqv_options, depth)
    trim(args.dastrim_options)
    patch(args.daspatch_options)
    new_db(args.dasedit_options)
    extract_reads()
    output_reads(read_name_dict, read_comment_dict)

    os.chdir(starting_dir)
    if not args.keep:
        delete_temp_dir(args.tempdir)

    print_blank_line()


def check_tools_exist():
    missing_tools = []

    def check_tool(tool_name):
        if shutil.which(tool_name) is None:
            return [tool_name]
        else:
            return []

    missing_tools += check_tool('fasta2DB')
    missing_tools += check_tool('DBsplit')
    missing_tools += check_tool('daligner')
    missing_tools += check_tool('REPmask')
    missing_tools += check_tool('datander')
    missing_tools += check_tool('TANmask')
    missing_tools += check_tool('DAScover')
    missing_tools += check_tool('DASqv')
    missing_tools += check_tool('DAStrim')
    missing_tools += check_tool('DASpatch')
    missing_tools += check_tool('DASedit')
    missing_tools += check_tool('DB2fasta')

    if missing_tools:
        sys.exit('Error: could not find tool' + ('' if len(missing_tools) == 1 else 's') +
                 ': ' + ', '.join(missing_tools))


def check_files_and_directories(args):
    if not os.path.isfile(args.input_reads):
        sys.exit('Error: input read file does not exist')
    if os.path.isdir(args.tempdir):
        sys.exit('Error: temporary directory already exists')


def get_arguments():
    parser = argparse.ArgumentParser(description='A wrapper tool for the DASCRUBBER pipeline for '
                                                 'scrubbing (trimming and chimera removal) of long '
                                                 'read sets (PacBio or ONT reads)',
                                     add_help=False)

    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument('-i', '--input_reads', type=str, required=True,
                               help='input set of long reads to be scrubbed')
    required_args.add_argument('-g', '--genome_size', type=str, required=True,
                               help='approximate genome size (examples: 3G, 5.5M or 800k), used '
                                    'to determine depth of coverage')

    directory_args = parser.add_argument_group('Optional arguments')
    directory_args.add_argument('-d', '--tempdir', type=str,
                                help='path of directory for temporary files (default: use a '
                                     'directory in the current location named dascrubber_temp_PID '
                                     'where PID is the process ID)')
    directory_args.add_argument('-k', '--keep', action='store_true',
                                help='keep the temporary directory (default: delete the temporary '
                                     'directory after scrubbing is complete)')
    directory_args.add_argument('-r', '--repeat_depth', type=float, default=3.0,
                                help='REPmask will be given a repeat threshold of this depth, '
                                     'relative to the overall depth (e.g. if 3, then regions with '
                                     '3x the base depth are considered repeats) (default: 3)')

    resource_args = parser.add_argument_group('Command options',
                                              description='You can specify additional options for '
                                                          'each of the Dazzler commands if you do '
                                                          'not want to use the defaults (example: '
                                                          '--daligner_options="-M80")')
    resource_args.add_argument('--dbsplit_options', type=str)
    resource_args.add_argument('--daligner_options', type=str)
    resource_args.add_argument('--repmask_options', type=str)
    resource_args.add_argument('--datander_options', type=str)
    resource_args.add_argument('--tanmask_options', type=str)
    resource_args.add_argument('--dascover_options', type=str)
    resource_args.add_argument('--dasqv_options', type=str)
    resource_args.add_argument('--dastrim_options', type=str)
    resource_args.add_argument('--daspatch_options', type=str)
    resource_args.add_argument('--dasedit_options', type=str)

    help_args = parser.add_argument_group('Help')
    help_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                           help='show this help message and exit')

    args = parser.parse_args()
    args.input_reads = os.path.abspath(args.input_reads)

    if args.tempdir is None:
        args.tempdir = 'dascrubber_temp_' + str(os.getpid())
    args.tempdir = os.path.abspath(args.tempdir)

    if args.repeat_depth <= 1.0:
        sys.exit('Error: repeat depth must be greater than 1')

    args.genome_size = parse_genome_size(args.genome_size)

    def process_extra_options(extra_options):
        if extra_options is not None:
            return extra_options.split(' ')
        else:
            return []

    args.dbsplit_options = process_extra_options(args.dbsplit_options)
    args.daligner_options = process_extra_options(args.daligner_options)
    args.repmask_options = process_extra_options(args.repmask_options)
    args.datander_options = process_extra_options(args.datander_options)
    args.tanmask_options = process_extra_options(args.tanmask_options)
    args.dascover_options = process_extra_options(args.dascover_options)
    args.dasqv_options = process_extra_options(args.dasqv_options)
    args.dastrim_options = process_extra_options(args.dastrim_options)
    args.daspatch_options = process_extra_options(args.daspatch_options)
    args.dasedit_options = process_extra_options(args.dasedit_options)

    return args


def parse_genome_size(genome_size_str):
    genome_size_str = genome_size_str.lower()
    try:
        last_char = genome_size_str[-1]
        if last_char == 'g':
            value_str = genome_size_str[:-1]
            multiplier = 1000000000
        elif last_char == 'm':
            value_str = genome_size_str[:-1]
            multiplier = 1000000
        elif last_char == 'k':
            value_str = genome_size_str[:-1]
            multiplier = 1000
        else:
            value_str = genome_size_str
            multiplier = 1
        if '.' in value_str:
            genome_size = int(round(float(value_str) * multiplier))
        else:
            genome_size = int(value_str) * multiplier
    except (ValueError, IndexError):
        sys.exit('Error: could not parse genome size')
    if genome_size < 1:
        sys.exit('Error: genome size must be a positive value')
    elif genome_size < 100:
        print_warning('genome size is very small (' + int_to_str(genome_size) + ' bases). '
                      'Did you mean to use a suffix (G, M, k)?')
    elif genome_size > 100000000000:
        print_warning('genome size is very large (' + int_to_str(genome_size) + ' bases). '
                      'Is that a mistake?')
    return genome_size


def create_temp_dir(tempdir):
    print_header('Creating temporary directory')
    print_command(['mkdir', tempdir])
    os.makedirs(tempdir)
    print_command(['cd', tempdir])
    os.chdir(tempdir)
    print_blank_line()


def rename_reads_with_fake_pacbio_names(read_filename, genome_size):
    print_header('Processing and renaming reads')
    files_before = list(os.listdir('.'))

    read_name_dict, read_comment_dict = {}, {}
    read_names = set()

    file_type = get_sequence_file_type(read_filename)
    open_func = get_open_func(read_filename)
    counter = itertools.count(start=0)
    base_total = 0

    with open_func(read_filename, 'rt') as seq_file:
        with open('renamed_reads.fasta', 'wt') as renamed_reads:
            read_num, update_interval = 0, 1

            for header in seq_file:
                try:
                    if file_type == 'FASTA':
                        assert header[0] == '>'
                    elif file_type == 'FASTQ':
                        assert header[0] == '@'
                    header = header[1:].strip()
                    assert len(header) > 0
                except AssertionError:
                    sys.exit('Error: failed to parse read header')

                header_parts = header.split(' ', 1)
                old_name = header_parts[0]
                if old_name in read_names:
                    sys.exit('Error: duplicate read name: ' + old_name)
                read_names.add(old_name)
                try:
                    old_comment = header_parts[1]
                except IndexError:
                    old_comment = None

                seq = next(seq_file).strip()
                read_num = next(counter)
                new_header = 'reads/' + str(read_num) + '/0_' + str(len(seq))

                if read_num == 0 or read_num % update_interval == 0:
                    print('\rReads: ' + int_to_str(read_num+1), file=sys.stderr, end='')
                    update_interval = random.randint(100, 110)

                read_name_dict[read_num] = old_name
                read_comment_dict[read_num] = old_comment

                renamed_reads.write('>')
                renamed_reads.write(new_header)
                renamed_reads.write('\n')
                renamed_reads.write(seq)
                renamed_reads.write('\n')

                base_total += len(seq)

                if file_type == 'FASTQ':
                    next(seq_file)
                    next(seq_file)

    depth = base_total / genome_size
    print('\rReads:', int_to_str(read_num+1), file=sys.stderr)
    print('Total bases:', int_to_str(base_total), file=sys.stderr)
    print('Depth of coverage:', float_to_str(depth, 1), file=sys.stderr)
    print_new_files(files_before)
    print_blank_line()
    return read_name_dict, read_comment_dict, depth


def create_dazzler_db(dbsplit_options):
    print_header('Creating Dazzler database')
    files_before = list(os.listdir('.'))

    run_command(['fasta2DB', 'reads.db', 'renamed_reads.fasta'])

    if any(x.startswith('-s') for x in dbsplit_options):
        split = []
    else:
        split = ['-s100']

    run_command(['DBsplit'] + split + dbsplit_options + ['reads'])

    print_new_files(files_before)
    print_blank_line()


def align_reads(daligner_options):
    print_header('Read overlap alignment with daligner')
    files_before = list(os.listdir('.'))

    print_command(['mkdir', 'align_temp'])
    os.makedirs('align_temp')

    run_command(['daligner', '-v', '-Palign_temp'] + daligner_options + ['reads', 'reads'])

    print_command(['rm', '-r', 'align_temp'])
    shutil.rmtree('align_temp')

    print_new_files(files_before)
    print_blank_line()


def mask_repeats_1(repmask_options, depth, repeat_depth):
    print_header('Masking repeats with REPmask')
    files_before = list(os.listdir('.'))

    if any(x.startswith('-c') for x in repmask_options):
        threshold = []
    else:
        threshold = ['-c' + str(int(round(depth * repeat_depth)))]

    run_command(['REPmask', '-v'] + threshold + repmask_options + ['reads', 'reads.reads.las'])

    print_new_files(files_before)
    print_blank_line()


def mask_repeats_2(datander_options):
    print_header('Finding tandem repeats with datander')
    files_before = list(os.listdir('.'))

    print_command(['mkdir', 'align_temp'])
    os.makedirs('align_temp')

    run_command(['datander', '-v', '-Palign_temp'] + datander_options + ['reads'])

    print_command(['rm', '-r', 'align_temp'])
    shutil.rmtree('align_temp')

    print_new_files(files_before)
    print_blank_line()


def mask_repeats_3(tanmask_options):
    print_header('Masking tandem repeats with TANmask')
    files_before = list(os.listdir('.'))

    run_command(['TANmask', '-v'] + tanmask_options + ['reads', 'TAN.reads'])

    print_new_files(files_before)
    print_blank_line()


def align_reads_with_masking(daligner_options):
    print_header('Read overlap alignment with daligner (with repeat masking)')
    files_before = list(os.listdir('.'))

    print_command(['mkdir', 'align_temp'])
    os.makedirs('align_temp')

    run_command(['daligner', '-v', '-Palign_temp', '-mrep', '-mtan'] + daligner_options +
                ['reads', 'reads'])

    print_command(['rm', '-r', 'align_temp'])
    shutil.rmtree('align_temp')

    print_new_files(files_before)
    print_blank_line()


def estimate_coverage(dascover_options):
    print_header('Computing estimated genome coverage with DAScover')
    files_before = list(os.listdir('.'))

    run_command(['DAScover', '-v'] + dascover_options + ['reads', 'reads.reads.las'])

    print_new_files(files_before)
    print_blank_line()


def intrinsic_quality(dasqv_options, depth):
    print_header('Finding intrinsic quality values with DASqv')
    files_before = list(os.listdir('.'))

    if any(x.startswith('-c') for x in dasqv_options):
        depth = []
    else:
        depth = ['-c' + str(int(round(depth)))]

    run_command(['DASqv', '-v'] + depth + dasqv_options + ['reads', 'reads.reads.las'])

    print_new_files(files_before)
    print_blank_line()


def trim(dastrim_options):
    print_header('Trimming reads and breaking chimeras with DAStrim')
    files_before = list(os.listdir('.'))

    run_command(['DAStrim', '-v'] + dastrim_options + ['reads', 'reads.reads.las'])

    print_new_files(files_before)
    print_blank_line()


def patch(daspatch_options):
    print_header('Patching low quality segments with DASpatch')
    files_before = list(os.listdir('.'))

    run_command(['DASpatch', '-v'] + daspatch_options + ['reads', 'reads.reads.las'])

    print_new_files(files_before)
    print_blank_line()


def new_db(dasedit_options):
    print_header('Building new database of scrubbed reads with DASedit')
    files_before = list(os.listdir('.'))

    run_command(['DASedit', '-v'] + dasedit_options + ['reads', 'patched_reads'])

    print_new_files(files_before)
    print_blank_line()


def extract_reads():
    print_header('Extracting scrubbed reads')
    files_before = list(os.listdir('.'))

    print_command(['mv', 'renamed_reads.fasta', 'temp.fasta'])
    os.rename('renamed_reads.fasta', 'temp.fasta')

    run_command(['DB2fasta', '-vU', 'patched_reads'])

    print_command(['mv', 'renamed_reads.fasta', 'scrubbed_reads.fasta'])
    os.rename('renamed_reads.fasta', 'scrubbed_reads.fasta')

    print_command(['mv', 'temp.fasta', 'renamed_reads.fasta'])
    os.rename('temp.fasta', 'renamed_reads.fasta')

    print_new_files(files_before)
    print_blank_line()


def output_reads(read_name_dict, read_comment_dict):
    print_header('Outputting scrubbed reads to stdout')

    def print_read(new_name, seq):
        parts = new_name.split('/')
        read_num = int(parts[1])
        read_range = parts[2]

        # Print read to stdout
        header = '>' + read_name_dict[read_num] + '/' + read_range
        if read_comment_dict[read_num] is not None:
            header += ' ' + read_comment_dict[read_num]
        print(header)
        print(seq)

        return len(seq)

    counter = itertools.count(start=0)
    base_total = 0
    new_read_num, update_interval = 0, 1

    with open('scrubbed_reads.fasta', 'rt') as scrubbed_reads:
        name = ''
        sequence = ''
        for line in scrubbed_reads:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new read
                if name:
                    new_read_num = next(counter)
                    base_total += print_read(name, sequence)
                    sequence = ''

                    if new_read_num == 0 or new_read_num % update_interval == 0:
                        print('\rReads: ' + int_to_str(new_read_num + 1), file=sys.stderr, end='')
                        update_interval = random.randint(100, 110)

                name = line[1:]
            else:
                sequence += line
        if name:
            new_read_num = next(counter)
            base_total += print_read(name, sequence)

    print('\rReads:', int_to_str(new_read_num+1), file=sys.stderr)
    print('Total bases:', int_to_str(base_total), file=sys.stderr)
    print_blank_line()


def delete_temp_dir(tempdir):
    print_header('Deleting temporary directory')
    print_command(['rm', '-r', tempdir])
    shutil.rmtree(tempdir)
    print_blank_line()


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for filetype, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = filetype
    if compression_type == 'bz2':
        sys.exit('Error: bzip2 format not supported')
    if compression_type == 'zip':
        sys.exit('Error: zip format not supported')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def get_sequence_file_type(filename):
    """
    Determines whether a file is FASTA or FASTQ.
    """
    if not os.path.isfile(filename):
        sys.exit('Error: could not find ' + filename)
    open_func = get_open_func(filename)

    with open_func(filename, 'rt') as seq_file:
        try:
            first_char = seq_file.read(1)
        except UnicodeDecodeError:
            first_char = ''

    if first_char == '>':
        return 'FASTA'
    elif first_char == '@':
        return 'FASTQ'
    else:
        raise ValueError('Error: file is neither FASTA or FASTQ')


END_FORMATTING = '\033[0m'
UNDERLINE = '\033[4m'
BOLD = '\033[1m'
DIM = '\033[2m'
BLUE = '\033[34m'
GREEN = '\033[32m'
RED = '\033[31m'


# Not all terminals support dim.
try:
    colours = int(subprocess.check_output(['tput', 'colors']).decode().strip())
    if colours <= 8:
        DIM = ''
except (ValueError, subprocess.CalledProcessError, FileNotFoundError, AttributeError):
    DIM = ''


def print_blank_line():
    print('', file=sys.stderr)


def print_header(text):
    print_blank_line()
    print(BOLD + UNDERLINE + text + END_FORMATTING, file=sys.stderr)


def print_command(command):
    text = ' '.join(command)
    print(GREEN + text + END_FORMATTING, file=sys.stderr)


def run_command(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print_command(command)
    full_output = []

    while True:
        output = process.stdout.readline().decode()
        full_output.append(output)
        if output == '' and process.poll() is not None:
            break
        output = output.rstrip()
        if output:
            print('  ' + DIM + output + END_FORMATTING, file=sys.stderr)

    return_code = process.poll()
    if return_code != 0:
        sys.exit('Error: command failed')
    return full_output


def print_warning(text):
    print_blank_line()
    print(RED + BOLD + 'WARNING: ' + text + END_FORMATTING, file=sys.stderr)
    print_blank_line()


def print_new_files(files_before):
    files_after = list(os.listdir('.'))
    new_files = sorted(f for f in files_after if f not in files_before)
    if new_files:
        if len(new_files) == 1:
            print('New file: ', file=sys.stderr, end='')
        else:
            print('New files: ', file=sys.stderr, end='')
        print(BLUE + ', '.join(new_files) + END_FORMATTING, file=sys.stderr)


def int_to_str(num):
    return '{:,}'.format(num)


def float_to_str(num, decimals):
    if decimals == 0:
        return int_to_str(int(round(num)))
    else:
        num_str = '%.' + str(decimals) + 'f'
        num_str = num_str % num
        parts = num_str.split('.')
        before_decimal = parts[0]
        after_decimal = parts[1]
        return int_to_str(int(before_decimal)) + '.' + after_decimal


if __name__ == '__main__':
    main()
