#!/usr/bin/env python
from __future__ import print_function, division
import sys
import numpy as np
import os
import argparse
import uuid
import subprocess
import re
import shlex


def read_xvg(file_handle):
    """Parses XVG file legends and data"""

    _ignored = {'legend', 'view'}
    _re_series = re.compile('s[0-9]+$')
    _re_xyaxis = re.compile('[xy]axis$')

    metadata = {}
    num_data = []

    metadata['labels'] = {}
    metadata['labels']['series'] = []

    for line in file_handle:
        line = line.strip()
        if line.startswith('@'):
            tokens = shlex.split(line[1:])
            if tokens[0] in _ignored:
                continue
            elif tokens[0] == 'TYPE':
                if tokens[1] != 'xy':
                    raise ValueError('Chart type unsupported: \'{0}\'. Must be \'xy\''.format(tokens[1]))
            elif _re_series.match(tokens[0]):
                metadata['labels']['series'].append(tokens[-1])
            elif _re_xyaxis.match(tokens[0]):
                metadata['labels'][tokens[0]] = tokens[-1]
            elif len(tokens) == 2:
                metadata[tokens[0]] = tokens[1]
            else:
                print('Unsupported entry: {0} - ignoring'.format(tokens[0]), file=sys.stderr)
        elif line[0].isdigit():
                num_data.append(map(float, line.split()))

    num_data = list(zip(*num_data))

    if not metadata['labels']['series']:
        for series in range(len(num_data) - 1):
            metadata['labels']['series'].append('')

    return metadata, num_data


def unquote(s):
    return s[1+s.find('"'):s.rfind('"')]


def uncomment(s):
    return s[3+s.find('/*'):s.rfind('*/')-1]


def make_lut(c, nb):
    color = c.split('/*')
    value = unquote(color[1])
    if value == "None":
        value = 0.0
    if value == "Present":
        value = 1.0
    tmp = unquote(color[0]).split(' c ')
    key = c[1:1+nb]
    color = tmp[1].strip()
    return key, (color, float(value))


def process_meta(metadata):
    uncommented = []
    for i in metadata:
        j = uncomment(i)
        if len(j) == 0:
            continue
        if j[0] == ' ':
            uncommented[len(uncommented)-1] += j
        elif ":" in j:
            uncommented.append(j)
    m_dict = dict()
    for i in uncommented:
        index = i.find(':')
        key = i[:index]
        value = i[index+1:].strip()
        if value[0] == '"':
            value = unquote(value)
        if key in m_dict:
            m_dict[key] += ' '+value
        else:
            m_dict[key] = value
    return m_dict


def read_xpm(fhandle, reverse=False):
    # Read in lines until we fidn the start of the array
    metadata = [fhandle.readline()]
    while not metadata[-1].startswith("static char *gromacs_xpm[]"):
        metadata.append(fhandle.readline())

    # The next line will contain the dimensions of the array
    dim = fhandle.readline()
    # There are four integers surrounded by quotes
    nx, ny, nc, nb = [int(i) for i in unquote(dim).split()]

    # The next dim[2] lines contain the color definitions
    # Each pixel is encoded by dim[3] bytes, and a comment
    # at the end of the line contains the corresponding value
    lut = dict([make_lut(fhandle.readline(), nb) for _ in range(nc)])

    colors = []
    values = []

    for i in fhandle:
        if i.startswith("/*"):
            metadata.append(i)
            continue
        j = unquote(i)
        colors.append([lut[j[k:k + nb]][0] for k in range(0, nx, nb)])
        values.append([lut[j[k:k + nb]][1] for k in range(0, nx, nb)])
    metadata = process_meta(metadata)

    if reverse:
        values.reverse()
        colors.reverse()

    return colors, values, metadata, lut


def check_write_ext(choices):
    class Act(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            extension = os.path.splitext(values)[1][1:]
            c = choices
            if type(c) == str:
                c = {c}
            if extension not in c:
                parser.error("file {} doesn't end with one of {}".format(values, choices))
            elif os.path.exists(values):
                if os.path.isdir(values):
                    parser.error("{} is a directory".format(values))
                elif not os.access(values, os.W_OK):
                    parser.error("file {} is not writeable".format(values))
                else:
                    setattr(namespace, self.dest, values)
            else:
                path = os.path.dirname(os.path.abspath(values))
                if not os.access(path, os.W_OK):
                    parser.error("the directory {} is not writeable".format(path))
                else:
                    setattr(namespace, self.dest, values)
    return Act


def delete_if_exists(file_name):
    if os.path.exists(file_name):
        os.remove(file_name)


def backup_if_exists(file_name):
    if os.path.exists(file_name):
        file_name = os.path.abspath(file_name)
        the_dir = os.path.dirname(file_name)
        base = os.path.basename(file_name)
        files = os.listdir(the_dir)
        backups = [i for i in files if i[:len(base)+1] == ('#' + base)]
        if len(backups) == 0:
            number = 1
        else:
            number = max([int(i.split('.')[-1][:-1]) for i in backups])+1
        new_name = os.path.join(the_dir, '#'+base+'.'+str(number)+'#')
        os.rename(file_name, new_name)
        file_name = os.path.join(the_dir, os.path.basename(file_name))
        new_name = os.path.join(the_dir, os.path.basename(new_name))
        sys.stdout.write('\nBack Off! I just backed up {} to {}\n\n'.format(file_name, new_name))


class Analysis:
    def __init__(self, directory):

        # read at init
        self.directory = directory
        self.selected_bonds = None
        self.group_index = None
        self.unique_meta = None
        self.unique_pair = None
        self.consensus = None
        self.time = None
        self.distance = None

        # filled later
        self.missing = set()
        self.found = []
        self.to_do = None
        self.time = None
        self.distances = None
        self.labels = None

        # TODO error management

        list_file = os.path.join(directory, 'hbond_list.txt')
        with open(list_file, 'r') as f_handle:
            self.selected_bonds = []
            for line in f_handle:
                if line[0] in {';', '#' '-'}:
                    continue
                drname = line[7:11].strip()
                drnum = line[12:19].strip()
                dnum = line[24:31].strip()
                arname = line[47:51].strip()
                arnum = line[52:59].strip()
                anum = line[64:71].strip()
                self.selected_bonds.append((drname + drnum + '_' + arname + arnum, (dnum, anum)))
        group_file = os.path.join(directory, 'group_index.ndx')
        with open(group_file, 'r') as f_handle:
            self.group_index = []
            index = 0
            pair = None
            for line in f_handle:
                line = line.strip()
                if line[0] == '[':
                    pair = line[2:-2]
                else:
                    atoms = tuple(line.split())
                    self.group_index.append((pair, atoms, index))
                    index += 1
        unique_file = os.path.join(directory, 'hbond_list_unique.txt')
        with open(unique_file, 'r') as f_handle:
            self.unique_meta = {}
            self.unique_pair = {}

            for line in f_handle:
                line = line.strip()
                if ':' in line:
                    param, value = tuple(map(lambda x: x.strip(), line.split(':')))
                    self.unique_meta[param] = value
                else:
                    donor, _, acceptor, occupancy = tuple(line[:-1].split())
                    self.unique_pair[donor + '_' + acceptor] = occupancy
        matrix_file = os.path.join(directory, 'hbmap.xpm')
        with open(matrix_file, 'r') as f_handle:
            _, matrix, _, _ = read_xpm(f_handle, reverse=True)
            self.consensus = np.array(matrix).transpose()
        distance_file = os.path.join(directory, 'hbond_dist.xvg')
        self.time = None
        self.distance = None
        with open(distance_file, 'r') as f_handle:
            _, distance_data = read_xvg(f_handle)
            self.distance = np.array(distance_data).transpose()
            self.time = self.distance[:, 0]
            self.distance = self.distance[:, 1:]

        # TODO assert doesn't work for arrays
        # assert(self.selected_bonds and self.group_index and self.unique_meta and self.unique_pair and self.consensus
        #        and self.time and self.distance)

    def calc_missing(self, uniques, pairs):
        # find missing
        for pair in uniques:
            if pair not in self.unique_pair:
                self.missing.add(pair)
        # look for already computed pairs
        missing = set()
        for pair1 in self.missing:
            found = False
            for pair2 in self.group_index:
                if pair2[0] == pair1:
                    self.found.append(pair2)
                    found = True
            if not found:
                missing.add(pair1)

        if missing:
            # compute missing pairs

            # 1) prepare index and selection
            index_file = os.path.join(self.directory, str(uuid.uuid4()) + '.ndx')
            distance_file = os.path.join(self.directory, str(uuid.uuid4()) + '.xvg')
            f_handle = open(index_file, 'w')
            selection_string = ""
            index = 0
            labels = []
            for pair in missing:
                assert(pair in pairs)
                for atoms in pairs[pair]:
                    f_handle.write('[ {} ]\n'.format(pair))
                    f_handle.write("{0:>6} {1:>6}\n".format(atoms[0], atoms[1]))
                    selection_string += '{}\n'.format(index)
                    labels.append(pair)
                    index += 1
            f_handle.close()
            ####################
            #                  #
            # run gmx distance #
            #                  #
            ####################
            command = ['gmx', 'distance', '-f', self.unique_meta['trajectory'], '-s',
                       self.unique_meta['topology'], '-oav', distance_file, '-n', index_file]
            if 'dt' in self.unique_meta:
                command.append('-dt')
                command.append(self.unique_meta['dt'])
            sys.stdout.write('Starting gmx distance\n\nRunning, please wait...\n\n')
            gmx = subprocess.run(command, input=selection_string, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 universal_newlines=True)
            if gmx.returncode:
                delete_if_exists(distance_file)
                sys.stderr.write('\nERROR: called Gromacs with line:\n{}\n\n'.format(' '.join(command)))
                sys.stderr.write('RETURN CODE: {}\n\n'.format(gmx.returncode))
                sys.stderr.write('STDOUT:\n{}\n\n'.format(gmx.stdout))
                sys.stderr.write('STERR:\n{}\n\n'.format(gmx.stderr))
                sys.exit(1)
            # sys.stdout.write('{}\n'.format(gmx.stdout))
            sys.stdout.write('Successfully completed gmx distance\n\n')
            ####################################
            #                                  #
            # read distance file and delete it #
            #                                  #
            ####################################
            with open(distance_file, 'r') as f_handle:
                _, distance_data = read_xvg(f_handle)
            os.remove(distance_file)
            os.remove(index_file)

            # TODO refactor read_xvg to avoid this transpose
            distances = np.array(distance_data).transpose()

            distance = distances[:, 1:]
            consensus = np.zeros(distance.shape, dtype='float32')
            self.write_distances(distance, consensus, labels)

    def calc_found(self):
        ####################################
        #                                  #
        # read distance file and delete it #
        #                                  #
        ####################################
        if self.found:
            labels = []
            columns = []
            for label, _, index in self.found:
                labels.append(label)
                columns.append(index)
            distance = self.distance[:, columns]
            consensus = self.consensus[:, columns]
            self.write_distances(distance, consensus, labels)

    def write_distances(self, distance, consensus, labels):
        ##########################################
        #                                        #
        # collapse hbonds between same residuals #
        #                                        #
        ##########################################
        non_unique_pairs = []  # list of [label, start_column, extent]
        for i, label in enumerate(labels):
            if len(non_unique_pairs) == 0:
                non_unique_pairs.append([label, i, 1])
            elif non_unique_pairs[-1][0] == label:
                non_unique_pairs[-1][2] += 1
            else:
                non_unique_pairs.append([label, i, 1])

        unique_pairs = [non_unique_pairs.pop(0)]  # list of [label, startcolumn1, extent1, startcolumn2, extent2, ...]
        while len(non_unique_pairs) > 0:
            pair1 = non_unique_pairs.pop(0)
            found = False
            for i, pair2 in enumerate(unique_pairs):
                if pair1[0] == pair2[0]:
                    unique_pairs[i].append(pair1[1])
                    unique_pairs[i].append(pair1[2])
                    found = True
                    break
            if not found:
                unique_pairs.append(pair1)

        nrow = distance.shape[0]
        ncol = len(unique_pairs)

        distances1 = np.zeros((nrow, ncol), dtype='float32')
        consensus1 = np.zeros((nrow, ncol), dtype='float32')

        for col, pair in enumerate(unique_pairs):
            n_chunks = (len(pair) - 1) // 2
            for row in range(nrow):
                dist = 1.0e20
                cons = 0.0
                for i in range(n_chunks):
                    start_column = pair[1 + i * 2]
                    extent = pair[2 + i * 2]
                    d1 = distance[row, start_column:start_column + extent].min()
                    dist = min(dist, d1)
                    cons = cons + consensus[row, start_column:start_column + extent].sum()
                distances1[row, col] = dist
                if cons > 0.0:
                    consensus1[row, col] = 1.0
                else:
                    consensus1[row, col] = 0.0

        #######################
        #                     #
        # write all the files #
        #                     #
        #######################
        for col, pair in enumerate(unique_pairs):
            r1, r2 = tuple(pair[0].split('_'))
            bond = r1 + '-' + r2
            f_name = os.path.join(self.directory, 'hbond-' + bond + '.csv')
            backup_if_exists(f_name)
            with open(f_name, 'w') as f_handle:
                f_handle.write('Time,Distance,Consensus\n')
                for row in range(nrow):
                    if consensus1[row, col]:
                        cons = 'Present'
                    else:
                        cons = 'None'
                    f_handle.write('{0:},{1:.3f},{2:}\n'.format(self.time[row], distances1[row, col], cons))
            sys.stdout.write('File {} written\n'.format(f_name))


def main():
    ######################
    #                    #
    # Parse command line #
    #                    #
    ######################
    p = argparse.ArgumentParser()
    p.add_argument('-dirs', required=True, dest='directories', nargs='*',
                   help='directories with different trajectories to analize')
    p.add_argument('-hdir', required=True, dest='h_dir',
                   help='subdirectory where the h bond analysis was performed (should be the same for all directories!')
    p.add_argument('-o', nargs='?', action=check_write_ext({'txt'}), dest='output_file',
                   help='comprehensive hbond report file .txt', default='hbonds.txt')
    args = p.parse_args()

    # TODO mancano tutti i controlli dei parametri

    directories = [os.path.join(os.path.abspath(directory), args.h_dir) for directory in args.directories]
    analyses = [Analysis(directory) for directory in directories]

    all_unique_pairs = set()
    all_pairs = dict()
    for analysis in analyses:
        pairs = set(analysis.unique_pair.keys())
        all_unique_pairs = all_unique_pairs.union(pairs)
        for pair, atoms, _ in analysis.group_index:
            if pair in all_pairs:
                all_pairs[pair].append(atoms)
            else:
                all_pairs[pair] = [atoms]

    for analysis in analyses:
        analysis.calc_missing(all_unique_pairs, all_pairs)
        analysis.calc_found()


if __name__ == '__main__':
    main()