import argparse
import os


def unquote(s):
    return s[1+s.find('"'):s.rfind('"')]


def uncomment(s):
    return s[3+s.find('/*'):s.rfind('*/')-1]


def read_axis(s):
    s1 = s[1+s.find(':'):s.rfind('*/')].split(" ")
    return [i for i in s1]


def col(c):
    color = c.split('/*')
    value = unquote(color[1])
    tmp = unquote(color[0]).split(' c ')
    key = tmp[0].strip()
    color = tmp[1].strip()
    return key, (color, value)


def process_meta(meta):
    uncommented = []
    for i in meta:
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


def read_xpm(name):

    # Open the xpm file for reading
    xpm = open(name)

    # Read in lines until we fidn the start of the array
    meta = [xpm.readline()]
    while not meta[-1].startswith("static char *gromacs_xpm[]"):
        meta.append(xpm.readline())

    # The next line will contain the dimensions of the array
    dim = xpm.readline()
    # There are four integers surrounded by quotes
    nx, ny, nc, nb = [int(i) for i in unquote(dim).split()]

    # The next dim[2] lines contain the color definitions
    # Each pixel is encoded by dim[3] bytes, and a comment
    # at the end of the line contains the corresponding value
    lut = dict([col(xpm.readline()) for _ in range(nc)])

    colors = []
    values = []

    for i in xpm:
        if i.startswith("/*"):
            meta.append(i)
            continue
        j = unquote(i)
        colors.append([lut[j[k:k + nb]][0] for k in range(0, nx, nb)])
        values.append([lut[j[k:k + nb]][1] for k in range(0, nx, nb)])
    meta = process_meta(meta)

    return colors, values, meta, lut


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


def check_read_ext(choices):
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
                elif not os.access(values, os.R_OK):
                    parser.error("file {} is not readable".format(values))
                else:
                    setattr(namespace, self.dest, values)
            else:
                parser.error("file {} doesn't exists".format(values))
    return Act


def main():
    p = argparse.ArgumentParser()
    p.add_argument('-f', required=True, dest='xpm_file', action=check_read_ext({'xpm'}),
                   help='work directory')
    p.add_argument('-o', required=True, dest='csv_file', action=check_write_ext({'csv'}),
                   help='work directory')
    args = p.parse_args()

    colors, values, meta, lut = read_xpm(args.xpm_file)
    a=1



if __name__ == '__main__':
    main()