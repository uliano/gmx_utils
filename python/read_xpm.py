#!/opt/anaconda/python

from PIL import Image
import numpy as np
import matplotlib.pyplot as plt


def make_image(colors):
    nrows = len(colors)
    ncols = len(colors[0])
    rgb_array = np.zeros((nrows, ncols, 3), 'uint8')
    for r in range(nrows):
        for c in range(ncols):
            value = colors[r][c]
            rgb_array[r, c, 0] = int(value[1:3], 16)
            rgb_array[r, c, 1] = int(value[3:5], 16)
            rgb_array[r, c, 2] = int(value[5:7], 16)
    img = Image.fromarray(rgb_array)
    return img


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

c, v, m, l = read_xpm('ss.xpm')

for i in m:
    print('{}: {}'.format(i,m[i]))

i = make_image(c)

x_axis = {str(i):str(float(i)/1000) for i in m['x-axis'].split(' ')}
print(x_axis)

plt.imshow(i)
for t in plt.xticks()[0]:
    print(t)

plt.xticks(plt.xticks()[0], [x_axis[str(int(t))] for t in plt.xticks()[0] if str(int(t)) in x_axis])

plt.savefig('ss.pdf', dpi=2400)


