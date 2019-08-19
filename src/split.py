#!/usr/bin/python
import sys
import getopt
import os

def read_a_fasta(filehandle):
    """
    read a fasta sequence from the filehandle
    
    usage:
    f = open('myfile.fa', 'r')
    cursor = read_a_fasta(f)
    for record in cursor:
        print record[0]  #fasta head
        print record[1]  #sequence
    """

    head = ''
    seq_line = []

    # skip empty lines, read first fasta head
    for line in filehandle:
        line = line.strip()
        if not line:
            continue
        elif not line.startswith('>'):
            sys.exit('Input file error.\n')
        else:
            head = line
            break

    # read all left contents
    for line in filehandle:
        line = line.strip()
        if not line:
            continue
        elif line.startswith('>'):
            # return one fasta head and its seq
            yield [head, ''.join(seq_line)]
            head = line
            seq_line = []
        else:
            seq_line.append(line)
    
    # return last fasta head and its
    yield [head, ''.join(seq_line)]

def write_a_fasta(filehandle, head, seq):
    """
    print a fasta sequence to the filehandle

    usage:
    f = open('out.fa', 'w')
    write_a_fasta(f, head, seq)

    or:
    t = (head, seq)
    write_a_fasta(f, *t)
    """
    filehandle.write('%s\n%s\n' % (head, seq))

def main(input_file, output_path, part_number):

    # generate output file handls
    out_files = []
    for i in range(1, part_number + 1):
        print(i)
        out_files.append(open(os.path.join(output_path, 'part%d' % i), 'w'))

    # split
    n = 0
    f = open(input_file)
    cursor = read_a_fasta(f)
    for fasta in cursor:
        f_current = out_files[n % part_number]
        write_a_fasta(f_current, fasta[0], fasta[1])
        n += 1

    f.close()
    for i in out_files:
        i.close()


def usage():
    print("""
Usage:

Essentail:
    -i, --input-file=<string>       Input fasta file. 
    -p, --part-number=<int>         Part number to split. (> 1)

Optional:
    -o, --output-path=<srting>      Directory to output. The directory should be exsits.
    -h, --help                      Show this information.
""")
    

if __name__ == '__main__':
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:o:p:h', ['input_file=','part-number=','output-path=','help'])
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    input_file  = ''
    output_path = os.path.abspath('.')
    part_number = 0

    flag = 0
    for o, a in opts:
        if o in ['-i', '--input_file']:
            input_file = os.path.abspath(a)
            flag += 1
        elif o in ['-o', '--output-path']:
            output_path = os.path.abspath(a)
        elif o in ['-p', '--part-number']:
            part_number = int(a)
            flag += 1
        elif o in ['-h', '--help']:
            usage()
            sys.exit(0)
        else:
            usage()
            sys.exit(1)

    if part_number <= 1 or flag < 2:
        usage()
        sys.exit(1)

    if not os.path.exists(output_path):
        os.makedirs(output_path)
    main(input_file, output_path, part_number)