from __future__ import absolute_import, division, print_function
import fnmatch
import gzip
import os.path
import sys


def clip_fastq(files, clip_length, expected_length=None):
    """
    Writes a GzipFile object to a new FASTQ.GZ file, with reads truncated at the
    5' end such that output is no shorter than (expected_length - clip_length)

    Parameters
    ----------
    files : I/O object, GzipFile object, list
        2-tuple of file objects: the first to be read from, the second to be
        written to.
    clip_length : int
        Maximum number of nucleotides to clip off from the 5' end of a read.
    expected_length : int
        Expected length of reads in nucleotides.
    """
    for infile, outfile in files:
        for i, line in enumerate(infile):
            if i % 2 == 0:
                outfile.write(line)
            else:
                if expected_length is not None:
                    # Do not clip off so much that the read would be shorter
                    # than expected_length - clip_length
                    clip_length -= expected_length - (len(line) - 1)
                outfile.write(line[max(0, clip_length):])


def find_fastq(top):
    """
    Queue all FASTQ files in selected directory.

    Parameters
    ----------
    top : str
        Directory in which to search for FASTQ files.

    Returns
    -------
    filepaths : generator
        Generator that returns a FASTQ file's path, relative to the present
        working directory.
    """
    for path, dirs, files in os.walk(top):
        for name in fnmatch.filter(files, '*.fastq*'):
            yield os.path.join(path, name)


def smart_open(filepaths):
    """
    Produce two file objects: one to be read in as input and another to be
    written to as output.

    Parameters
    ----------
    filepaths : str, list
        Path to file with respect to the present working directory.

    Returns
    -------
    files : generator
        Generator that returns a tuple of I/O objects or GzipFile objects: the
        first to be read from, the second to be written to.
    """
    for filepath in filepaths:
        with gzip.open(filepath, 'rb') as file:
            try:
                file.readline()
            except (IOError, OSError):
                infile = open(filepath, 'rb')
            else:
                infile = gzip.open(filepath, 'rb')
            outfile = gzip.open(filepath.replace('grouped', 'clipped'), 'wb', 6)
            yield infile, outfile
            infile.close(), outfile.close()


if __name__ == '__main__':
    # Script can be run from command line using:
    #   python <directory containing FASTQs> <clip length> <expected length>
    args = sys.argv
    
    # Import clip length
    assert args[2].isdigit(), 'Clip length must be integer'
    clip_length = int(args[2])
    
    # Import expected length
    try:
        assert args[3].isdigit(), 'Expected length must be integer'
    except IndexError:
        expected_length = None
    else:
        expected_length = int(args[3])

    # Run clip algorithm for all FASTQs in directory
    filepaths = find_fastq(args[1])
    files = smart_open(filepaths)
    if not os.path.exists(args[1].replace('grouped', 'clipped')):
        os.makedirs(args[1].replace('grouped', 'clipped'))
    clip_fastq(files, clip_length, expected_length)
        