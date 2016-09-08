import gzip
import os.path


def clip_fastq(fq_in, fq_out, clip_length):
    with smart_open(fq_in, 'rb') as f_in, smart_open(fq_out, 'wb') as f_out:
        for i, line in enumerate(f_in):
            if i % 2 == 0:
                f_out.write(line)
            else:
                f_out.write(line[clip_length:])


def file_is_gzipped(fname):
    with gzip.open(fname, 'rb') as f:
        try:
            f.readline()
        except OSError:
            return False
        else:
            return True


def smart_open(fname, mode):
    # Check file extension if opening in write mode or file does not yet exist
    if mode.startswith('w') or not os.path.isfile(fname):
        if fname.endswith('.gz'):
            return gzip.open(fname, mode, 6)
        else:
            return open(fname, mode)

    # Otherwise check if contents have actually been compressed by gzip
    if file_is_gzipped(fname):
        return gzip.open(fname, mode)
    else:
        return open(fname, mode)


if __name__ == '__main__':
    fq_in = 'raw/an_hytb1_t24_1_zm_S10_L001_R1_001.fastq.gz'
    fq_out = 'cut/out3.fq.gz'
    clip_length = 8
    clip_fastq(fq_in, fq_out, clip_length)
