# adapter-clipper

### Function

Cuts the first <clip_length> nucleotides from each RNA-Seq read, but not such
that the cut read is shorter than <expected_length> - <clip_length> nucleotides
long.

### Terminal command

`python <directory_containing_FASTQs> <clip_length> <expected_length>`
