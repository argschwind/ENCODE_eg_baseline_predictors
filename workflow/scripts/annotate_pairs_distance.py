## process CRE-TSS/gene pairs line-by-line and compute distance to TSS and distance to feature (i.e.
## distance to gene body for cre-gene pairs)

# import modules
import sys, argparse, gzip

# define functions ---------------------------------------------------------------------------------

# function to check if input file is gzipped
def is_gzipped(path):
  with open(path, 'rb') as f:
    return f.read(2) == b'\x1f\x8b'

# function to process one line in the input and compute distance to TSS and feature for that pair
def annotate_pair(line, start_tss, end_tss, chr_tss):
  
  # split line into different fields 
  entry = line.rstrip().split('\t')
  
  # get tss coordinates
  if entry[chr_tss] == "+":
    tss_bed = int(entry[start_tss])
    tss = tss_bed + 1
  else:
    tss_bed = int(entry[end_tss])
    tss = tss_bed
    
  # compute distance to tss
  cre_center = (int(entry[1]) + int(entry[2])) / 2
  dist_to_tss = round(abs(cre_center - tss))
  
  # compute distance to feature (e.g. gene body)
  feature = [int(entry[start_tss]) + 1, int(entry[end_tss])]
  if (cre_center >= feature[0]) and (cre_center <= feature[1]):
    dist_to_feature = 0
  else:
    dist_to_feature = round(min([abs(x - cre_center) for x in feature]))
  
  # create and return output entry  
  out_entry = entry + [str(tss_bed), str(dist_to_tss), str(dist_to_feature)]
  return out_entry

# execute if run as main program -------------------------------------------------------------------

# parse command line arguments
if __name__ == '__main__':

  # parse command line arguments
  parser = argparse.ArgumentParser(description = 'Annotate CRE-TSS/pairs with distance.')
  parser.add_argument('-i', '--inputfile', help = 'File containing CRE-TSS/pairs')
  parser.add_argument('-c', '--crecols', help = 'Number of CRE columns in the input file. Used to \
    separate CRE and gene columns.', type = int, default = 6)
  args = parser.parse_args()
  
  # get TSS coordinate fields based on gene column argument
  start_tss = args.crecols + 1
  end_tss = args.crecols + 2
  chr_tss = args.crecols + 5
  
  # check whether input is a (gzipped) file or passed via standard in and process accordingly
  if args.inputfile:
    if is_gzipped(args.inputfile):
      with gzip.open(args.inputfile, 'rt') as f:
        for line in f:
          annot_entry = annotate_pair(line, start_tss, end_tss, chr_tss)
          print('\t'.join(annot_entry))
    else:
      with open(args.inputfile, 'rt') as f:
        for line in f:
          annot_entry = annotate_pair(line, start_tss, end_tss, chr_tss)
          print('\t'.join(annot_entry))
  else:
    for line in sys.stdin:
      annot_entry = annotate_pair(line, start_tss, end_tss, chr_tss)
      print('\t'.join(annot_entry))
