## Create distance to TSS or gene baseline predictor from CRE-TSS or CRE-gene pairs

# import modules
import sys, argparse, gzip

# define functions ---------------------------------------------------------------------------------

# function to check if input file is gzipped
def is_gzipped(path):
  with open(path, 'rb') as f:
    return f.read(2) == b'\x1f\x8b'
  
# helper function to print header
def print_header():
  colnames = ['chr', 'start', 'end', 'name', 'class', 'TargetGene', 'TargetGeneEnsemblID', 
    'TargetGeneTSS', 'CellType', 'Score', 'DistanceToTSS']
  print('\t'.join(colnames))

# function to process one line in the input and compute distance to TSS and feature for that pair
def get_distance(line, cell_type, dist_index):
  
  # split line into different fields 
  entry = line.rstrip().split('\t')
  
  # get distance to TSS/gene predictor score
  score = entry[dist_index]
  
  # create and return output entry
  output = entry[0:4] + ['NA'] + [entry[9]] + ['NA'] + [entry[12]] + [cell_type] + [score] + [entry[13]]
  return output

# execute if run as main program -------------------------------------------------------------------

# parse command line arguments
if __name__ == '__main__':

  # parse command line arguments
  parser = argparse.ArgumentParser(description = 'Compute distance to TSS/gene predictor from \
    CRE-TSS/gene pairs.')
  parser.add_argument('-i', '--inputfile', help = 'File containing CRE-TSS/gene pairs.')
  parser.add_argument('-c', '--cell_type', help = 'Cell type for CRE-TSS/gene pairs.', default = 'NA')
  parser.add_argument('-t', '--type', help = 'Should predictor be computed based on distance to TSS or gene?', default = 'tss')
  args = parser.parse_args()
  
  # process type argument
  if args.type == 'tss':
    dist_index = 13
  elif args.type == 'gene':
    dist_index = 14
  else:
    sys.exit('Invalid type argument.')
  
  # check whether input is a (gzipped) file or passed via standard in and process accordingly
  if args.inputfile:
    if is_gzipped(args.inputfile):
      with gzip.open(args.inputfile, 'rt') as f:
        print_header()
        for line in f:
          annot_entry = get_distance(line, args.cell_type, dist_index)
          print('\t'.join(annot_entry))
    else:
      with open(args.inputfile, 'rt') as f:
        print_header()
        for line in f:
          annot_entry = get_distance(line, args.cell_type, dist_index)
          print('\t'.join(annot_entry))
  else:
    print_header()
    for line in sys.stdin:
      annot_entry = get_distance(line, args.cell_type, dist_index)
      print('\t'.join(annot_entry))
