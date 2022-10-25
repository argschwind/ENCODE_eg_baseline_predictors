## Create is nearest TSS or gene baseline predictor from CRE-TSS or CRE-gene pairs

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
def compute_nearest_pred(line, cell_type):
  
  # split line into different fields 
  entry = line.rstrip().split('\t')
  
  # compute nearest TSS/gene predictor score
  score = entry[9] == entry[15]
  score = str(score * 1)
  
  # create and return output entry
  output = entry[0:4] + ['NA'] + [entry[9]] + ['NA'] + [entry[12]] + [cell_type] + [score] + [entry[13]]
  return output

# execute if run as main program -------------------------------------------------------------------

# parse command line arguments
if __name__ == '__main__':

  # parse command line arguments
  parser = argparse.ArgumentParser(description = 'Compute is nearest TSS or gene predictor from \
    CRE-TSS/gene pairs.')
  parser.add_argument('-i', '--inputfile', help = 'File containing CRE-TSS/gene pairs.')
  parser.add_argument('-c', '--cell_type', help = 'Cell type for CRE-TSS/gene pairs.', default = 'NA')
  args = parser.parse_args()
  
  # check whether input is a (gzipped) file or passed via standard in and process accordingly
  if args.inputfile:
    if is_gzipped(args.inputfile):
      with gzip.open(args.inputfile, 'rt') as f:
        print_header()
        for line in f:
          annot_entry = compute_nearest_pred(line, args.cell_type)
          print('\t'.join(annot_entry))
    else:
      with open(args.inputfile, 'rt') as f:
        print_header()
        for line in f:
          annot_entry = compute_nearest_pred(line, args.cell_type)
          print('\t'.join(annot_entry))
  else:
    print_header()
    for line in sys.stdin:
      annot_entry = compute_nearest_pred(line, args.cell_type)
      print('\t'.join(annot_entry))
