# MitoMut
MitoMut is a tool to call mitochondrial deletions from next generation sequencing (NGS) data 

# Dependencies
To run MitoMut, you need the following dependencies:
  1. python2.7+ or python3.0+
  2. PySam (Python Module that must be downloaded into your version of Python)
  3. Samtools (must be on environment path)
  4. BLAT (must be on environment path)

To check your dependencies, run python check_dependencies.py on your environments command line or terminal.

# Sample Runs
MitoMut Usage:
  python MitoMut.py [-f/-c/-s/-si/-e/-l] *.bam
  
  *Note*– Flags may go in any order. If they require an argument, it must
  come directly after the flag.
 
  -h: This flag shows the help screen with all optional flags and usage information.\n
  -c: This flag specifies the header for the mitochondrial section of the genome.
  -f: This flag specifies the path to fasta reference genome 
  -d: Sets the write directory for all intermediate and result files
  -e: Use this flag if the reference genome is a whole genome reference rather than
      mitochondrial only.
  -q: Sets the minimum quality score to be considered a deletion
  -s: Sets the minimum number of reads required to support a deletion
  -l: Sets the length of the mitochondrial genome (only use if not human)
  -si: sets the minimum number of deleted bases required to be qualified as a deleiton
  
  Flag Defaults:
  
  -c: MT
  -f: mt.fasta (supplied when downloading MitoMut)
  -q: 30
  -s: 5
  -si: 5
  -l: 16569
  -d: .

Sample Runs:
  (Do not type what is in parentheses) 
  
  python MitoMut.py test.bam (A standard run with no extra configuration)
  python MitoMut.py -c chrM test.bam (The bam file's mitochondrial header is chrM instead of MT)
  python MitoMut.py -d /Users/example_user/ test.bam (Writing to a different directory than where the bam file is)
  python MitoMut.py -e -f genome.fasta test.bam (extracting the mitochondrial portion of the reference genome)
  python MitoMut.py -s 20 test.bam (Only deletions with at least 20 supporting reads will be reported)
