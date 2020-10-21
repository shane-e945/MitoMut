# MitoMut
MitoMut identifies mitochondrial deletions from paired-end next generation sequencing (NGS) data.
</br>(Paper: <a href="https://dl.acm.org/doi/10.1145/3307339.3342158">Mitomut</a>)

# Citation
If our work is useful to you, please cite:
```Bash
@inproceedings{inproceedings,
  author = {Elder, C. Shane and Welsh, Catherine},
  year = {2019},
  month = {09},
  pages = {177-182},
  title = {MitoMut: An Efficient Approach to Detecting Mitochondrial DNA Deletions from Paired-end Next-generation Sequencing Data},
  isbn = {978-1-4503-6666-3},
  doi = {10.1145/3307339.3342158}
}
```

# Dependencies
To run MitoMut, you need the following dependencies: </br>
  1. python2.7+ or python3.0+ </br>
  2. PySam (Python Module that must be downloaded into your version of Python) </br>
  3. Samtools (must be on environment path) </br>
  4. BLAT (must be on environment path) (available at http://hgdownload.soe.ucsc.edu/admin/exe/) </br>
</br>
<strong>To check the presence of dependencies, cd into MitoMut's directory and run python check_dependencies.py on your environment's command line or terminal.</strong>

# Instructions
MitoMut Usage:</br>
  python MitoMut.py [-f/-c/-s/-si/-e/-l] *.bam</br></br>
  
  <strong><em>*Note*– Flags may go in any order. If they require an argument, it must</br>
  come directly after the flag.</em></strong></br>
 
  -h: This flag shows the help screen with all optional flags and usage information.</br>
  -c: This flag specifies the header for the mitochondrial section of the genome.</br>
  -f: This flag specifies the path to fasta reference genome </br>
  -d: Sets the write directory for all intermediate and result files</br>
  -e: Use this flag if the reference genome is a whole genome reference rather than</br>
      mitochondrial only.</br>
  -q: Sets the minimum quality score to be considered a deletion</br>
  -s: Sets the minimum number of reads required to support a deletion</br>
  -l: Sets the length of the mitochondrial genome (only use if not human)</br>
  -si: sets the minimum number of deleted bases required to be qualified as a deletion</br>
  </br>
  Flag Defaults:</br>
  </br>
  -c: MT</br>
  -f: mt.fasta (supplied when downloading MitoMut)</br>
  -q: 30</br>
  -s: 5</br>
  -si: 5</br>
  -l: 16569</br>
  -d: . </br>

# Sample Runs
  (Do not type what is in parentheses) </br>
  </br>
  python MitoMut.py test.bam (A standard run with no extra configuration)</br>
  python MitoMut.py -c chrM test.bam (The bam file's mitochondrial header is chrM instead of MT)</br>
  python MitoMut.py -d /Users/example_user/ test.bam (Writing to a different directory than where the bam file is)</br>
  python MitoMut.py -e -f genome.fasta test.bam (extracting the mitochondrial portion of the reference genome)</br>
  python MitoMut.py -s 20 test.bam (Only deletions with at least 20 supporting reads will be reported)</br>
  
# Test Installation (Works on any Unix based operating system)

Once unzipped, test your installation with the following steps: </br>
  1. cd into MitoMut's download directory </br>
  2. Check your dependencies: python ./check_dependencies.py </br>
  3. If all dependencies are installed: python ./MitoMut.py -f mt.fasta -c chrM test.bam </br>
  4. Check the output files: diff test_results.txt test.bam_results.txt </br>
  5. If diff displays nothing, all is successful! </br>
