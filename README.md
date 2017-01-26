<<<<<<< HEAD
SeedSQrrL: A new tool for mitochondrial Genome Reconstruction in Non-Model Organisms

The tool consists of three programs: 1) MitoDBMaker, and 2) MitoDBRelativeMaker, and 3) MitoDBExtractor. They are to be run in that order. The first program creates an SQLite database, populates it with the gene sequence and metadata (from the NCBI nucleotide database) for a list of organisms that you provide then generates a list of unmatched genes for those organisms. If a list of genes is not provided via command line, the default list is used. The next program will find relatives for the resulting list of unmatched genes using the NCBI Taxonomy DB as a guide and then populate the database with the genes for those relatives. The final program extracts those genes into individual files to be run with mitobim. From there, the genes are extended using mitobim and contigs can be placed together (e.g. using Geneious).


How to Install

Use Anaconda Python 2.7 distribution. Install the necessary bs4 and biopython (for mitobim) modules using the following command in a terminal:
conda install beautifulsoup4 biopython
Copy the three SeedSQrrL programs and the automate_mitobim.py program into the folder where you are going to run them.
Install mitobim and its necessary packages (mira, perl if not already installed but probably is). For mitobim, just copy MITObim_1.8.pl to /usr/local/bin or elsewhere in your $PATH and chmod a+x. For mira, you can download the platform-appropriate compiled binaries (mira_4.0.2_darwin…tar.bz2 for OSX, mira_4.0.2_linux-gnu...tar.bz2 for linux), and copy everything from the tar.bz2 file’s bin/ folder into the same place as mitobim above (e.g. /usr/local/bin). Copy everything from the mira library (lib) folder into usr/local/lib
Notes: use ls -l to make sure that the permissions are correct (-rwxr-xr-x). If there is an @ at the end of the permissions, use xattr -d com.apple.quarantine /usr/local/bin/filename on the file. This will likely need to be done on the mira executable and the library files.

How to run

To Run SeedSQrrL:

1. Run mitoDBMaker
python mitoDBmaker.py [samplelist].csv [optional genelist]

Example, let’s say you wanted ND5 in addition to the six default genes:
python mitoDBmaker.py Zaher.csv [\'CO1\',\'ND2\',\'12S\',\'16S\',\'COX1\',\'trnV\',\'ND5\']
2. Run mitoDBRelativeMaker
python mitoDBRelativemaker.py [samplelist]NeedReference.csv

3. Create folder called “seeds”. Run mitoDBExtractor to Generate seed files.






To run Mitobim
Organize your samples in a folder called RawReads, with individual samples underneath, in folders named named Sample_PXXXX_FG_IXXXX. The (still gzipped) sample files go inside. For each sample, if there is a folder called redo inside, samples therein will be used instead:

Make sure all seeds are in a folder called seeds, and are named by the sample id (without the Sample_ part):

Run python automate_mitobim.py







Trouble shooting

Be sure there are no spaces between items in the sample list.

If an organism is not found in the taxonomy database or a gene is not found, check that the spelling is correct.

Currently gene synonyms for six default genes are written into the program using information found in a human gene database. If there are a different set of genes you wish to use, you should add potential synonyms into the mitoDBmaker and mitoDBRelativeMaker programs. I would like to eventually automate this but I have to find a suitable database (e.g. not restricted to human mitochondrial gene synonyms).
Example:
if gene == "ND5":
gene_synonym = ["Mitochondrially Encoded NADH:Ubiquinone Oxidoreductase Core Subunit 5", "NADH Dehydrogenase Subunit 5", "EC 1.6.5.3", "MTND5", "ND5", "Mitochondrially Encoded NADH Dehydrogenase 5","NADH Dehydrogenase", Subunit 5 (Complex I)", "NADH-Ubiquinone Oxidoreductase Chain 5", "Complex I ND5 Subunit", "NADH Dehydrogenase 5", "NADH5"]


Notes on improvements:
Want to run select genes at a time, select the highest rank, and error with mira.
=======
# SeedSQrrL
>>>>>>> 5e98bdc8bb1f1cae1839e904c133dfc0332c0fb2
