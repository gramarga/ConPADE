######################################################################
#                                                                    #
#         ConPADE: Contig Ploidy and Allele Dosage Estimation        #
#                                                                    #
######################################################################
#                                                                    #
# Version 1.0-4                                                      #
# Last updated: 10/30/2015                                           #
#                                                                    #
######################################################################


For licensing information, please read file LICENSE.txt distributed with this package


######################################################################


Installing from source:

- ConPADE is written in C# and can be built with any compatible compiler (Visual Studio, Xamarin Studio or MonoDevelop)
- To compile ConPADE from source, you need to link against the .NET Bio libraries. These can be found at https://github.com/dotnetbio/bio

After successfully building the binaries, make sure you also have the error model files (errorModel.bin and substModel.bin) from the main distribution folder.


######################################################################


First steps:

To run a sample data set, simply open a command line, navigate to the ConPADE folder and type

\YourWorkingDirectory\RunTestData.bat

or issue the command

\YourWorkingDirectory\ConPADE -bamName TestData.bam

Make sure all downloaded files are in the folder.


######################################################################


Result files:

Default behavior is to produce three files from the input BAM file.
- ploidy: one line per contig, with the second column indicating the most likely ploidy, followed by the log-likelihoods for each evaluated ploidy
- readStats: read usage statistics, a table containing information on numbers of aligned reads and base pairs for each contig
- SNP: a table with identified variants, one SNP per line

Optionally, argument -splitContigs can be used to produce four files for each individual contig.
- logLikelihoods: contains log-likelihoods for each evaluated ploidy
- ploidy: a single integer indicating the most likely ploidy
- readStats: read usage statistics
- SNP: a table with identified variants


######################################################################


To get help on usage and detailed information about arguments, open a command line, navigate to the ConPADE folder and type

YourWorkingDirectory\ConPADE


######################################################################


Version changes:

* 1.0-1
- Using updated .Net Bio version 2.0 to fix BAM parsing
- Output is now combined for all contigs (there is an option to split files for individual contigs)
- Changed argument from -bamNames to -bamName

* 1.0-2
- Bug fix: properly ignore not aligned reads

* 1.0-3
- Bug fix: added support for dummy reads in BAM file
- Bug fix: hard-clipped sequences are now correctly parsed

* 1.0-4
- First source code release
- Removed dependency on Escience.dll