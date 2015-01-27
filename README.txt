######################################################################
#                                                                    #
#         ConPADE: Contig Ploidy and Allele Dosage Estimation        #
#                                                                    #
######################################################################
#                                                                    #
# Version 1.0-0                                                      #
# Last updated: 01/27/2015                                           #
#                                                                    #
######################################################################


For licensing information, please read file LICENSE.txt distributed with this package


######################################################################


First steps:

To run a sample data set, simply open a command line, navigate to the ConPADE folder and type

\YourWorkingDirectory\RunTestData.bat

or issue the command

\YourWorkingDirectory\ConPADE -bamNames TestData.bam

Make sure all downloaded files are in the folder.


######################################################################


Result files:

Output consists of four files for each contig in each BAM file.
- logLikelihoods: contains log-likelihoods for each evaluated ploidy
- ploidy: a single integer indicating the most likely ploidy
- readStats: read usage statistics
- SNP: a table with identified variants


######################################################################


To get help on usage and detailed information about arguments, open a command line, navigate to the ConPADE folder and type

YourWorkingDirectory\ConPADE


######################################################################
