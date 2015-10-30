using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Diagnostics;
using Bio.Util;
using Bio.Util.Distribute;
using Bio.Util.ArgumentParser;
using Bio.IO.SAM;
using Bio.IO.BAM;


namespace ConPADE
{
    /// <summary>
    /// Class to store the two alleles, read depths, best allele
    /// dosage and posterior probability for each putative SNP.
    /// </summary>
    class Best_Dose
    {
        public long position { get; private set; }
        public byte nuc_one { get; private set; }
        public int count_one { get; private set; }
        public byte nuc_two { get; private set; }
        public int count_two { get; private set; }
        public byte[] best_dose { get; private set; }
        public double[] SNP_posterior { get; private set; }

        /// <summary>
        /// Standard constructor for class Best_Dose.
        /// </summary>
        /// <param name="position">Zero-based position in the reference sequence.</param>
        /// <param name="nuc_one">Reference allele.</param>
        /// <param name="count_one">Reference allele depth.</param>
        /// <param name="nuc_two">Alternative allele.</param>
        /// <param name="count_two">Alternative allele depth.</param>
        /// <param name="best_dose">Most likely allele dosage for each possible ploidy.</param>
        /// <param name="SNP_posterior">Posterior probability of there being a SNP for each possible ploidy.</param>
        public Best_Dose(long position, byte nuc_one, int count_one, byte nuc_two, int count_two,
            byte[] best_dose, double[] SNP_posterior)
        {
            this.position = position;
            this.nuc_one = nuc_one;
            this.count_one = count_one;
            this.nuc_two = nuc_two;
            this.count_two = count_two;
            this.best_dose = best_dose;
            this.SNP_posterior = SNP_posterior;
        }
    }

    /// <summary>
    /// Main class of ConPADE.
    /// </summary>
    class ConPADE : SelfDistributable
    {
        private string _jobName = "ConPADE";
        public override string JobName
        {
            get { return _jobName; }
        }

        /// <summary>
        /// Path to a sorted BAM file.
        /// </summary>
        [Parse(ParseAction.Required, typeof(InputFile))]
        public InputFile bamName = null;

        /// <summary>
        /// Path to error model file.
        /// </summary>
        [Parse(ParseAction.Optional, typeof(InputFile))]
        public InputFile modelFile = new InputFile();
        
        /// <summary>
        /// Path to substitution model file.
        /// </summary>
        [Parse(ParseAction.Optional, typeof(InputFile))]
        public InputFile substFile = new InputFile();

        /// <summary>
        /// Maximum ploidy to evaluate.
        /// </summary>
        [Parse(ParseAction.Optional, typeof(int))]
        public int max_ploidy = 4;

        /// <summary>
        /// Phred-like threshold for outputting a SNP.
        /// </summary>
        [Parse(ParseAction.Optional, typeof(int))]
        public int SNPthres = 40;

        /// <summary>
        /// SNP density.
        /// </summary>
        [Parse(ParseAction.Optional, typeof(int))]
        public int snpDens = 200;

        /// <summary>
        /// Store results for different contigs in separate files.
        /// </summary>
        [Parse(ParseAction.Optional, typeof(bool))]
        public bool splitContigs = false;

        public override void RunTasks(RangeCollection tasksToRun, long taskCount)
        {
            RunFile(bamName.FullName);
        }

        public override void Cleanup(long taskCount)
        {

        }


        // Calculate log(x + y), given log(x) and log(y)
        private static double LogSum(double log_x, double log_y)
        {
            double result;

            if (log_x > (log_y + 40))
            {
                result = log_x;
            }
            else if (log_y > (log_x + 40))
            {
                result = log_y;
            }
            else
            {
                result = Math.Log(Math.Exp(log_x - log_y) + 1) + log_y;
            }

            return result;
        }


        // For each ploidy to be evaluated, determine all possible nucleotide proportions/probabilities.
        // Store the logarithm of each probability.
        private double[][][] Nuc_Props(int min_ploidy, int number_of_ploidies)
        {
            double[][][] nuc_props = new double[number_of_ploidies][][];
            for (int i = 0; i < number_of_ploidies; i++)
            {
                int ploidy = i + min_ploidy;
                nuc_props[i] = new double[ploidy + 1][];

                for (int j = 0; j <= ploidy; j++)
                {
                    double prob = (double)j / ploidy;
                    nuc_props[i][j] = new double[2]
                    {
                        Math.Log(prob),
                        Math.Log(1-prob)
                    };
                }
            }
            return nuc_props;
        }


        // Set the model for the probability of each allele dosage.
        // We currently use uniform probabilities for each heterozygous genotype.
        private double[][] Dose_Probs(int min_ploidy, int number_of_ploidies, double SNP_density, double no_SNP_prob)
        {
            double[][] dose_probs = new double[number_of_ploidies][];

            // For a ploidy of 1: no SNP allowed.
            dose_probs[0] = new double[2];
            dose_probs[0][0] = dose_probs[0][1] = Math.Log(0.5);
            
            for (int i = 1; i < number_of_ploidies; i++)
            {
                int ploidy = i + min_ploidy;
                dose_probs[i] = new double[ploidy + 1];

                dose_probs[i][0] = no_SNP_prob;
                dose_probs[i][ploidy] = no_SNP_prob;

                double other_probs = Math.Log(SNP_density / (ploidy - 1));
                for (int j = 1; j < ploidy; j++)
                {
                    dose_probs[i][j] = other_probs;
                }
            }
            return dose_probs;
        }

        // Get sequencing error probabilities from the error model file.
        // Error model file contains double values nested according to the following sequence:
        // GG precedes - quality score - neighboring quality score - true nucleotide - error/no error
        private double[, , , ,] Error_Probs()
        {
            if (modelFile.FullName == null)
            {
                modelFile.FullName = "errorModel.bin";
            }
            if (!File.Exists(modelFile.ToString()))
            {
                throw new FileNotFoundException(String.Format("File {0} not found.", modelFile.FullName.ToString()));
            }

            double[, , , ,] log_probs = new double[2, 40, 40, 4, 2];
            bool got_file_handle = false;
            while (!got_file_handle)
            {
                got_file_handle = true;
                try
                {
                    using (FileStream model_stream = new FileStream(modelFile.ToString(), FileMode.Open))
                    {
                        using (BinaryReader model_reader = new BinaryReader(model_stream))
                        {
                            for (int is_GG = 0; is_GG < 2; is_GG++)
                            {
                                for (int qual = 0; qual < 40; qual++)
                                {
                                    for (int neigh_qual = 0; neigh_qual < 40; neigh_qual++)
                                    {
                                        for (int nuc = 0; nuc < 4; nuc++)
                                        {
                                            for (int is_correct = 0; is_correct < 2; is_correct++)
                                            {
                                                log_probs[is_GG, qual, neigh_qual, nuc, is_correct] =
                                                    model_reader.ReadDouble();
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                catch (IOException)
                {
                    got_file_handle = false;
                }
            }

            return log_probs;
        }

        // Get nucleotide substitution probabilities from the substitution model file.
        // Substitution model file contains double values nested according to the following sequence:
        // true nucleotide - GG precedes - observed nucleotide
        private double[, ,] Subst_Probs()
        {
            if (substFile.FullName == null)
            {
                substFile.FullName = "substModel.bin";
            }
            if (!File.Exists(substFile.ToString()))
            {
                throw new FileNotFoundException(String.Format("File {0} not found.", substFile.FullName.ToString()));
            }

            double[, ,] log_subst_probs = new double[4, 2, 4];
            bool got_file_handle = false;
            while (!got_file_handle)
            {
                got_file_handle = true;
                try
                {
                    using (FileStream subst_stream = new FileStream(substFile.ToString(), FileMode.Open))
                    {
                        using (BinaryReader subst_reader = new BinaryReader(subst_stream))
                        {
                            for (int real_nuc = 0; real_nuc < 4; real_nuc++)
                            {
                                for (int is_GG = 0; is_GG < 2; is_GG++)
                                {
                                    for (int obs_nuc = 0; obs_nuc < 4; obs_nuc++)
                                    {
                                        log_subst_probs[real_nuc, is_GG, obs_nuc] =
                                            subst_reader.ReadDouble();
                                    }
                                }
                            }
                        }
                    }
                }
                catch (IOException)
                {
                    got_file_handle = false;
                }
            }
            return log_subst_probs;
        }

        // Search the BAM file for the next valid read aligned against the current contig.
        // Update read/base pairs statistics.
        private void Search_Reads(BAMParser parser, ref SAMAlignedSequence next_alignment, string contig_name,
            ref long number_of_aligned_reads, ref long number_of_aligned_base_pairs, ref long number_of_used_reads,
            ref long number_of_used_base_pairs, Queue<Padded_Read> read_queue, long current_position)
        {
            while (next_alignment != null && 
                !next_alignment.IsDummyRead &&
                next_alignment.RName == contig_name &&
                (next_alignment.Pos - 1) == current_position)
            {
                // The next alignment overlaps with current position, so continue.
                number_of_aligned_reads++;
                number_of_aligned_base_pairs += next_alignment.QuerySequence.Count;

                // Maybe we should let the mininum alignment quality be a parameter.
                // We currently leave it for the user to pre-filter the BAM file.
                if (next_alignment.MapQ > 0)
                {
                    number_of_used_reads++;
                    number_of_used_base_pairs += next_alignment.QuerySequence.Count;
                    read_queue.Enqueue(new Padded_Read(next_alignment));
                }
                
                #region Parse BAM file until next alignment is found
                if (!parser.IsEOF())
                {
                    next_alignment = parser.GetAlignedSequence(true);
                    
                    while ((next_alignment == null || next_alignment.RName == "*" || next_alignment.IsDummyRead) && !parser.IsEOF())
                    {
                        next_alignment = parser.GetAlignedSequence(true);
                    }
                }
                else
                {
                    next_alignment = null;
                }
                #endregion Parse BAM file until next alignment is found
            }
        }

        // For each read in the queue that overlaps the current position,
        // extract information on the nucleotide and quality score.
        private void Extract_Read_Info(Queue<Padded_Read> read_queue, long current_position, out byte[] obs_nucs,
            out byte[] is_GG, out bool[] reverse, out int[] quality_scores, out int[] neigh_quality_scores,
            out int[] scores, out int[] counts, out int k)
        {
            obs_nucs = new byte[read_queue.Count];
            is_GG = new byte[read_queue.Count];
            reverse = new bool[read_queue.Count];
            quality_scores = new int[read_queue.Count];
            neigh_quality_scores = new int[read_queue.Count];
            scores = new int[4];
            counts = new int[4];
            k = 0;

            foreach (Padded_Read curRead in read_queue)
            {
                // The check for overlap below may not be necessary depending on how we clear read cache.
                if ((curRead.alignment.Pos + curRead.alignment_length - 2) >= current_position)
                {
                    byte nuc = curRead.padded_sequence[curRead.cur_pos_ind];
                    if (nuc < 4)
                    {
                        obs_nucs[k] = nuc;
                        counts[nuc]++;
                        scores[nuc] += quality_scores[k] =
                            curRead.padded_quality_scores[curRead.cur_pos_ind];
                        
                        // Quality scores are allowed in the [2,41] range.
                        if (quality_scores[k] < 2 )
                        {
                            quality_scores[k] = 3;
                        }
                        else if (quality_scores[k] > 41)
                        {
                            quality_scores[k] = 41;
                        }

                        #region Get neighboring quality scores
                        int max_ind = curRead.located_sequence.Length - 1;
                        int left_found = 0;
                        int right_found = 0;
                        int total_qual = 0;
                        int no_qual = -10;

                        // We currently define the neighboring nucleotides as 5 on each side.
                        int running_ind = curRead.cur_pos_ind - 1;
                        while (left_found < 5 && running_ind > 0)
                        {
                            int temp_qual = curRead.padded_quality_scores[running_ind--];
                            if (temp_qual != no_qual)
                            {
                                left_found++;
                                total_qual += temp_qual;
                            }
                        }

                        running_ind = curRead.cur_pos_ind + 1;
                        while (right_found < 5 && running_ind < max_ind)
                        {
                            int temp_qual = curRead.padded_quality_scores[running_ind++];
                            if (temp_qual != -10)
                            {
                                right_found++;
                                total_qual += temp_qual;
                            }
                        }

                        neigh_quality_scores[k] = (int)Math.Round((double)total_qual / (left_found + right_found));
                        if (neigh_quality_scores[k] < 2)
                        {
                            neigh_quality_scores[k] = 3;
                        }
                        else if (neigh_quality_scores[k] > 41)
                        {
                            neigh_quality_scores[k] = 41;
                        }
                        #endregion Get neighboring quality scores

                        #region Get preceding nucleotides and update values if read is reversed
                        reverse[k] = curRead.is_reverse;
                        int nucs_found = 0;
                        if (reverse[k])
                        {
                            running_ind = curRead.cur_pos_ind + 1;
                            while (nucs_found < 2 && running_ind < max_ind)
                            {
                                byte temp_nuc = curRead.padded_sequence[running_ind++];
                                // If read is reversed, we look for C nucleotides after the current nucleotide.
                                if (temp_nuc == 1)
                                {
                                    nucs_found++;
                                }
                                else
                                {
                                    if (temp_nuc <= 4)
                                    {
                                        break;
                                    }
                                }
                            }
                        }
                        else
                        {
                            running_ind = curRead.cur_pos_ind - 1;
                            while (nucs_found < 2 && running_ind > 0)
                            {
                                byte temp_nuc = curRead.padded_sequence[running_ind--];
                                // If read is not reversed, we look for G nucleotides before the current nucleotide.
                                if (temp_nuc == 2)
                                {
                                    nucs_found++;
                                }
                                else
                                {
                                    if (temp_nuc <= 4)
                                    {
                                        break;
                                    }
                                }
                            }
                        }
                        if (nucs_found == 2)
                        {
                            is_GG[k] = 1;
                        }
                        #endregion Get preceding nucleotides and update values if read is reversed

                        k++;
                    }

                    // Update cur_pos_ind in the padded read to ignore insertions.
                    while (curRead.located_sequence[++curRead.cur_pos_ind] == curRead.located_sequence[curRead.cur_pos_ind - 1]) ;
                }
            }
        }

        // Find two most abundant nucleotides for a given position.
        // We currently use A as a neutral nucleotide, based on its substitution model.
        private static void Get_Two_Nucs(int[] scores, out byte nuc_one, out byte nuc_two)
        {
            if (scores[1] > scores[0])
            {
                nuc_one = 1;
                nuc_two = 0;
            }
            else
            {
                nuc_one = 0;
                nuc_two = 1;
            }
            for (byte i = 2; i < 4; i++)
            {
                if (scores[i] > scores[nuc_one])
                {
                    nuc_two = nuc_one;
                    nuc_one = i;
                }
                else if (scores[i] > scores[nuc_two])
                {
                    nuc_two = i;
                }
            }
        }

        // For a set of aligned nucleotides and their auxiliary information,
        // return the probability of the observed values.
        // This includes both the sequencing error and the substitution models.
        private double[][] Obs_Probs(double[, , , ,] log_probs, double[, ,] log_subst_probs, byte[] obs_nucs,
            byte[] is_GG, bool[] reverse, int[] quality_scores, int[] neigh_quality_scores, int[] counts, int k,
            byte nuc_one, byte nuc_two)
        {
            byte nuc_one_reverse = (byte)(3 - nuc_one);
            byte nuc_two_reverse = (byte)(3 - nuc_two);

            int reads_to_use = counts[nuc_one] + counts[nuc_two];

            double[][] log_nuc_probs = new double[reads_to_use][];
            int l = 0;
            for (int i = 0; i < k; i++)
            {
                if (obs_nucs[i] == nuc_one)
                {
                    byte nuc_one_to_use = nuc_one;
                    byte nuc_two_to_use = nuc_two;
                    if (reverse[i])
                    {
                        nuc_one_to_use = nuc_one_reverse;
                        nuc_two_to_use = nuc_two_reverse;
                    }
                    log_nuc_probs[l++] = new double[2] {
                                    log_probs[is_GG[i], quality_scores[i]-2, neigh_quality_scores[i]-2, nuc_one_to_use, 1],
                                    log_probs[is_GG[i], quality_scores[i]-2, neigh_quality_scores[i]-2, nuc_two_to_use, 0] +
                                    log_subst_probs[nuc_two_to_use, is_GG[i], nuc_one_to_use]
                                };
                }
                else if (obs_nucs[i] == nuc_two)
                {
                    byte nuc_one_to_use = nuc_one;
                    byte nuc_two_to_use = nuc_two;
                    if (reverse[i])
                    {
                        nuc_one_to_use = nuc_one_reverse;
                        nuc_two_to_use = nuc_two_reverse;
                    }
                    log_nuc_probs[l++] = new double[2]  {
                                    log_probs[is_GG[i], quality_scores[i]-2, neigh_quality_scores[i]-2, nuc_one_to_use, 0] + 
                                    log_subst_probs[nuc_one_to_use, is_GG[i], nuc_two_to_use],
                                    log_probs[is_GG[i], quality_scores[i]-2, neigh_quality_scores[i]-2, nuc_two_to_use, 1]
                                };
                }
            }
            return log_nuc_probs;
        }

        // For each ploidy to be evaluated, calculate the likelihood of the genotypes given the observed data.
        private static double[][] Log_Likelihoods(int min_ploidy, int max_ploidy, double[][] log_nuc_probs,
            double[][][] nuc_props)
        {
            int number_of_ploidies = max_ploidy - min_ploidy + 1;

            double[][] log_likelihoods = new double[number_of_ploidies][];

            // Likelihoods of the homozygous genotypes are the same for all ploidies.
            double log_like_P0 = 0;
            double log_like_P1 = 0;
            for (int read = 0; read < log_nuc_probs.Length; read++)
            {
                log_like_P0 += LogSum(
                    log_nuc_probs[read][0] + nuc_props[0][0][0],
                    log_nuc_probs[read][1] + nuc_props[0][0][1]);

                log_like_P1 += LogSum(
                    log_nuc_probs[read][0] + nuc_props[0][min_ploidy][0],
                    log_nuc_probs[read][1] + nuc_props[0][min_ploidy][1]);
            }

            for (int i = 0; i < number_of_ploidies; i++)
            {
                int ploidy = i + min_ploidy;
                log_likelihoods[i] = new double[ploidy + 1];

                log_likelihoods[i][0] = log_like_P0;
                log_likelihoods[i][ploidy] = log_like_P1;

                // Calculate the likelihood of heterozygous genotypes.
                for (int j = 1; j < ploidy; j++)
                {
                    for (int read = 0; read < log_nuc_probs.Length; read++)
                    {
                        log_likelihoods[i][j] += LogSum(
                            log_nuc_probs[read][0] + nuc_props[i][j][0],
                            log_nuc_probs[read][1] + nuc_props[i][j][1]);
                    }
                }
            }

            return log_likelihoods;
        }

        // Calculate the likelihood of each ploidy for the current position, given genotype likelihoods.
        // Use this value to update the global likelihood of each ploidy.
        // Also, find the most likely allele dosage and store the corresponding posterior
        // probability of there being a SNP in the current position, for each ploidy.
        private void Global_Likelihood_Keep_Dose(int min_ploidy, int number_of_ploidies, double[][] dose_probs,
            double[] global_log_like, Queue<Best_Dose> dose_queue, long current_position, int[] counts,
            byte nuc_one, byte nuc_two, double[][] log_likelihoods)
        {
            byte[] best_dose = new byte[number_of_ploidies];
            double[] SNP_posterior = new double[number_of_ploidies];

            for (int i = 0; i < number_of_ploidies; i++)
            {
                int ploidy = i + min_ploidy;

                best_dose[i] = 0;
                double best_posterior = double.NegativeInfinity;
                double aggregate = double.NegativeInfinity;
                for (byte j = 0; j <= ploidy; j++)
                {
                    double posterior = dose_probs[i][j] + log_likelihoods[i][j];
                    aggregate = LogSum(
                        aggregate, posterior);
                    if (posterior > best_posterior)
                    {
                        best_dose[i] = j;
                        best_posterior = posterior;
                    }
                }
                global_log_like[i] += aggregate;

                SNP_posterior[i] = LogSum(
                       dose_probs[i][0] + log_likelihoods[i][0],
                       dose_probs[i][ploidy] + log_likelihoods[i][ploidy]) - aggregate;
            }

            // Store information that will be part of the called SNPs output.
            dose_queue.Enqueue(new Best_Dose(current_position, nuc_one, counts[nuc_one], nuc_two,
                counts[nuc_two], best_dose, SNP_posterior));
        }

        /// <summary>
        /// Run ConPADE on each contig of the input BAM file.
        /// </summary>
        /// <param name="bamName">Name of the input BAM file.</param>
        public void RunFile(string bamName)
        {
            // Current implementation requires that minimum ploidy be 1
            int min_ploidy = 1;
            int number_of_ploidies = max_ploidy - min_ploidy + 1;

            // Set nucleotide proportions (genotypes)
            double[][][] nuc_props = Nuc_Props(min_ploidy, number_of_ploidies);

            // Set dosage probabilities
            double SNP_density = (double)1 / snpDens;
            double no_SNP_prob = Math.Log((1 - SNP_density) / 2);
            double[][] dose_probs = Dose_Probs(min_ploidy, number_of_ploidies, SNP_density, no_SNP_prob);

            // Set HiSeq error model
            double[, , , ,] log_probs = Error_Probs();

            // Set substitution model
            double[, ,] log_subst_probs = Subst_Probs();

            // Set SNP calling probability
            double log_SNP_thres = SNPthres * Math.Log(10) / -10;

            Stopwatch clock = new Stopwatch();

            Console.WriteLine("Program started at {0}\n", DateTime.Now);

            Stream bam_stream = new FileStream(bamName, FileMode.Open, FileAccess.Read);
            BAMParser parser = new BAMParser();
            SAMAlignmentHeader header = parser.GetHeader(bam_stream);
            string temp = Path.GetFileNameWithoutExtension(bamName);

            // Find first valid alignment in BAM file
            SAMAlignedSequence next_alignment = parser.GetAlignedSequence(true);
            while (next_alignment == null || next_alignment.RName == "*" || next_alignment.IsDummyRead)
            {
                next_alignment = parser.GetAlignedSequence(true);
            }

            TextWriter writer_log_like = null;
            TextWriter writer_SNP = null;
            TextWriter writer_ploidy = null;
            TextWriter writer_reads = null;
            
            // Create global output files and write headers.
            if (!splitContigs)
            {
                string SNP_file = temp + "_SNP.txt";
                writer_SNP = new StreamWriter(SNP_file);
                writer_SNP.WriteLine("Contig\tPosition\tAlleles\tCounts\tDosage\tPhredQuality");

                string ploidy_file = temp + "_ploidy.txt";
                writer_ploidy = new StreamWriter(ploidy_file);
                writer_ploidy.Write("Contig\tBestPloidy");
                for (int i = 0; i < number_of_ploidies; i++)
                {
                    writer_ploidy.Write("\tlogLike_M{0}", i + min_ploidy);
                }
                writer_ploidy.WriteLine("");

                string reads_file = temp + "_readStats.txt";
                writer_reads = new StreamWriter(reads_file);
                writer_reads.WriteLine("Contig\tAlignedReads\tAlignedBases\tUsedReads\tUsedBases");
            }

            // Run over each contig in input BAM file.
            int contig_ind = -1;
            while (next_alignment != null && next_alignment.RName != "*" && !next_alignment.IsDummyRead)
            {
                string contig_name = next_alignment.RName;

                Console.WriteLine("Started contig {0} at {1}",
                    contig_name, DateTime.Now);

                clock.Restart();

                #region Variables and file handles for current contig
                long number_of_aligned_reads = 0;
                long number_of_aligned_base_pairs = 0;
                long number_of_used_reads = 0;
                long number_of_used_base_pairs = 0;

                // Create individual output files for the current contig.
                if (splitContigs)
                {
                    string name = temp + "_" + contig_name;
                    string log_like_file = name + "_log_likelihoods.txt";
                    writer_log_like = new StreamWriter(log_like_file);

                    string SNP_file = name + "_SNP.txt";
                    writer_SNP = new StreamWriter(SNP_file);

                    string ploidy_file = name + "_ploidy.txt";
                    writer_ploidy = new StreamWriter(ploidy_file);

                    string reads_file = name + "_readStats.txt";
                    writer_reads = new StreamWriter(reads_file);
                }

                double[] global_log_like = new double[number_of_ploidies];

                while (header.ReferenceSequences[++contig_ind].Name != contig_name) ;
                long contig_length = header.ReferenceSequences[contig_ind].Length;

                // Create a queue to include all reads that overlap with a given position.
                Queue<Padded_Read> read_queue = new Queue<Padded_Read>();

                // Create a queue to include best doses for each tested position.
                Queue<Best_Dose> dose_queue = new Queue<Best_Dose>((int)contig_length);
                #endregion Variables and file handles  for current contig

                int positions_to_compute = 0;
                long current_position = 0;

                #region Run over every position in contig
                while (current_position < contig_length)
                {
                    if ((current_position % 1000000) == 0 && current_position != 0)
                    {
                        Console.WriteLine("At position {0} of {1}", current_position + 1, contig_length);
                    }

                    // Search for reads starting at current position.
                    Search_Reads(parser, ref next_alignment, contig_name, ref number_of_aligned_reads,
                        ref number_of_aligned_base_pairs, ref number_of_used_reads, ref number_of_used_base_pairs,
                        read_queue, current_position);

                    if (read_queue.Count > 0)
                    {
                        positions_to_compute++;

                        // Extract information from each read in queue.
                        byte[] obs_nucs;
                        byte[] is_GG;
                        bool[] reverse;
                        int[] quality_scores;
                        int[] neigh_quality_scores;
                        int[] scores;
                        int[] counts;
                        int k;

                        Extract_Read_Info(read_queue, current_position, out obs_nucs, out is_GG, out reverse,
                            out quality_scores, out neigh_quality_scores, out scores, out counts, out k);

                        // Find two most abundant nucleotides for this position.
                        byte nuc_one;
                        byte nuc_two;
                        Get_Two_Nucs(scores, out nuc_one, out nuc_two);

                        // Calculate Pr(obs|allele1) and Pr(obs|allele2).
                        double[][] log_nuc_probs = Obs_Probs(log_probs, log_subst_probs, obs_nucs, is_GG, reverse,
                            quality_scores, neigh_quality_scores, counts, k, nuc_one, nuc_two);

                        // Calculate log_likelihoods of genotypes for current position.
                        double[][] log_likelihoods = Log_Likelihoods(min_ploidy, max_ploidy, log_nuc_probs, nuc_props);

                        // Calculate log_likelihood of each ploidy and keep most likely allele dosage.
                        Global_Likelihood_Keep_Dose(min_ploidy, number_of_ploidies, dose_probs, global_log_like,
                            dose_queue, current_position, counts, nuc_one, nuc_two, log_likelihoods);
                    }

                    // Remove finished reads from queue. Finished reads no longer overlap with current position.
                    Padded_Read read_to_remove;
                    if (read_queue.Count > 0)
                    {
                        read_to_remove = read_queue.First();
                    }
                    else
                    {
                        read_to_remove = null;
                    }

                    while (read_to_remove != null &&
                           (read_to_remove.alignment.Pos + read_to_remove.alignment_length - 2) < current_position)
                    {
                        read_queue.Dequeue();
                        if (read_queue.Count > 0)
                        {
                            read_to_remove = read_queue.First();
                        }
                        else
                        {
                            read_to_remove = null;
                        }
                    }

                    ++current_position;
                }
                #endregion Run over every position in contig

                // Output log_likelihoods.
                int best_log_like = 0;
                for (int i = 0; i < number_of_ploidies; i++)
                {
                    if (global_log_like[i] > global_log_like[best_log_like])
                    {
                        best_log_like = i;
                    }

                    if (splitContigs)
                    {
                        writer_log_like.WriteLine("Ploidy {0} - log_likelihood {1}", i + min_ploidy, global_log_like[i]);
                    }
                }

                // Output most likely ploidy.
                int best_ploidy = best_log_like + min_ploidy;
                if (splitContigs)
                {
                    writer_ploidy.WriteLine(best_ploidy);
                }
                else
                {
                    writer_ploidy.Write("{0}\t{1}", contig_name, best_ploidy);
                    for (int i = 0; i < number_of_ploidies; i++)
                    {
                        writer_ploidy.Write("\t{0}", global_log_like[i]);
                    }
                    writer_ploidy.WriteLine("");
                }

                // Output SNPs.
                if (splitContigs)
                {
                    writer_SNP.WriteLine("Position\tAlleles\tCounts\tDosage\tPhredQuality");
                }
                char[] nuc_chars = new char[4] { 'A', 'C', 'G', 'T' };
                foreach (Best_Dose cur_doses in dose_queue)
                {
                    double cur_SNP_posterior = cur_doses.SNP_posterior[best_log_like];
                    if (cur_SNP_posterior <= log_SNP_thres)
                    {
                        int cur_best_dose = cur_doses.best_dose[best_log_like];
                        if (cur_best_dose != best_ploidy && cur_best_dose != 0)
                        {
                            if (splitContigs)
                            {
                                writer_SNP.WriteLine("{0}\t{1}|{2}\t{3}|{4}\t{5}\t{6}", cur_doses.position + 1,
                                    nuc_chars[cur_doses.nuc_one], nuc_chars[cur_doses.nuc_two], cur_doses.count_one,
                                    cur_doses.count_two, cur_best_dose, -10 * cur_SNP_posterior / Math.Log(10));
                            }
                            else
                            {
                                writer_SNP.WriteLine("{0}\t{1}\t{2}|{3}\t{4}|{5}\t{6}\t{7}", contig_name,
                                    cur_doses.position + 1, nuc_chars[cur_doses.nuc_one], nuc_chars[cur_doses.nuc_two],
                                    cur_doses.count_one, cur_doses.count_two, cur_best_dose,
                                    -10 * cur_SNP_posterior / Math.Log(10));
                            }
                        }
                    }
                }
                
                // Output read statistics.
                if (splitContigs)
                {
                    writer_reads.WriteLine("\nNumber of aligned reads: {0}", number_of_aligned_reads);
                    writer_reads.WriteLine("Number of aligned base pairs: {0}", number_of_aligned_base_pairs);
                    writer_reads.WriteLine("\nNumber of used reads: {0}", number_of_used_reads);
                    writer_reads.WriteLine("Number of used base pairs: {0}", number_of_used_base_pairs);
                }
                else
                {
                    writer_reads.WriteLine("{0}\t{1}\t{2}\t{3}\t{4}", contig_name, number_of_aligned_reads,
                        number_of_aligned_base_pairs, number_of_used_reads, number_of_used_base_pairs);
                }

                if (splitContigs)
                {
                    writer_log_like.Close();
                    writer_SNP.Close();
                    writer_ploidy.Close();
                    writer_reads.Close();
                } 
                
                clock.Stop();
                Console.WriteLine("Time to run contig: {0} s\n", (double)clock.ElapsedMilliseconds / 1000);
            }

            if (!splitContigs)
            {
                writer_SNP.Close();
                writer_ploidy.Close();
                writer_reads.Close();
            } 

            parser.Dispose();
            Console.WriteLine("Finished at {0}\n", DateTime.Now);
        }

        static void Main(string[] args)
        {
            CommandArguments.ConstructAndRun<ConPADE>(args);
        }
    }
}
