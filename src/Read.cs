using System;
using System.Linq;
using Bio;
using Bio.IO.SAM;

namespace ConPADE
{
    /// <summary>
    /// Class to build and store padded reads.
    /// The CIGAR string is used to pad insertions and deletions.
    /// Read positions can then be accessed based on the alignment to the reference.
    /// </summary>
    class Padded_Read
    {
        public SAMAlignedSequence alignment { get; private set; }
        public byte[] padded_sequence;
        public int[] padded_quality_scores;
        public int[] located_sequence;
        public bool is_reverse;
        public int cur_pos_ind;
        public int[] numbers;
        public char[] letters;
        public int vector_ind;
        public int alignment_length;

        /// <summary>
        /// Padded_Read standard constructor.
        /// </summary>
        /// <param name="alignment">An input object of type SAMAlignedSequence.</param>
        public Padded_Read(SAMAlignedSequence alignment)
        {
            this.alignment = alignment;
            this.is_reverse = alignment.Flag.HasFlag(SAMFlags.QueryOnReverseStrand);
            this.Pad_Read();
        }

        // Split CIGAR string into paired "numbers" and "letters".
        private void Split_CIGAR()
        {
            this.numbers = new int[1000];
            this.letters = new char[1000];
            int cigar_ind = -1;
            this.vector_ind = 0;
            while (cigar_ind < (this.alignment.CIGAR.Length - 1))
            {
                int cur_number = 0;
                int current;
                while (Int32.TryParse(this.alignment.CIGAR[++cigar_ind].ToString(), out current))
                {
                    cur_number = 10 * cur_number + current;
                }
                this.numbers[this.vector_ind] = cur_number;
                this.letters[this.vector_ind++] = this.alignment.CIGAR[cigar_ind];
            }
        }

        // Count the total length of deletions to pad.
        // Currently considering deletions ("D") and soft padding ("P")
        private int Count_Padding()
        {
            int extras = 0;
            for (int i = 0; i < this.numbers.Length; i++)
            {
                if (this.letters[i] == 'D' || this.letters[i] == 'P')
                {
                    extras += this.numbers[i];
                }
            }
            return extras;
        }

        /// <summary>
        /// Pad deletions in an aligned read.
        /// Create an index of insertions/deletions so that aligned
        /// positions may be accessed based on reference index.
        /// </summary>
        public void Pad_Read()
        {
            byte[] sequence = this.alignment.QuerySequence.ToArray();
            int[] quality_scores = (this.alignment.QuerySequence as QualitativeSequence).GetQualityScores();

            this.Split_CIGAR();

            int extras = Count_Padding();

            // size includes deletions and one spacer on both sides
            int size = sequence.Length + extras + 2;
            this.padded_sequence = new byte[size];
            this.padded_quality_scores = new int[size];
            this.located_sequence = new int[size];

            // Values to use for spacers
            int no_qual = -10;
            byte spacer_nuc = 5;
            int spacer_ind = -1;

            int last_ind = size - 1;
            this.padded_sequence[0] = spacer_nuc;
            this.padded_sequence[last_ind] = spacer_nuc;
            this.padded_quality_scores[0] = no_qual;
            this.padded_quality_scores[last_ind] = no_qual;
            this.located_sequence[0] = spacer_ind;
            this.located_sequence[last_ind] = spacer_ind;

            int k = 1;
            int l = 0;
            int m = 1;

            //         Nucleotide : Index
            //                  A : 0
            //                  C : 1
            //                  G : 2
            //                  T : 3
            //                  N : 4
            // Others (incl dels) : 5
            // We may want to use a separate index for deletions, if we want to call indels in the future

            int last = spacer_ind;
            for (int i = 0; i < this.letters.Length; i++)
            {
                char cur_letter = this.letters[i];
                if (cur_letter == 'D' ||
                    cur_letter == 'P')
                {
                    if (cur_letter == 'D')
                    {
                        for (int j = 0; j < this.numbers[i]; j++)
                        {
                            this.located_sequence[m++] = ++last;
                            this.padded_quality_scores[k] = no_qual;
                            this.padded_sequence[k++] = spacer_nuc;
                        }
                    }
                }
                else
                {
                    if (cur_letter == 'M' ||
                        cur_letter == '=' ||
                        cur_letter == 'X')
                    {
                        for (int j = 0; j < this.numbers[i]; j++)
                        {
                            this.located_sequence[m++] = ++last;
                        }
                    }
                    else if (cur_letter == 'I')
                    {
                        if (last == spacer_ind)
                        {
                            ++last;
                        }
                        for (int j = 0; j < this.numbers[i]; j++)
                        {
                            this.located_sequence[m++] = last;
                        }
                    }
                    else if (cur_letter == 'S')
                    {
                        for (int j = 0; j < this.numbers[i]; j++)
                        {
                            this.located_sequence[m++] = spacer_ind;
                        }
                    }

                    if (cur_letter != 'H') {
                        for (int j = 0; j < this.numbers[i]; j++)
                        {
                            this.padded_quality_scores[k] = quality_scores[l];
                            switch (sequence[l++])
                            {
                                case 65:
                                    this.padded_sequence[k++] = 0;
                                    break;
                                case 67:
                                    this.padded_sequence[k++] = 1;
                                    break;
                                case 71:
                                    this.padded_sequence[k++] = 2;
                                    break;
                                case 84:
                                    this.padded_sequence[k++] = 3;
                                    break;
                                case 78:
                                    this.padded_sequence[k++] = 4;
                                    break;
                                default:
                                    this.padded_sequence[k++] = 5;
                                    break;
                            };
                        }
                    }
                }
            }

            // Set cur_pos_ind to the position with the first aligned nucleotide.
            // This value will be updated by other methods.
            this.cur_pos_ind = 0;
            while (this.located_sequence[++this.cur_pos_ind] == spacer_ind) ;
            
            this.alignment_length = this.alignment.RefEndPos - this.alignment.Pos + 1;
        }
    }
}