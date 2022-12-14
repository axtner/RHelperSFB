% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demultiplex.R
\name{demultiplex}
\alias{demultiplex}
\title{Demultiplex double indexed Illumina sequencing reads.}
\source{
(1) \code{bcl2fastq}: 
\url{https://emea.support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html}

(2) \code{cutadapt}:
\url{https://doi.org/10.14806/ej.17.1.200}

(3) \code{adapterremoval}:  
\url{https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-016-1900-2}

(4) \code{vsearch}: 
\url{https://github.com/torognes/vsearch}
}
\usage{
demultiplex(
  seq_dir = NA,
  seq_run = NA,
  out_dir = NA,
  loci = "0100",
  cycles = 300,
  index = 5,
  tag_dir = "sample_tags",
  seq_machine = "M01108"
)
}
\arguments{
\item{seq_dir}{Run folder of an Illumina sequencing run.}

\item{seq_run}{Sequencing run name (e.g. seq_run = "LeechSeq001").}

\item{out_dir}{Custom name of the created output folder that is created in the run folder. Default is "DerepOut_" followed by current date and time (e.g. "DerepOut_20221219_110204")}

\item{loci}{String decoding the genetic markers that were sequenced as a binary string with each locus coded as either present (1) or absent (0) in the order: 12S, 16S, CytB, COI. Default is set to 16S only loci= "0100", 12S and CytB would be 1010.}

\item{cycles}{Integer giving the number of sequencing cycles. Default is cycles = 300.}

\item{index}{Integer giving length of the indices. Default is index = 5.}

\item{tag_dir}{Directory containing a text file with a list of sample names and the sample tags for each plate tag. File name format: PlateLabel.txt (e.g. "p5001.txt"). PlateLabel must match Sample_ID in SampleSheet.csv extactly. File contents format: sample_name fwdTag revTag}

\item{seq_machine}{Name of the sequencing machine, needed for the bcl2fastq script that is deployed to convert the base calls files to fastq files. default is set to "M01108", the Illumina MiSeq machine of the Leibniz-IZW.}
}
\description{
demultiplex
}
\details{
Function deploys the read_preprocessing.sh script from the screenforbio pipeline of Alex Crampton-Platt. The original bash script has been edited to eliminate bugs. The bash script uses \code{bcl2fastq} (1), \code{cudadapt} (2), \code{adapterremoval} (3) and \code{vsearch} (4). All four programmes must be installed on the computer running the analysis and their superior directory where their executable files are saved must be added to the PATH variable of the environment (e.g. use base::Sys.setenv() to add them in R).
}
\examples{
Simple example to analyse 16S sequencing reads with using th default settings:
demultiplex(seq_dir = "/SequencingData/221114_M01108_0115/", seq_run = "SeqRun001")

}
