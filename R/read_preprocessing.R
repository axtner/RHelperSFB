# system("samtools view filename.bam | awk '{OFS=\"\\t\"; print \">\"$1\"\\n\"$10}' -> filename.fasta")

# bash /home/bioadmin/Projects/BiDoup_2022/221114_M01108_0115/read_preprocessing_v2.sh /home/bioadmin/Projects/BiDoup_2022/221114_M01108_0115/Data/Intensities/BaseCalls/ M01108 WaterSeq001_results 300 5 SampleSheet.csv 0100 sample_tags

paste0("read_preprocessing_v2.sh")

