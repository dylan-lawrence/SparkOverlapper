# Typhoon
Typhoon: cloud based de novo genome assembler using Spark

# Getting Started
$ spark-submit --class "typhoon" --master local[4] target/scala-2.10/typhoon_2.10-1.0.jar data/10reads_overlap_beforeTransitiveEdgeReduction_edge_list.txt
