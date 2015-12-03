import org.apache.spark.SparkConf

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._

import org.apache.spark.rdd.RDD


import scala.util.control.Breaks._
import scala.io.Source


object overlapper{

  val usage = """
    Usage: overlapper fileName
  """

  def main(args: Array[String]) {

    // debug mode can be true or false
    val debug = true

    // start transitive edge reduction
    println(s"====================================================")
    println(s"Overlapper Started")
    println(s"")

    if (args.length == 0) {
      println(usage)
      System.exit(1)
    }

    // We pass the SparkContext constructor a SparkConf object which contains 
    // information about our application
    val conf = new SparkConf().setAppName("overlapper")
    // Initialize a SparkContext as part of the program
    val sc = new SparkContext(conf)

    // Intialize readsListFile
    // val readsListFile = "10reads_forward.fasta"
    def makeSubStr(k:string) : ArrayBuffer[String] = {
    	val substrings = ArrayBuffer[String]
    	for {start <- 0 to k.length; end <- (start + 25) to s.length} yield //loop over all possible substrings of length 25 or greater
    	substrings += k.substring(start, end) //append the substring
	 return substrings
    }
    val readsListFile = sc.textFile(args(0))
    //turn into key value RDD where the read is the key and the number of the read is the value (ATCG, 1)
    val keyRDD = readsListFile.map(line => (line.split(" ")(1), line.split(" ")(0))
    //Now we need to generate a rdd for every K, V that is of type (SubString K, V) for every substring of at least len 25
    val Substrings = keyRDD.flatMap(read: (String, Int) => (read._1.makeSubStr))
    //Then simply group by key and filter for values with len of 2 or more
    val overlaps = Substrings.groupByKey().filter(reads =>  reads.length >= 2)
	
	//val source = the number of the source
	//val destination = the number of the destination
	
	var overlapLen = 0
	var sourceLen = 35
	var sourceStart = 0
	var sourceEnd = 0
	var destLen = 35
	var destStart = 0
	var destEnd = 0
	var reads : Array[String] = new Array[String](10)
	var readsInd = 0
	var done = 0
// edge has 9 attributes Ex) 3,F,33,0,0,2,34,0,32
// Col1: overlap orientation
// 0 = u<--------<v      reverse of u to reverse of v  
//   => This case is handled in DOT file preprocessing step and changed to 3 (u>-->v)
// 1 = u<-------->v      reverse of u to forward of v
// 2 = u>--------<v      forward of u to reverse of v
// 3 = u>-------->v      forward of u to forware of v
// Col2: overlap property F:forward, 
//                        FRC::read1 overlaps with the reverse complement of read2
// Col3~9: overlap length, substitutions, edits, start1, stop1, start2, stop2
// Properties (String, Boolean)

    println(s"====================================================")
    println(s"Find Overlaps: done

  }


}
