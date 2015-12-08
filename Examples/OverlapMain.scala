import org.apache.spark.SparkConf
import scala.collection.mutable.ArrayBuffer
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

    //set the file to the one passed from the commandline
    val readsListFile = sc.textFile(args(0))
    //turn into key value RDD where the read is the key and the number of the read is the value (ATCG, 1)
    val keyRDD = readsListFile.map(line => (line.split(" ")(1), line.split(" ")(0)))
    //Now we need to generate a rdd for every K, V that is of type (SubString K, V) for every substring of at least len 25
    val Substrings = keyRDD.flatMap((read: (String, String)) => makeSubStr(read._1, read._2))
    //Then simply group by key and filter for values with len of 2 or more
    val overlaps = Substrings.groupByKey().filter((reads: (String, Iterable[String])) =>  reads._2.size >= 2)
	
    println(s"====================================================")
    println(s"Find Overlaps: done")

  }
  def makeSubStr(k:String, v:String) : ArrayBuffer[(String, String)] = {
    val substrings = ArrayBuffer[(String, String)]()
    for {start <- 0 to (k.length-25); end <- (start + 25) to k.length} yield //loop over all possible substrings of length 25 or greater
    substrings += ((k.substring(start, end), v)) //append the substring array
    return substrings
  }

}
