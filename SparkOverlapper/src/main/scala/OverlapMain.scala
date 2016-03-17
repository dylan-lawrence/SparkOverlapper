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
  // minimum overlap length
  val minovl = 30


  def main(args: Array[String]) {

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

    // set the file to the one passed from the commandline
    val readsListFile = sc.textFile(args(0))
    // turn into key value RDD where the read is the key and the number of the read is the value (ATCG, 1)
    val keyRDD = readsListFile.map(line => (line.split(" ")(1), line.split(" ")(0)))
    // Now we need to generate a rdd for every K, V that is of type (SubString K, V) for every substring of at least len 25
    val Substrings = keyRDD.flatMap((read: (String, String)) => makeSubStr(read._1, read._2))
    // Then simply group by key and filter for values with len of 2 or more
    val overlaps = Substrings.groupByKey().filter((reads: (String, Iterable[String])) =>  reads._2.size >= 2)
    // print overlaps
    val ovrArray = overlaps.collect()//put into array
    val overlapProp = "[0-9]+".r
    for (i <- 0 to (ovrArray.length-1)){
     	 val subArray = ovrArray(i).toString.split(",C")
	 val props = overlapProp.findAllIn(subArray(1)).toList
         var source = props(0)
	 var sourceLen = props(1)
	 var overlapLen = props(2)
         var sourceStart = props(3)
	 var sourceEnd = props(4)
	 for(x <- Range(start = 5, end = props.length, step = 5)){
 	    var dest = props(x)
	    var destLen = props(x+1)
	    var destStart = props(x+3)
	    var destEnd = props(x+4)
	    if(destStart != sourceStart){
	        println(source + "\t" + dest + "\t" + "0," +overlapLen + ",0,0," + sourceLen + "," + sourceStart + "," + sourceEnd + "," + destLen + "," + destStart + "," + destEnd + ",NA")
	    }
	 }	    
    }	
    println(s"====================================================")
    println(s"Find Overlaps: done")

  }
  def makeSubStr(k:String, v:String) : ArrayBuffer[(String, String)] = {
    var keyLength = k.length
    val substrings = ArrayBuffer[(String, String)]()
    //loop over all possible substrings of length minovl or greater
    for {start <- 0 to (keyLength-minovl); end <- (start + minovl) to keyLength} yield   
    if(start==0||end==keyLength){
	var subS = (k.substring(start, end))
	substrings += ((subS,"(" + v + ":" + keyLength + ":" + subS.length + ":" + start + ":" + end + ")")) //append the substring array
   }	 
   return substrings
  }
  //def formatOut(

}
