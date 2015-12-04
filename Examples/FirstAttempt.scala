import org.apache.spark.SparkConf

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._

import org.apache.spark.rdd.RDD

import org.apache.spark.graphx._
import org.apache.spark.graphx

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
    val readsListFile = sc.textFile(args(0))
	
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
	

    // Load reads into spark rdd
	for (line <- Source.fromFile(readsListFile).getLines()) {
		if line !startsWith(">"){
			reads(readsInd) = line
			readsInd++
		}
	}
	for (firstRead <- 0 until reads.length){
		for (secondRead <- 0 until reads.length){
			if !(firstRead == secondRead){ //makes sure were not comparing the same read
				firstArray = reads(firstSeq).toCharArray
				secondArray = reads(secondSeq).toCharArray
				for (firstChar <- 0 until firstArray.length){\
					var tempOverlapLen = 0
					if ((firstChar + 26 ) <= 35){ //The potential overlap is len 25 or greater
						i = firstChar
						for(secondChar <- 0 until secondArray.length){
							if(secondArray(secondChar) == firstArray(firstChar)){
								i++
								tempOverlapLen++
								destEnd = secondChar	
							}
						}
					if (tempOverlapLen >= 25){
						overlapLen = tempOverlapLen
						sourceStart = firstChar
						sourceEnd = firstChar + overlapLen
						destStart = destEnd - overlapLen
						done = 1
					}
					
					if(done == 1) break
					}
				}
			}
		}
	}

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
