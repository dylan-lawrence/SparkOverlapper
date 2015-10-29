import org.apache.spark.SparkConf

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._

import org.apache.spark.rdd.RDD

import org.apache.spark.graphx._
import org.apache.spark.graphx

import scala.util.control.Breaks._


object typhoon{

  val usage = """
    Usage: typhoon filename
  """

  def main(args: Array[String]) {

    // debug mode can be true or false
    val debug = true

    // start transitive edge reduction
    println(s"====================================================")
    println(s"Typhoon assembler started!")
    println(s"")

    if (args.length == 0) {
      println(usage)
      System.exit(1)
    }

    // We pass the SparkContext constructor a SparkConf object which contains 
    // information about our application
    val conf = new SparkConf().setAppName("typhoon")
    // Initialize a SparkContext as part of the program
    val sc = new SparkContext(conf)

    // Intialize edgeListFile
    // val edgeListFile = "data/3reads_overlap_beforeTransitiveEdgeReduction_edge_list.txt"
    val edgeListFile = args(0)

    // Load edge list into graph
    val edges = sc.textFile(edgeListFile).flatMap { line =>
      if (!line.isEmpty && line(0) != '#') {
        val lineArray = line.split("\\s+")
        if (lineArray.length < 2) {
          None
        } else {
          val srcId = lineArray(0).toLong
          val dstId = lineArray(1).toLong
          val attr  = lineArray(2)
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
          Some(Edge(srcId, dstId, attr+":"))
        }
      } else {
        None
      }
    }

    // construct graph from edges of DOT file
    val dotGraph = Graph.fromEdges(edges, 1)

    println(s"====================================================")
    println(s"Load edge list and generate graph: Done")
    println(s"Number of edges of the graph: ${dotGraph.edges.count()}")
    println(s"")

    // transtivie edge reduction
    // val transitiveEdgeReductionGraph: Graph[Int, String] = TransitiveEdgeReduction(dotGraph, debug)
    val transitiveEdgeReductionGraph = TransitiveEdgeReduction(dotGraph, debug)

    // composite edge contraction
    val compositeEdgeContractionGraph = CompositeEdgeContraction(transitiveEdgeReductionGraph, debug)
  }


}
