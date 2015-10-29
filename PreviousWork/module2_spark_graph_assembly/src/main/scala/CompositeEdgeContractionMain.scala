import org.apache.spark.SparkConf

import org.apache.spark.SparkContext
import org.apache.spark.SparkContext._

import org.apache.spark.rdd.RDD

import org.apache.spark.graphx._
import org.apache.spark.graphx

import scala.util.control.Breaks._
import collection.mutable.HashMap


object CompositeEdgeContraction{
  def apply(dotGraph: Graph[Int, String], debug: Boolean) = { 

    // start composite edge contraction
    println(s"====================================================")
    println(s"Composite edge contraction started!")
    println(s"Number of edges of the graph before composite edge contraction: ${dotGraph.edges.count()}")
    println(s"")

    // debug
    if (debug) {
      // print vertices and edges
      println(s"====================================================")
      println(s"dotGraph vertices:")
      dotGraph.vertices.collect.foreach(println(_))
      println(s"====================================================")
      println(s"dotGraph edges:")
      dotGraph.edges.collect.foreach(println(_))

      // the Edge class has srcId, dstId, attr members
      println(s"====================================================")
      println(s"dotGraph edges orientation info:")
      for (e <- dotGraph.edges.collect) {
        println(s"${e.srcId} --> ${e.dstId} has ${e.attr} edge properties")
      }   
    }   

    // mapReduceTriplets will be deprecated. Use aggregateMessages
    /*
    // in and out degree count and save to inOutVertexRDD (read orientation degree into in(Long), out(Long))
    val inOutVertexRDD: VertexRDD[(Long,Long)] = dotGraph.mapReduceTriplets[(Long,Long)](
      // map function: for each edge send a message to the src and dst vertices with the edge attribute
      triplet => {
        val orientation = triplet.attr.split(",")(0).toInt
        if (orientation == 1) {       // 1 = u<-------->v      reverse of u to forward of v
          Iterator((triplet.srcId, (1,0)), (triplet.dstId, (1,0)))
        }
        else if (orientation == 2) {  // 2 = u>--------<v      forward of u to reverse of v
          Iterator((triplet.srcId, (0,1)), (triplet.dstId, (0,1)))
        }
        else if (orientation == 3) {  // 3 = u>-------->v      forward of u to forware of v
          Iterator((triplet.srcId, (0,1)), (triplet.dstId, (1,0)))
        }
        else {
          Iterator.empty
        }
      },  
      // reduce function: sum in and out degree count
      (a, b) => (a._1+b._1, a._2+b._2)
    )   
    */

    // in and out degree count and save to inOutVertexRDD (read orientation degree into in(Long), out(Long))
    val inOutVertexRDD: VertexRDD[(Long,Long)] = dotGraph.aggregateMessages[(Long,Long)](
      // map function: for each edge send a message to the src and dst vertices with the edge attribute
      ctx => {
        // get orientation
        val orientation = ctx.attr.split(",")(0).toInt

        if (orientation == 1) {       // 1 = u<-------->v      reverse of u to forward of v
          //Iterator((triplet.srcId, (1,0)), (triplet.dstId, (1,0)))
          ctx.sendToSrc((1,0))
          ctx.sendToDst((1,0))
        }
        else if (orientation == 2) {  // 2 = u>--------<v      forward of u to reverse of v
          //Iterator((triplet.srcId, (0,1)), (triplet.dstId, (0,1)))
          ctx.sendToSrc((0,1))
          ctx.sendToDst((0,1))
        }
        else if (orientation == 3) {  // 3 = u>-------->v      forward of u to forware of v
          //Iterator((triplet.srcId, (0,1)), (triplet.dstId, (1,0)))
          ctx.sendToSrc((0,1))
          ctx.sendToDst((1,0))
        }
        //else {
        //  Iterator.empty
        //}
      },  
      // reduce function: sum in and out degree count
      (a, b) => (a._1+b._1, a._2+b._2)
    )   

    // debug
    if (debug) {
      println(s"====================================================")
      println(s"inOutVertexRDD VertexRDD generated:")
      inOutVertexRDD.collect.foreach(println(_))

      // check one in and one out degree vertices
      println(s"====================================================")
      println(s"Check who has one in and one out degree:")
      inOutVertexRDD.filter {
        case (id, attr) => (attr._1 == 1 && attr._2 == 1)
      }.collect.foreach {
        case (id, attr) => println(s"Vertex ${id} has one in and out degree with properties: ${attr}")
      }
    }

    // create initial inOutGraph
    var inOutGraph: Graph[(Long,Long),String] = dotGraph.mapVertices{ case (id,default) => (0,0) }

    // Merge the inOutVertexRDD values into the inOutGraph
    inOutGraph = inOutGraph.outerJoinVertices(inOutVertexRDD) {
       (vid, old, newOpt) => newOpt.getOrElse(old)
    }

    // debug
    if (debug) {
      println(s"====================================================")
      println(s"inOutGraph was generated:}")
      inOutGraph.vertices.collect.foreach(println(_))
      inOutGraph.edges.collect.foreach(println(_))
    }

    // Mark the edge to be removed and contracted.
    // If destination node has (1,1) in out degrees, then mark the edge.
    // If the edge was marked, then the edge will be removed and contracted later
    val markedGraph: Graph[(Long,Long),(String, Boolean)] = inOutGraph.mapTriplets(
      triplet => if (triplet.dstAttr._1 == 1 && triplet.dstAttr._2 == 1) (triplet.attr, true)
                 else (triplet.attr, false)
    )

    // debug
    if (debug) {
      println(s"====================================================")
      println(s"markedGraph was generated:}")
      markedGraph.vertices.collect.foreach(println(_))
      markedGraph.edges.collect.foreach(println(_))
    }

    // Restric graph to contractable (mark=true) edges
    val contractableGraph = markedGraph.subgraph(epred = e => e.attr._2 == true)

    // debug
    if (debug) {
      println(s"====================================================")
      println(s"contractableGraph was generated:}")
      contractableGraph.vertices.collect.foreach(println(_))
      contractableGraph.edges.collect.foreach(println(_))
    }

    // Compute connected component id for each vertex from the contractableGraph
    val conectedComponentVertices = contractableGraph.connectedComponents.vertices.cache()

    // debug
    if (debug) {
      println(s"====================================================")
      println(s"conectedComponentVertices (connected component IDs of contractableGraph) were generated:}")
      conectedComponentVertices.collect.foreach(println(_))
    }

    // Convert the IDs of the same connected component (_._2) to be same to merge
    val duplicateIdVertices = markedGraph.vertices.innerJoin(conectedComponentVertices) { (id, attr, cc) => (cc, attr) }.map(_._2)

    // debug
    if (debug) {
      println(s"====================================================")
      println(s"duplicateIdVertices was generated:}")
      duplicateIdVertices.collect.foreach(println(_))
    }

    // New contracted vertices after aggregating same index
    val contractedVertices: VertexRDD[(Long,Long)] = markedGraph.vertices.aggregateUsingIndex(duplicateIdVertices, (a,b) => (a._1+b._1, a._2+b._2))

    // debug
    if (debug) {
      println(s"====================================================")
      println(s"contractedVertices was generated:}")
      contractedVertices.collect.foreach(println(_))
    }

    // Restrict graph to remained (mark=false) edges
    val remainedGraph = markedGraph.subgraph(epred = e => e.attr._2 != true)

    // debug
    if (debug) {
      println(s"====================================================")
      println(s"remainedGraph was generated:}")
      remainedGraph.vertices.collect.foreach(println(_))
      remainedGraph.edges.collect.foreach(println(_))
    }

    // Contracted edges
    val contractedEdges = remainedGraph.outerJoinVertices(conectedComponentVertices) {
      (id, _, ccOpt) => ccOpt.get }
    .triplets.map { e => Edge(e.srcAttr, e.dstAttr, e.attr) }

    // debug
    if (debug) {
      println(s"====================================================")
      println(s"contractedEdges was generated:}")
      contractedEdges.collect.foreach(println(_))
    }

    // Contracted graph
    val contractedGraph = Graph(contractedVertices, contractedEdges)

    // debug
    if (debug) {
      println(s"====================================================")
      println(s"contractedGraph was generated:}")
      contractedGraph.vertices.collect.foreach(println(_))
      contractedGraph.edges.collect.foreach(println(_))
    }

    // end composite edge contraction
    println(s"====================================================")
    println(s"Composite edge contraction done!")
    println(s"Number of edges of the graph after composite edge contraction: ${contractedGraph.edges.count()}")
    println(s"")

  }
}
