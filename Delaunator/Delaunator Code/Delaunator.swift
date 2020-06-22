//
//  triangulator.swift
//  Delaunator_Swift
//
//  Created by Z Chameleon on 4/5/20.
//  Copyright © 2020 Z Cha. All rights reserved.
//
import Foundation

// Define statics
struct Static {
  static var Epsilon:Double = pow(2.0, -52)
  static var Tolerance:Double = 1.0e-12 // Arbitrary number
  static let Threshold:Double = 3.3306690738754716e-16
}

// Write to the error stream
final class StandardErrorOutputStream: TextOutputStream {
  func write(_ string: String) {
    FileHandle.standardError.write(Data(string.utf8))
  }
}
var errorStream = StandardErrorOutputStream()

// Start with some global constants
var EDGE_STACK  = Array(repeating:0, count:512)


// Define a simple point struct
struct Point: Hashable, Codable {
  var x, y: Double
}


// For recording
struct Triangulation: Hashable, Codable, Identifiable {

  var points: Array<Point>
  var triangles: Array<Int>
  var halfEdges: Array<Int>
  var hull: Array<Int>
  var numberEdges: Int
  
  // Extras for display code
  var id:UUID? = UUID()
  
  
  // Provided init
  init(using delaunay:Delaunator_Swift, with points:Array<Point>) {
    triangles   = delaunay.triangles
    halfEdges   = delaunay.halfEdges //.map({Index(id:$0)})
    hull        = delaunay.hull //.map({Index(id:$0)})
    numberEdges = delaunay.numberEdges
    self.points = points
  }
  
  // Empty triangulation
  init () {
    triangles = [Int]()
    halfEdges = [Int]()
    hull = [Int]()
    points = [Point]()
    numberEdges = 0
  }
}

/*:
 ## Delaunator algorithm
 
 The  triangulation is based on the [Sweep Circle](http://cglab.ca/~biniaz/papers/Sweep%20Circle.pdf) algorithm of
 Biniaz & Dastghaibyfard (2011).
 
 The idea is for points to be included in the triangulation as they enter an expanding circle centred somewhere on the plane.
 ![Sweep Circle Example](./Documentation/Sweep Circle Algorithm.png)
 
 The new triangle is made locally delaunay - through edge flippping, and then once this is done
 the advancing front of the conex hull must be updated. This can be done by searching both
 anti-clockwise and clockwise from the new point, adding extra edges as needed. To facilitate this
 a list of points is maintained sorted by angle from an initial point inside  the convex hull.
 
 This initial seed point need not be the same as the radial origin, but can be if the radial point is
 within the convex hull. To [ensure this](https://cs.stackexchange.com/questions/63362/delaunay-sweep-circle-initialization)
 a good seed point is one is to take a point from the set, its nearest neighbour, and a 3rd point that creates the smallest circumcircle with them will always work
 as this is a delaunay triangle (nearest-neighbours are always delaunay edges, and the 2 triangles based on this edge have to create an empty circumcircle, so
 one of them has to be the smallest) which makes the circumcircle empty of other points, so using its origin as the radial origin is correct.
 But (in principle) the circumcentre might be outside the convex hull; so using the incircle centre would work to for the radial origin & the angular origin.
 
 
 */

// The main structure
struct Delaunator_Swift {
  var triangles: Array<Int>
  var halfEdges: Array<Int>
  var hull: Array<Int>
  
  var numberEdges: Int
  fileprivate var numberPoints, maxTriangles,  hashSize, hullStart, hullSize: Int
  fileprivate  var hullTri: Array<Int>
  fileprivate var hullPrev: Array<Int>
  fileprivate var centre: Point
  fileprivate  var coords = [Double]()
  
  // Use a provided init
  init(from points: Array<Point>) {
    numberPoints = points.count
    coords = Array(repeating:0.0, count:2*numberPoints)
    
    // arrays that will store the triangulation graph
    maxTriangles = max(2 * numberPoints - 5, 0)
    numberEdges = 0
    triangles = Array(repeating:0, count:3 * maxTriangles)
    halfEdges = Array(repeating:0, count:3 * maxTriangles)
    
    hullTri  = Array(repeating:0, count:numberPoints) // edge to adjacent triangle
    hullPrev = Array(repeating:0, count:numberPoints) // edge to prev edge
    hull     = Array(repeating:0, count:numberPoints)
    
    // The size allowed for the hash arrau ordered by angle
    // Since it is for  the hull it is approx. the square root of the total
    // number of points
    hashSize = Int(Double(numberPoints).squareRoot().rounded(.up))
    hullStart = 0
    hullSize = 0
    centre = Point(x:0.0, y:0.0)
    
    // Back door escape if numberPoints is zero
    if 0 == numberPoints {return}
    
    var hullNext = Array(repeating:0, count:numberPoints) // edge to next edge
    var hullHash = Array(repeating:-1, count:hashSize) // angular edge hash
    
    // Distances from the origin plus index ids
    var dists = Array(repeating:0.0, count:numberPoints)
    var ids   = Array(repeating:0, count:numberPoints)
    
    // populate an array of point indices; calculate input data bbox
    var minX =  Double.infinity
    var maxX = -Double.infinity
    var minY =  Double.infinity
    var maxY = -Double.infinity
    
    // Compute bounds and set the coords array
    for (i, p) in points.enumerated() {
      // Set each co-ordinate
      coords[2 * i]     = p.x
      coords[2 * i + 1] = p.y
      if (p.x < minX)  {minX = p.x}
      if (p.y < minY)  {minY = p.y}
      if (p.x > maxX)  {maxX = p.x}
      if (p.y > maxY)  {maxY = p.y}
      ids[i] = i
    }
    
    // Initially find the point and its nearest neighbour closest the centre
    // of the bounds and then the third point that forms s triangle with the
    // smallest circum-radius. This is a locally delaunay triangle and forms
    // the origin both fo rth eradial distances and azimuthal angles
    //
    // Another option would be to use the centre of this triangle's in-circle.
    //
    // Centre
    let c = Point(x: 0.5 * (minX + maxX), y: 0.5 * (minY + maxY))
    
    var minDist = Double.infinity
    var i0 = 0, i1 = 0, i2 = 0
    
    // pick a seed point close to the centre
    for i in 0..<numberPoints {
      let d = dist(ax: c.x, ay: c.y, bx: coords[2 * i], by: coords[2 * i + 1])
      if d < minDist {
        i0 = i
        minDist = d
      }
    }
    let i0x = coords[2 * i0]
    let i0y = coords[2 * i0 + 1]
    
    minDist = Double.infinity
    
    // find the point closest to the seed
    for i in 0..<numberPoints {
      if i == i0 { continue }
      let d = dist(ax: i0x, ay: i0y, bx: coords[2 * i], by: coords[2 * i + 1])
      
      // d is the distance squared - so if it is very small th epoints are very close
      if (d < minDist && d > Static.Epsilon) {
        i1 = i
        minDist = d
      }
    }
    var i1x = coords[2 * i1]
    var i1y = coords[2 * i1 + 1]
    
    var minRadius = Double.infinity
    
    // find the third point which forms the smallest circumCircle with the first two
    for i in 0..<numberPoints {
      if (i == i0) || (i == i1) { continue }
      let r = circumRadius(ax: i0x, ay: i0y, bx: i1x, by: i1y,
                           cx: coords[2 * i], cy: coords[2 * i + 1])
      if r < minRadius {
        i2 = i
        minRadius = r
      }
    }
    var i2x = coords[2 * i2]
    var i2y = coords[2 * i2 + 1]
    
    
    // All the possible circumcircles have a
    // a centre at infinity - so all the points are colinear
    if minRadius == Double.infinity {
      // order colinear points by dx (or dy if all x are identical)
      // and return the list as a hull
      for i in 0..<numberPoints {
        // Save the x-distance unless its zero and then save the y-difference instead
        let deltaX = coords[2 * i] - coords[0]
        if (is_near_zero(x: deltaX)) {
          dists[i] = coords[2 * i + 1] - coords[1]
        } else {
          dists[i] = deltaX
        }
      }
      
      // This is a quicksort
      // This sorts the indices in place
      quicksortRandom(&ids, using:dists, low:0, high:numberPoints - 1)
      var j = 0
      var d0 = -Double.infinity
      
      
      /*
       This seems to do this
       Which I don't quite understand
       
       for id in ids {
       if (dists[id] > d0) {
       hull.append(id)
       d0 = dists[id]
       }
       }
       
       */
      
      for i in 0..<numberPoints {
        let id = ids[i]
        if (dists[id] > d0) {
          hull[j] = id
          j += 1
          d0 = dists[id]
        }
      }
      hull.removeLast(numberPoints - j)
      triangles = []
      halfEdges = []
      return
    }
    
    // swap the order of the seed points for anticlockwise orientation
    if (!orient(rx: i0x, ry:i0y, qx: i1x, qy: i1y, px: i2x, py: i2y)) {
      // Swap i1 & i2
      (i1, i2) = (i2, i1)
      (i1x, i2x) = (i2x, i1x)
      (i1y, i2y) = (i2y, i1y)
    }
    
    // All the work above was to initialize the seed triangle
    // The circum-centre of this traingle is the origin for the
    // radial and azimuthal sorts
    
    // Get the circumcentre of this triangle
    centre = circumCentre(first:  Point(x:i0x, y:i0y),
                          second: Point(x:i1x, y:i1y),
                          third:  Point(x:i2x, y:i2y))
    
    // Get the radial distance of each point from th ecentre
    for i in 0..<numberPoints {
      dists[i] = dist(ax: coords[2 * i], ay:coords[2 * i + 1],
                      bx:centre.x, by:centre.y)
    }
    
    // sort the points by distance from the seed triangle circumCentre
    quicksortRandom(&ids, using:dists, low:0, high:numberPoints - 1)
    
    // set up the seed triangle as the starting hull
    hullStart = i0
    hullSize = 3
    
    // hullNext[i] -> points to next index on the hull in an anti-clockwise direction
    // hullPrev[i] -> points to next index on the hull in an      clockwise direction
    hullNext[i0] = i1; hullPrev[i2] = i1
    hullNext[i1] = i2; hullPrev[i0] = i2
    hullNext[i2] = i0; hullPrev[i1] = i0
    
    // The initial triangle defines the hull
    hullTri[i0] = 0
    hullTri[i1] = 1
    hullTri[i2] = 2
    
    // The hash is based on a pseudo-angle so it increases with azimuth
    // So
    hullHash[hashKey(p:Point(x: i0x, y: i0y))] = i0
    hullHash[hashKey(p:Point(x: i1x, y: i1y))] = i1
    hullHash[hashKey(p:Point(x: i2x, y: i2y))] = i2
    
    // Add the first triangle!
    numberEdges = 0 // This is not needed...
    
    // The first trianle is th econvex hull so all linked edges ar (-1)
    let _ = addTriangle(i0, i1, i2, -1, -1, -1)
    
    // Now start building the triangulation
    // The previous point
    var xp = 0.0, yp = 0.0
    
    // Labels! Just like good old FORTRAN IV
    projectPoints: for (k, i) in ids.enumerated() {
      //for k in 0..<ids.count {
      //let i = ids[k]
      let x = coords[2 * i]
      let y = coords[2 * i + 1]
      
      // skip near-duplicate points
      if (k > 0 && is_near_zero(x:x-xp) && is_near_zero(x:y-yp)) {continue projectPoints}
      // if (k > 0 && dist(ax: x,  ay: y, bx: xp, by: yp) < Static.shortestEdgeSquared) {continue}
      
      // Save previous point
      xp = x
      yp = y
      
      // skip seed triangle points
      if (i == i0 || i == i1 || i == i2) {continue projectPoints}
      
      // find a visible edge on the convex hull using edge hash
      // This is the projection step in the algorithm
      // Get the current pseudo-angle hash key
      var start = 0
      let key = hashKey(p:Point(x: x, y: y))
      
      // The projection of the point on the convex hull
      // will cross an edge - find the start of this edge
      for j in 0..<hashSize {
        // Loop over all the hash keys
        start = hullHash[(key + j) % hashSize]
        
        // Eventually get to a non-empty hash entry
        // Then check the next entry - if it is identical keep going because
        // that means the point was removed
        // When one is found it is the next point on the hull with a greater azimuth
        if (start != -1 && start != hullNext[start]) {break}
      }
      
      // get the previous point - we have bracketed  the projected point
      // azimuth(e) <= azimuth(i) <= azimuth(q)
      var e = hullPrev[start]
      var q = start
      start = e
      
      // Check for positive (anti-clockwise) orientation of the triangle [i, q, e] where the projected point is i
      // So if [i, q, e] is clockwise this test fails and the triangle is replaced in the loop
      // So
      while (!orient(rx: x, ry: y, qx: coords[2 * q], qy: coords[2 * q + 1], px: coords[2 * e], py: coords[2 * e + 1])) {
        if (q == start) {
          // Can't add this point! All colinear input?
          break projectPoints
        }
        
        // Move around  the hull by one triangle until it [i, e, q] is clockwise
        //
        //       i            i                 i
        //       |           /|                /|
        //       |          / |               / |
        // n <-- q     =>  n--q  relabel =>  q--e
        //       |
        //       |
        //       e <--
        //   +ve <== Direction
        e = q
        q = hullNext[e]
      }
            
      // add the first triangle from the point    i
      // This is actually just [e, i, q]         / \
      //                                        q---e
      //
      // The linked edges are -1 (outside the hull) for (e->i) and (i->q) and hillTri[e] for (q->e)
      var t = addTriangle(e, i, q, -1, -1, hullTri[e])

      // Number of edges has increased by 3; t is the old number of triangles
      // Make the new triangle locally delaunany
      hullTri[i] = legalize(check: t + 2)
      hullTri[e] = t  // keep track of boundary triangles on the hull
      
      // We added one triangle - so the hull size goes up by one
      //
      //     i
      //    / \
      //   q---e        The hull edge (q-e)  has been replaced by two new edges  e (e-i) and i (i-q)
      //    \ /
      //     p
      //
      // Replaced hullTri[e] and added hullTri[i]
      // Notice that hullTri has lots of zero entries - its size is not being changed here, just some zeros are being replaced
      // It might be better if they were initialized with a negative number
      hullSize += 1
      
      // Move forward around the hull
      //
      //        i             i
      //       / \           / \
      //      /   \         /   \
      //    q ---- n   <=  q ---- e
      //        Anticlockwise
      var n = q
      q = hullNext[q]

      // Walk around anti-clockwise; add the triangle [n, i, q]
      //
      //
      // This orientation test is on the triangle [n, i, q] => true if anticlockwise
      while (orient(rx: x, ry: y, qx: coords[2 * q], qy: coords[2 * q + 1], px: coords[2 * n], py: coords[2 * n + 1])) {

        // Add this triangle to the triangulation
        t = addTriangle(n, i, q, hullTri[i], -1, hullTri[n])
        
        // Fix locally non-delaunay edges
        // Returns hull edge i on return (the new edge)
        hullTri[i] = legalize(check: t + 2)

        // This indicates the edge n is no longer on the hull
        hullNext[n] = n // mark as removed

        // Reduced hullSize by one
        hullSize -= 1

        // Continue going anti-clockwise
        n = q
        q = hullNext[n]
      }
      
      // walk backward from the other side, adding more triangles
      // This is a clockwise search
      // Only carried out when the projected point i split
      // the edge e->q to create an anti-clockwise triangle [e, i, q]
      // Because otherwise we have
      //
      //       i
      //       |
      //       |
      // n <-- q
      //       |
      //       |
      //       e <--
      //   +ve <== Direction
      //
      // And so the point e is not visible from i
      if (e == start) {
        
        // Move backward around the hull
        //
        //        i             i
        //       / \           / \
        //      /   \         /   \
        //    q ---- e   =>  e ---- q
        //         Clockwise
        q = hullPrev[e]
        
        // This orientation test is on the triangle [i, e, q] => true if anticlockwise
        while (orient(rx: x, ry: y, qx: coords[2 * e], qy: coords[2 * e + 1], px: coords[2 * q], py: coords[2 * q + 1])) {
          
          // Add this triangle to the triangulation
          t = addTriangle(q, i, e, -1, hullTri[e], hullTri[q])
          
          // Fix locally non-delaunay edges
          // Ignore return value (which is edge i - not on the hull)
          let _ = legalize(check:t + 2)
          
          // Add new hull edge q
          hullTri[q] = t
          
          // This indicates the edge is no longer on the hull
          hullNext[e] = e // mark as removed
          
          // Reduced hullSize by one
          hullSize -= 1
          
          // Continue going clockwise
          e = q
          q = hullPrev[e]
        }
      }
      
      // update the hull indices
      hullStart   = e
      hullPrev[i] = e ; hullNext[e] = i
      hullPrev[n] = i ; hullNext[i] = n

      // save the two new edges in the hash table
      hullHash[hashKey(p:Point(x: x, y: y))] = i
      hullHash[hashKey(p:Point(x:coords[2 * e], y:coords[2 * e + 1]))] = e
    }
    
    // Query how does hullTri differ from hull?
    hull = Array(repeating: 0, count: hullSize)
    var e = hullStart
    for i in 0..<hullSize {
      hull[i] = e
      e = hullNext[e]
    }
    
    // Remove unused entries
    triangles.removeLast(triangles.count - numberEdges)
    halfEdges.removeLast(halfEdges.count - numberEdges)
  }
  
  
  mutating func legalize(check edge: Int) -> Int {
    var i = 0, ar = 0
    var a = edge
    
    // recursion eliminated with a fixed-size stack
    // Maintains a stack of edges for flipping
    // If an edge needs flipping adds it to the stack - retains edge
    
    flipEdge: while (true) {
      let b:Int = halfEdges[a]
      
      /* if the pair of triangles doesn't satisfy the Delaunay condition
       * (p1 is inside the circumCircle of [p0, pl, pr]), flip them,
       * then do the same check/flip recursively for the new pair of triangles
       *
       *           pl                    pl
       *          /||\                  /  \
       *       al/ || \bl            al/    \a
       *        /  ||  \              /      \
       *       /  a||b  \    flip    /___ar___\
       *     p0\   ||   /p1   =>   p0\---bl---/p1
       *        \  ||  /              \      /
       *       ar\ || /br             b\    /br
       *          \||/                  \  /
       *           pr                    pr
       *
       *
       *    p0 -> p0 or q3    a  -> a     b  -> b
       *    pr -> p1 or q2    al -> a1    bl -> b1
       *    p1 -> p3 or q0    ar -> a2    br -> b2
       *    pl -> p2 or q1
       *
       *    flip
       *
       *    triangles[al] -> p3 or q0
       *    triangles[b1] -> p0 or q3
       *    halfEdges[a]  -> halfEdges[b1]
       *    halfEdges[b]  -> halfEdges[a1]
       *    halfEdges[a1] <-> halfEdges[b1]
       *
       
       */
      
      
      // New triangle is triangle [a]
      // No way to know which of the [0, 1, 2] edges of b were cut in the projection
      // so complex remainder trickery to get the edge ids
      // Once flipped _al_ and _b_ are boundary edges - so no need to check delaunay status
      // From symmetry surely both _a_ and _br_ should be labelled as non-delaunay?
      
      

      
      // First case is if (a <-> b) is the convex hell edge
      // In this case no flip can occur
      if (b == -1) {
        
        // If no edges left on stack all edges are locally delaunay
        if (i == 0) {break flipEdge}
        
        // Get the next pending edge
        i -= 1
        a = EDGE_STACK[i]
        continue flipEdge
      }
      
      // Get the
      let a0 = a - a % 3
      ar = a0 + (a + 2) % 3
      
      let b0 = b - b % 3
      let al = a0 + (a + 1) % 3
      let bl = b0 + (b + 2) % 3
      
      let p0 = triangles[ar]
      let pr = triangles[a]
      let pl = triangles[al]
      let p1 = triangles[bl]
      
      let illegal = inCircle(
        ax: coords[2 * p0], ay: coords[2 * p0 + 1],
        bx: coords[2 * pr], by: coords[2 * pr + 1],
        cx: coords[2 * pl], cy: coords[2 * pl + 1],
        px: coords[2 * p1], py: coords[2 * p1 + 1])
      
      if (illegal) {
        triangles[a] = p1
        triangles[b] = p0
        
        let hbl = halfEdges[bl]
        
        // edge swapped on the other side of the hull (rare); fix the halfEdge reference
        // Repairing the hull when a flip has disturbed it
        // _bl_ is the only edge this can happen with - because _br_ is not swapped, and _ar_ and _al_ are on
        if (hbl == -1) {
          var e = hullStart
          repairHull: repeat {
            if (hullTri[e] == bl) {
              hullTri[e] = a
              break repairHull
            }
            e = hullPrev[e]
          } while (e != hullStart)
        }
        
        link(a, hbl)
        link(b, halfEdges[ar])
        link(ar, bl)
        
        // let br = b0 + (b + 1) % 3
        
        // don't worry about hitting the cap: it can only happen on extremely degenerate input
        // Fixme
        if (i < EDGE_STACK.count) {
          EDGE_STACK[i] = b0 + (b + 1) % 3 // br
          i += 1
        }
      } else {
        // Get here when an edge (a) is locally delaunay
        
        
        if (i == 0) {break flipEdge}
        i -= 1
        a = EDGE_STACK[i]
      }
    }
    
    // Only gets here when no flip - because stack is not empty after a flip
    return ar
  }
  
  // Get a hash key - need to do  lot of type casting
  func hashKey(p:Point) -> Int {
    return Int((Double(hashSize) * pseudoAngle(dx:p.x - centre.x, dy:p.y - centre.y)).rounded(.down)) % hashSize
  }
  
  
  mutating func link(_ a: Int, _ b: Int) {
    halfEdges[a] = b
    //logTriangles(e: a)
    
    if (b != -1) {
      halfEdges[b] = a
      //logTriangles(e: b)
    }
    
  }
  
  // add a new triangle given vertex indices and adjacent half-edge ids
  mutating func addTriangle(_ i0: Int, _ i1: Int, _ i2: Int,
                            _  a: Int, _  b: Int, _  c: Int) -> Int {
    
    // Save the value of numberEdges
    let t = numberEdges
    
    // Add the triangle [i0, i1, i2]
    // These actually indicate halfEdge ids
    triangles[t] = i0
    triangles[t + 1] = i1
    triangles[t + 2] = i2
    
    // Connect the paired halfEdges
    link(t, a);
    link(t + 1, b)
    link(t + 2, c)
    
    // Added three extra halfEdges
    numberEdges += 3
    
    // Return the original number of halfEdges
    return t
  }
  
  func logTriangles(e: Int) {
    let t = triangles[e]
    
    if 0 == e % 3 {
      print(String(format:"\nTriangle %4d", e/3), to: &errorStream)
    }
    print (String(format:"  halfEdge[%04d] => %04d triangle[%04d] => %04d point => [%7.2f, %7.2f]",
                  e, halfEdges[e], e, t, coords[2*t], coords[2*t+1]), to: &errorStream)
  }
  
}


// Check for a valid triangulation which conserves area and has a convex hull
func testEqual(_ phrase: String,
               _ points: Array<Point>,
               _ t: Delaunator_Swift? = nil,
               _ array: Array<Int>?) -> Bool {
  let d = t ?? Delaunator_Swift(from: points)
  assert(d.triangles == [], "Incorrect Triangulation")
  assert(d.hull == array ?? [], "Incorrect hull")
  
  // All ok if we got here
  print(phrase, to:&errorStream)
  
  return true
}

// Print out some simple summaries
func logTriangulation(_ phrase: String,
                      _ points: Array<Point>,
                      _ t: Delaunator_Swift? = nil,
                      _ array: Array<Int>?) -> Bool {
  // Get delaunay triangulation
  let d = t ?? Delaunator_Swift(from: points)
  
  // Log triangles and half edges
  print (String(format: "Number Half Edges => %04d", d.numberEdges), to: &errorStream)
  for (e, _) in d.halfEdges.enumerated()  {
    d.logTriangles(e: e)
  }
  
  print (String(format: "Number Points => %04d", points.count), to: &errorStream)
  for (i, p) in points.enumerated() {
    print (String(format: "Point[%4d] => [%7.2g, %7.2g]", i, p.x, p.y),
           to: &errorStream)
  }
  return true
}

// Check for a valid triangulation which conserves area and has a convex hull
func validateTriangulation(_ phrase: String,
                           _ points: Array<Point>,
                           _ t: Delaunator_Swift? = nil,
                           _ array: Array<Int>? = nil) -> Bool {
  let d = t ?? Delaunator_Swift(from: points)
  
  // validate halfedges
  for (i, i2) in d.halfEdges.enumerated() {
    if (i2 != -1 && d.halfEdges[i2] != i) {
      print("Invalid halfEdge connexion", to:&errorStream)
      print(String(format:"\thalfEdge[%04d] => %04d", i, i2), to:&errorStream)
      print(String(format:"\thalfEdge[%04d] => %04d", i2, d.halfEdges[i2]), to:&errorStream)
      
      return false
    }
  }
  print("halfEdges are valid", to: &errorStream)
  
  // validate triangulation
  // This library enforces a convex hull
  var hullAreas = [Double]()
  var j = d.hull.count - 1
  for (i, pi) in d.hull.enumerated() {
    let pj = d.hull[j]
    hullAreas.append((points[pi].x - points[pj].x) * (points[pi].y + points[pj].y))
    
    // Surely i == j + 1 and why is it (j + 3) here and not (j +2)??
    let c = orient(rx:points[d.hull[j]].x, ry:points[d.hull[j]].y,
                   qx:points[d.hull[(j + 1) % d.hull.count]].x,
                   qy:points[d.hull[(j + 1) % d.hull.count]].y,
                   px:points[d.hull[(j + 3) % d.hull.count]].x,
                   py:points[d.hull[(j + 3) % d.hull.count]].y)
    if !c {
      print(String(format:"hull is not convex at %d", j), to: &errorStream)
      return false
    }
    j = i
  }
  let hullArea = sum(x:hullAreas)
  
  var triangleAreas = [Double]()
  for i in stride(from:0, to:d.triangles.count, by:3) {
    let a = points[d.triangles[i]]
    let b = points[d.triangles[i + 1]]
    let c = points[d.triangles[i + 2]]
    let area = abs((b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y))
    triangleAreas.append(area)
  }
  let trianglesArea = sum(x:triangleAreas)
  
  let err = abs((hullArea - trianglesArea) / hullArea)
  if (err <= Static.Epsilon) {
    print(String(format:"triangulation is valid => %f error", err), to:&errorStream)
  } else {
    print(String(format:"triangulation is invalid => %f error", err), to:&errorStream)
    return false
  }
  
  // All ok if we got here
  print(phrase, to:&errorStream)
  return true
}


/* Helper Functions
 dist:
 inCircle:
 circumCentre:
 circumRadius:
 
 */
func dist(ax: (Double), ay: (Double),
          bx: (Double), by: (Double)) -> Double {
  let dx = ax - bx
  let dy = ay - by
  return dx * dx + dy * dy
}

func inCircle(ax: (Double), ay: (Double),
              bx: (Double), by: (Double),
              cx: (Double), cy: (Double),
              px: (Double), py: (Double)) -> Bool {
  let dx = ax - px
  let dy = ay - py
  let ex = bx - px
  let ey = by - py
  let fx = cx - px
  let fy = cy - py
  
  let ap = dx * dx + dy * dy
  let bp = ex * ex + ey * ey
  let cp = fx * fx + fy * fy
  
  let result = dx * (ey * cp - bp * fy) -
    dy * (ex * cp - bp * fx) +
    ap * (ex * fy - ey * fx)
  return result < -Static.Tolerance
}


func circumRadius(ax: Double, ay: Double,
                  bx: Double, by: Double,
                  cx: Double, cy: Double) -> Double {
  let dx = bx - ax
  let dy = by - ay
  let ex = cx - ax
  let ey = cy - ay
  
  let bl = dx * dx + dy * dy
  let cl = ex * ex + ey * ey
  let d = 0.5 / (dx * ey - dy * ex)
  
  let x = (ey * bl - dy * cl) * d
  let y = (dx * cl - ex * bl) * d
  let r = x * x + y * y
  
  // Is zero present in this list?
  if ([bl, cl, d, r].contains {is_near_zero(x: $0)} ) {
    return 0.0
  }
  
  return r
}

// monotonically increases with real angle, but doesn't need expensive trigonometry
func pseudoAngle(dx: Double, dy: Double) -> Double {
  let p = dx / (abs(dx) + abs(dy))
  return (dy > 0.0 ? 3.0 - p : 1.0 + p) / 4.0 // [0..1]
}
//
//// return 2d orientation sign if we're confident in it through J. Shewchuk's error bound check
//func orientIfSure(px: (Double), py: (Double),
//                  rx: (Double), ry: (Double),
//                  qx: (Double), qy: (Double)) -> Double {
//  let l = (ry - py) * (qx - px);
//  let r = (rx - px) * (qy - py);
//  return abs(l - r) >= Static.Threshold * abs(l + r) ? l - r : 0.0
//}
//
//// a more robust orientation test that's stable in a given triangle (to fix robustness issues)
//func orient(rx: (Double), ry: (Double),
//            qx: (Double), qy: (Double),
//            px: (Double), py: (Double)) -> Bool {
//  // First test
//  var sense = orientIfSure(px: px, py: py, rx: rx, ry: ry, qx: qx, qy: qy)
//  if (sense != 0.0) {return sense < -Static.Tolerance}
//
//  // Second test
//  sense = orientIfSure(px: rx, py: ry, rx: qx, ry: qy, qx: px, qy: py)
//  if (sense != 0.0) {return sense < -Static.Tolerance}
//
//  // Final test
//  sense = orientIfSure(px: qx, py: qy, rx: px, ry: py, qx: rx, qy: ry)
//
//  // Final case
//  return (sense < -Static.Tolerance)
//}

// return 2d orientation sign if we're confident in it through J. Shewchuk's error bound check
//
// Is a point p to the left or right of a vector r->q
// Or more symmetrically are the points (r->q->p) in an anti-clockwise order? +ve yes -ve no
//
/**
 [Robust Predicates](http://www.cs.cmu.edu/~quake/robust.html)
 */
func orientIfSure(px: (Double), py: (Double),
                  rx: (Double), ry: (Double),
                  qx: (Double), qy: (Double)) -> Double {
  // Points a, b, c
  // p => c
  // r => b
  // q => a
  
  // The left and right determinants
  let detLeft  = (ry - py) * (qx - px)
  let detRight = (rx - px) * (qy - py)
  
  // We only need to determine the sign +ve/zero/-ve
  // Positive - anticlockwise
  // Zero     - colinear
  // Negative - clockwise
  let det = detLeft - detRight
  return abs(det) >= Static.Threshold * abs(detLeft + detRight) ? det : 0
}



// a more robust orientation test that's stable in a given triangle (to fix robustness issues)
// This returns true  (anticlockwise)
//              false (clockwise)
func orient(rx: (Double), ry: (Double),
            qx: (Double), qy: (Double),
            px: (Double), py: (Double)) -> Bool {
  // First test triangle [p, r, q]
  var sense = orientIfSure(px: px, py: py, rx: rx, ry: ry, qx: qx, qy: qy)
  if (sense != 0) {return sense > 0}
  
  // Second test triangle [r, q, p]
  sense = orientIfSure(px: rx, py: ry, rx: qx, ry: qy, qx: px, qy: py)
  if (sense != 0) {return sense > 0}
  
  // Final test triangle [q, p, r]
  // If sense == 0 it was so for all the tests - so in this case assume true if zero
  sense = orientIfSure(px: qx, py: qy, rx: px, ry: py, qx: rx, qy: ry)
  
  // Final case
  return (sense > 0)
}


/*
 Lomuto's partitioning algorithm.
 This is conceptually simpler than Hoare's original scheme but less efficient.
 The return value is the index of the pivot element in the new array. The left
 partition is [low...p-1]; the right partition is [p+1...high], where p is the
 return value.
 The left partition includes all values smaller than or equal to the pivot, so
 if the pivot value occurs more than once, its duplicates will be found in the
 left partition.
 */
func partitionLomuto<T: Comparable>(_ index: inout [Int], _ a: [T], low: Int, high: Int) -> Int {
  // We always use the highest item as the pivot.
  let pivot = a[index[high]]
  
  // This loop partitions the array into four (possibly empty) regions:
  //   [low  ...      i] contains all values <= pivot,
  //   [i+1  ...    j-1] contains all values > pivot,
  //   [j    ... high-1] are values we haven't looked at yet,
  //   [high           ] is the pivot value.
  var i = low
  for j in low..<high {
    if a[index[j]] <= pivot {
      (index[i], index[j]) = (index[j], index[i])
      
      i += 1
    }
  }
  
  // Swap the pivot element with the first element that is greater than
  // the pivot. Now the pivot sits between the <= and > regions and the
  // array is properly partitioned.
  (index[i], index[high]) = (index[high], index[i])
  
  return i
}

/*
 Uses a random pivot index. On average, this results in a well-balanced split
 of the input array.
 */
func quicksortRandom<T: Comparable>(_ index: inout [Int], using a: [T], low: Int, high: Int) {
  if low < high {
    // Create a random pivot index in the range [low...high].
    let pivotIndex = Int.random(in:low...high)

    
    // Because the Lomuto scheme expects a[high] to be the pivot entry, swap
    (index[pivotIndex], index[high]) = (index[high], index[pivotIndex])
    
    let p = partitionLomuto(&index, a, low: low, high: high)
    quicksortRandom(&index, using: a, low: low, high: p - 1)
    quicksortRandom(&index, using: a, low: p + 1, high: high)
  }
}

// Near zero
func is_near_zero(x: Double) -> Bool {
  return (x <= Static.Tolerance) && (x >= -Static.Tolerance)
}

// Nearly equal
func is_near(x: Double, y: Double) -> Bool {
  let z = x - y
  return (z <= Static.Tolerance) && (z >= -Static.Tolerance)
}
// Kahan and Babuska summation, Neumaier variant; accumulates less FP error
func sum(x: [Double]) ->  Double {
  var sum = 0.0
  var err = 0.0
  for k in x {
    let m = sum + k
    
    // Accumulate error separately
    err += abs(sum) >= abs(k) ? sum - m + k : k - m + sum
    sum = m
  }
  
  // Adjust final sum with accumulated error
  return sum + err
}

func circumCentre(first a:Point, second b:Point, third c:Point) -> Point {
  let ad = a.x * a.x + a.y * a.y
  let bd = b.x * b.x + b.y * b.y
  let cd = c.x * c.x + c.y * c.y
  let D = a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y)
  return Point(x: 0.5 / D * (ad * (b.y - c.y) + bd * (c.y - a.y) + cd * (a.y - b.y)),
               y: 0.5 / D * (ad * (c.x - b.x) + bd * (a.x - c.x) + cd * (b.x - a.x)))
}


/* Some Simple helper functions */
func nextHalfEdge(edge e:Int) -> Int { (e % 3 == 2) ? e - 2 : e + 1 }
func prevHalfEdge(edge e:Int) -> Int { (e % 3 == 0) ? e + 2 : e - 1 }


/* Triangle functions */
func edgesOf(triangle t: Int) -> Array<Int> { [3 * t, 3 * t + 1, 3 * t + 2] }
func triangleOf(edge e: Int) -> Int { (e / 3) } // e is an integer
func pointsOf(triangle t:Int, using triangles:Array<Int>) -> Array<Int> {edgesOf(triangle: t).map {triangles[$0]}}

func triangleCentre(triangle t:Int, using points:Array<Point>, using triangles:Array<Int>) -> Point  {
  let vertices:Array<Point> = pointsOf(triangle: t, using: triangles).map({points[$0]})
  return circumCentre(first:  vertices[0],
                      second: vertices[1],
                      third:  vertices[2])
}
