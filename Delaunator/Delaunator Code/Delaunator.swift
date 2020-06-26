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
  static var Epsilon:Double = pow(2.0, -53)

  // Shewchuk error bounds
  static let orientThreshold:Double =  (3.0 + 16.0 * Epsilon) * Epsilon
  static let circleThreshold:Double = (10.0 + 96.0 * Epsilon) * Epsilon

}

// Write to the error stream
final class StandardErrorOutputStream: TextOutputStream {
  func write(_ string: String) {
    FileHandle.standardError.write(Data(string.utf8))
  }
}
var errorStream = StandardErrorOutputStream()

// Start with some global constants
var Edge_Stack:Array<Int> = []

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
    halfEdges   = delaunay.halfEdges
    hull        = delaunay.hull
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
  fileprivate var hashFactor: Double
  fileprivate var hullTri: Array<Int>
  fileprivate var hullPrev: Array<Int>
  fileprivate var centre: Point
  fileprivate var coords = [Double]()
  
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
    hashFactor = Double(numberPoints).squareRoot().rounded(.up)
    hashSize = Int(hashFactor)
    hashFactor *= 0.25
    
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
    
    var i0 = 0, i1 = 0, i2 = 0
    
    // pick a seed point close to the centre
    var minDist = Double.infinity
    for i in 0..<numberPoints {
      let d = dist(ax: c.x, ay: c.y, bx: coords[2 * i], by: coords[2 * i + 1])
      if d < minDist {
        i0 = i
        minDist = d
      }
    }
    let i0x = coords[2 * i0]
    let i0y = coords[2 * i0 + 1]
    
    // find the point closest to the seed
    minDist = Double.infinity
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
    
    // find the third point which forms the smallest circumCircle with the first two
    var minRadius = Double.infinity
    for i in 0..<numberPoints {
      if (i == i0) || (i == i1) { continue }
      let r = circumRadius(i0x, i0y,
                           i1x, i1y,
                           coords[2 * i], coords[2 * i + 1])
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
      
      // This sorts the indices in place
      // Will either be the x or y component of the actual separations
      quicksortRandom(&ids, using:dists, low:0, high:numberPoints - 1)
      
      // The hull is just the sorted ids
      hull = ids
      hull.removeLast(numberPoints - hull.count)
      triangles = []
      halfEdges = []
      return
    }
    
    // swap the order of the seed points
    // to ensure seed triangle has
    // an anticlockwise orientation
    if (!orient(i0x, i0y, i1x, i1y, i2x, i2y)) {
      // Swap i1 & i2
      (i1, i2) = (i2, i1)
      (i1x, i2x) = (i2x, i1x)
      (i1y, i2y) = (i2y, i1y)
    }
    
    // All the work above was to initialize the seed triangle
    // The circum-centre of this traingle is the origin for the
    // radial and azimuthal sorts
    
    // Get the circumcentre of this triangle
    centre = circumCentre(i0x, i0y, i1x, i1y, i2x, i2y)
    
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
    // The first trianle is the convex hull so all linked edges ar (-1)
    let _ = addTriangle(i0, i1, i2, -1, -1, -1)
    
    // Now start building the triangulation
    // The previous point
    var xp:Double = 0, yp:Double = 0

    // Labels! Just like good old FORTRAN IV
    projectPoints: for (k, i) in ids.enumerated() {
      let x = coords[2 * i]
      let y = coords[2 * i + 1]
      
      // skip near-duplicate points
      if (k > 0 && is_near_zero(x:x-xp) && is_near_zero(x:y-yp)) {continue projectPoints}
      
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
      while (!orient(x, y,
                     coords[2 * q], coords[2 * q + 1],
                     coords[2 * e], coords[2 * e + 1])) {
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
      hullTri[i] = legalize(edge: t + 2)
      hullTri[e] = t  // keep track of boundary triangles on the hull
      
      // We added one triangle - so the hull size goes up by one
      //
      //     i
      //    / \
      //   q---e        The hull edge (e->q) has been replaced by two new edges (e->i) and (i->q)
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
      while (orient(x, y,
                    coords[2 * q], coords[2 * q + 1],
                    coords[2 * n], coords[2 * n + 1])) {

        // Add this triangle to the triangulation
        t = addTriangle(n, i, q, hullTri[i], -1, hullTri[n])
        
        // Fix locally non-delaunay edges
        // Returns hull edge i on return (the new edge)
        hullTri[i] = legalize(edge: t + 2)

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
        while (orient(x, y,
                      coords[2 * e], coords[2 * e + 1],
                      coords[2 * q], coords[2 * q + 1])) {
          
          // Add this triangle to the triangulation
          t = addTriangle(q, i, e, -1, hullTri[e], hullTri[q])
          
          // Fix locally non-delaunay edges
          // Ignore return value (which is edge i - not on the hull)
          let _ = legalize(edge:t + 2)
          
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
    
    // Record the hull
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
  
  
  private mutating func legalize(edge e: Int) -> Int {
    var a  = e,
        a2 = 0
    
    // recursion eliminated with a fixed-size stack
    // Maintains a stack of edges for flipping
    // If an edge needs flipping adds it to the stack - retains edge
    
    flipEdge: while (true) {
      // A subtlety here is that if an edge [a] is flipped
      // then no new edge is popped from the stack so the edge [a] (the old [b2])
      // is tested again - so in fact two edges ([a] and [b1]) are checked
      let b:Int = halfEdges[a]
      
      /* if the pair of triangles doesn't satisfy the Delaunay condition flip it
       *
       *    Input edge is a
       *
       *
       *            n                     n
       *          /| |\                  /  \
       *      a1 / | | \b2            a1/    \a
       *        /  | |  \              /      \
       *       /   | |   \    flip    /___a2___\
       *      i   a| |b  p     =>    i__________p
       *       \   | |   /            \   b2   /
       *        \  | |  /              \      /
       *       a2\ | | /b1             b\    /b1
       *          \| |/                  \  /
       *            q                     q
       *
       *
       */
      
      // Old triangles are [a] triangle => [a, a1, a2] halfEdges => [q, n, i]
      //               and [b] triangle => [b, b1, b2] halfEdges => [n, q, p]

      // New triangles are [a] triangle => [a, a1, a2] halfEdges => [p, n, i]
      //               and [b] triangle => [b, b1, b2] halfEdges => [i, q, p]
      //
      // So changes    are triangles[a]  => p
      //                   triangles[b]  => i
      //                   halfEdges[a]  <=> halfEdges[b2]
      //                   halfEdges[a2] <=> halfEdges[b2]
      //                   halfEdges[b]  <=> halfEdges[a2]

      // No way to know which of the [0, 1, 2] edges of b were cut in the projection
      // so complex remainder trickery to get the edge ids
      
      // Get the edges associated with each triangle
      // Triangle a - need a2 for special case
      let a0 = a - a % 3
      a2 = a0 + (a + 2) % 3
      
      // First case is if (a <-> b) is the convex hell edge
      // In this case no flip can occur
      if (b == -1) {
        // If no edges left on stack all edges are locally delaunay
        // Get the next pending edge
        a = Edge_Stack.popLast() ?? -1
        
        // All done if this is negative
        if (a == -1) {break flipEdge}
        continue flipEdge
      }
      
      // Finish Triangle a & Triangle b - don't actually need b1 unless added to the stack
      let a1 = a0 + (a + 1) % 3
      let b0 = b - b % 3
      let b2 = b0 + (b + 2) % 3

      // Now the four points
      let n = triangles[a1]
      let i = triangles[a2]
      let q = triangles[a]
      let p = triangles[b2]

      // We need to flip the edges a-b if the point p is inside the circum-circle
      // of the triangle [n, i, q]

      // Use Shewchuk's robust version of inCircumCircle
      let illegal = inCircumCircle(
        ax: coords[2 * n], ay: coords[2 * n + 1],
        bx: coords[2 * i], by: coords[2 * i + 1],
        cx: coords[2 * q], cy: coords[2 * q + 1],
        px: coords[2 * p], py: coords[2 * p + 1])
      
      // We have to flip this edge
      if (illegal) {
        triangles[a] = p
        triangles[b] = i
        
        // edge swapped on the other side of the hull (rare); fix the halfEdge reference
        // Repairing the hull when a flip has disturbed it
        // _b2_ is the only edge this can happen with - because _b1_ is not swapped,
        // and _a1_ and _a2_ are on the hull already
        let hb2 = halfEdges[b2]
        if (hb2 == -1) {
          // Repair the hull
          // Simply patch in the edge a at point p
          hullTri[p] = a
        }
        
        // Link the half edges
        //   halfEdges[a]  <=> halfEdges[b2]
        //   halfEdges[a2] <=> halfEdges[b2]
        //   halfEdges[b]  <=> halfEdges[a2]
        link(a, hb2)
        link(b, halfEdges[a2])
        link(a2, b2)
                
        // Add the edge [b1] to the stack
        // The edge [a] is not popped so will be tested again
        // What about [a1] and [b] ? If they are on the hull
        // this is ok - but after some flips is this still true?
        Edge_Stack.append(b0 + (b + 1) % 3) // b1
      } else {
        // Get here when an edge (a) is locally delaunay
        
        // Pop the last element
        a = Edge_Stack.popLast() ?? -1
        if (a == -1) {break flipEdge}
      }
    }
    
    // Only gets here when no flip - because stack is not empty after a flip
    return a2
  }
  
  // Get a hash key - slightly simplified
  private func hashKey(p:Point) -> Int {
    return Int(hashFactor * pseudoAngle(dx:p.x - centre.x, dy:p.y - centre.y).rounded(.down)) % hashSize
  }
  
  
  private mutating func link(_ a: Int, _ b: Int) {
    halfEdges[a] = b
    //logTriangles(e: a)
    
    if (b != -1) {
      halfEdges[b] = a
      //logTriangles(e: b)
    }
    
  }
  
  // add a new triangle given vertex indices and adjacent half-edge ids
  private mutating func addTriangle(_ i0: Int, _ i1: Int, _ i2: Int,
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
  
  
  // check if an edge is locally Delaunay
  func isDelaunay(edge a:Int) -> Bool {
    // The opposing halfEdge
    let b:Int = halfEdges[a]
    
    // If this is on the convex hull it is locally Delaunay
    // Also only consider internal edges for which a > b
    // (to prevent testing the same four four points twice)
    if (a > b) {
      // This will always be true when e is on the hull since b == -1
      return true
    }
    
    // Get the edges associated with each triangle
    // Triangle a
    let a0 = a - a % 3
    let a1 = a0 + (a + 1) % 3
    let a2 = a0 + (a + 2) % 3
    
    // Triangle b
    let b0 = b - b % 3
    let b2 = b0 + (b + 2) % 3

    /* if the pair of triangles doesn't satisfy the Delaunay condition return false
     *
     *    Input edge is a
     *
     *            n
     *          /| |\
     *      a1 / | | \b2
     *        /  | |  \
     *       /   | |   \
     *      i   a| |b   p
     *       \   | |   /
     *        \  | |  /
     *       a2\ | | /b1
     *          \| |/
     *            q
     *
     */
    
    // Now need the four points making up the quadrilateral
    // Now the four points
    let n = triangles[a1]
    let i = triangles[a2]
    let q = triangles[a]
    let p = triangles[b2]
    
    // two possible tests - rounding error can cause an endless cycle of flips
    // This is overkill it seems
    let firstTest = inCircumCircle(
      ax: coords[2 * n], ay: coords[2 * n + 1],
      bx: coords[2 * i], by: coords[2 * i + 1],
      cx: coords[2 * q], cy: coords[2 * q + 1],
      px: coords[2 * p], py: coords[2 * p + 1])
    let otherTest = inCircumCircle(
      ax: coords[2 * p], ay: coords[2 * p + 1],
      bx: coords[2 * n], by: coords[2 * n + 1],
      cx: coords[2 * i], cy: coords[2 * i + 1],
      px: coords[2 * q], py: coords[2 * q + 1])
    
    // If tests agree - don't flip otherwise use first
    // because the points are essentially on the circumcircle
    // and therefore no real difference in which edge is used
    return otherTest == firstTest ? true : !firstTest
  }
}

/* Helper Functions
 dist:
 inCircumCircle:
 circumCentre:
 circumRadius:
 
 */
func dist(ax: (Double), ay: (Double),
          bx: (Double), by: (Double)) -> Double {
  let dx = ax - bx
  let dy = ay - by
  return dx * dx + dy * dy
}

// Is point p inside the triangle abc
// Documentation - Shewchuk
func inCircumCircle(ax: (Double), ay: (Double),
              bx: (Double), by: (Double),
              cx: (Double), cy: (Double),
              px: (Double), py: (Double)) -> Bool {
  let adx = ax - px
  let ady = ay - py
  let bdx = bx - px
  let bdy = by - py
  let cdx = cx - px
  let cdy = cy - py
  
  let bdxcdy = bdx * cdy;
  let cdxbdy = cdx * bdy;
  let alift = adx * adx + ady * ady;
  
  let cdxady = cdx * ady;
  let adxcdy = adx * cdy;
  let blift = bdx * bdx + bdy * bdy;
  
  let adxbdy = adx * bdy;
  let bdxady = bdx * ady;
  let clift = cdx * cdx + cdy * cdy;
  
  let det = alift * (bdxcdy - cdxbdy)
    + blift * (cdxady - adxcdy)
    + clift * (adxbdy - bdxady);
  
  let permanent = (abs(bdxcdy) + abs(cdxbdy)) * alift
    + (abs(cdxady) + abs(adxcdy)) * blift
    + (abs(adxbdy) + abs(bdxady)) * clift;
  let errbound = Static.circleThreshold * permanent
  return  det < -errbound
}

// Get the circum-radius
func circumRadius(_ ax: Double, _ ay: Double,
                  _ bx: Double, _ by: Double,
                  _ cx: Double, _ cy: Double) -> Double {
  // Points a, b, c
  // Form a triangle
  
  let dx = bx - ax
  let dy = by - ay
  let ex = cx - ax
  let ey = cy - ay
  
  let bl = dx * dx + dy * dy
  let cl = ex * ex + ey * ey
  let d = 0.5 / (dx * ey - dy * ex)
 
  // Is zero present in this list?
  if ([bl, cl, d].contains {is_near_zero(x: $0)} ) {
    return 0.0
  }
  
  // Compute r squared
  let x = (ey * bl - dy * cl) * d
  let y = (dx * cl - ex * bl) * d
  
  return x * x + y * y // r
}

// Get the centre of the circumcircle
func circumCentre(_ ax: (Double), _ ay: (Double),
                  _ bx: (Double), _ by: (Double),
                  _ cx: (Double), _ cy: (Double)) -> Point {
  let ad = ax * ax + ay * ay
  let bd = bx * bx + by * by
  let cd = cx * cx + cy * cy
  let D = ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)
  return Point(x: 0.5 / D * (ad * (by - cy) + bd * (cy - ay) + cd * (ay - by)),
               y: 0.5 / D * (ad * (cx - bx) + bd * (ax - cx) + cd * (bx - ax)))
}

// monotonically increases with real angle, but doesn't need expensive trigonometry
func pseudoAngle(dx: Double, dy: Double) -> Double {
  let p = dx / (abs(dx) + abs(dy))
  return dy > 0 ? 3 - p : 1 + p // [0..4]
}

// return 2d orientation sign if we're confident in it through J. Shewchuk's error bound check
//
// Is a point c to the left or right of a vector a->b
// Or more symmetrically are the points (a->b->c) in an anti-clockwise order? +ve yes -ve no
// The magnitude of this quantity is twice the triangle area
//
/**
 [Robust Predicates](http://www.cs.cmu.edu/~quake/robust.html)
 */
func orientIfSure(_ ax: (Double), _ ay: (Double),
                  _ bx: (Double), _ by: (Double),
                  _ cx: (Double), _ cy: (Double)) -> Double {
  
  // Points a, b, c
  // Form a triangle
  
  // The left and right determinants
  let detLeft  = (ay - cy) * (bx - cx)
  let detRight = (ax - cx) * (by - cy)
  
  // We only need to determine the sign +ve/zero/-ve
  // Positive - anticlockwise
  // Zero     - colinear
  // Negative - clockwise
  let det = detLeft - detRight
  return abs(det) >= Static.orientThreshold * abs(detLeft + detRight) ? det : 0
}

// a more robust orientation test that's stable in a given triangle (to fix robustness issues)
// This returns true  (anticlockwise)
//              false (clockwise)
func orient(_ ax: (Double), _ ay: (Double),
            _ bx: (Double), _ by: (Double),
            _ cx: (Double), _ cy: (Double)) -> Bool {
  // First test triangle [a, b, c]
  var sense = orientIfSure(ax, ay,
                           bx, by,
                           cx, cy)
  if (sense != 0) {return sense > 0}
  
  // Second test triangle [b, c, a]
  sense = orientIfSure(bx, by,
                       cx, cy,
                       ax, ay)
  if (sense != 0) {return sense > 0}
  
  // Final test triangle [c, a, b]
  // If sense == 0 it was so for all the tests - so in this case assume true if zero
  sense = orientIfSure(cx, cy,
                       ax, ay,
                       bx, by)
  
  // Final case
  return (sense > 0)
}

// Get the area of a triangle [a, b, c]
func getTriangleArea(_ ax: (Double), _ ay: (Double),
                     _ bx: (Double), _ by: (Double),
                     _ cx: (Double), _ cy: (Double)) -> Double { 0.5 * orientIfSure(ax, ay, bx, by, cx, cy) }


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
  return (x <= Static.Epsilon) && (x >= -Static.Epsilon)
}

// Nearly equal
func is_near(x: Double, y: Double) -> Bool {
  let z = x - y
  return (z <= Static.Epsilon) && (z >= -Static.Epsilon)
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

/* Some Simple helper functions */
func nextHalfEdge(edge e:Int) -> Int { (e % 3 == 2) ? e - 2 : e + 1 }
func prevHalfEdge(edge e:Int) -> Int { (e % 3 == 0) ? e + 2 : e - 1 }

/* Triangle functions */
func edgesOf(triangle t: Int) -> Array<Int> { [3 * t, 3 * t + 1, 3 * t + 2] }
func triangleOf(edge e: Int) -> Int { (e / 3) } // e is an integer
func pointsOf(triangle t:Int, using triangles:Array<Int>) -> Array<Int> {edgesOf(triangle: t).map {triangles[$0]}}

func triangleCentre(triangle t:Int, using points:Array<Point>, using triangles:Array<Int>) -> Point  {
  let vertices:Array<Point> = pointsOf(triangle: t, using: triangles).map({points[$0]})
  return circumCentre(vertices[0].x, vertices[0].y,
                      vertices[1].x, vertices[1].y,
                      vertices[2].x, vertices[2].y)
}
