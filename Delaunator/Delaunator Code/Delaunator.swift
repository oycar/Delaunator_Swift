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
  static var Epsilon:Double = pow(2.0, -51)
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

// The main structure
struct Delaunator_Swift {
  var triangles: Array<Int>
  var halfEdges: Array<Int>
  var hull: Array<Int>
  
  var numberEdges: Int
  fileprivate var n, maxTriangles,  hashSize, hullStart, hullSize: Int
  fileprivate  var hullTri: Array<Int>
  fileprivate var hullPrev: Array<Int>
  fileprivate var centre: Point
  fileprivate  var coords = [Double]()
  
  // Use a provided init
  init(from points: Array<Point>) {
    n = points.count
    coords = Array(repeating:0.0, count:2*n)
    
    // Set the coords array
    for (i, p) in points.enumerated() {
      // Set each co-ordinate
      coords[2 * i]     = p.x
      coords[2 * i + 1] = p.y
    }
    
    // Update n
    n = coords.count >> 1
    
    // arrays that will store the triangulation graph
    maxTriangles = max(2 * n - 5, 0)
    numberEdges = 0
    triangles = Array(repeating:0, count:3 * maxTriangles)
    halfEdges = Array(repeating:0, count:3 * maxTriangles)
    
    hullTri  = Array(repeating:0, count:n) // edge to adjacent triangle
    hullPrev = Array(repeating:0, count:n) // edge to prev edge
    hull     = Array(repeating:0, count:n)
    
    hullTri  = Array(repeating:0, count:n) // edge to adjacent triangle
    hullPrev = Array(repeating:0, count:n) // edge to prev edge
    hull     = Array(repeating:0, count:n)
    
    hashSize = Int(Double(n).squareRoot().rounded(.up))
    hullStart = 0
    hullSize = 0
    centre = Point(x:0.0, y:0.0)
    
    // Back door
    if 0 == n {return}

    var hullNext = Array(repeating:0, count:n) // edge to next edge
    var hullHash = Array(repeating:-1, count:hashSize) // angular edge hash
    
    var dists = Array(repeating:0.0, count:n)
    var ids   = Array(repeating:0, count:n)
    
    // populate an array of point indices; calculate input data bbox
    var minX =  Double.infinity
    var maxX = -Double.infinity
    var minY =  Double.infinity
    var maxY = -Double.infinity
    
    // Compute bounds
    for i in 0..<n {
      let x = coords[2 * i]
      let y = coords[2 * i + 1]
      if (x < minX)  {minX = x}
      if (y < minY)  {minY = y}
      if (x > maxX)  {maxX = x}
      if (y > maxY)  {maxY = y}
      ids[i] = i
    }
    
    // Centre
    let c = Point(x: 0.5 * (minX + maxX), y: 0.5 * (minY + maxY))
        
    var minDist = Double.infinity
    var i0 = 0, i1 = 0, i2 = 0
    
    // pick a seed point close to the centre
    for i in 0..<n {
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
    for i in 0..<n {
      if i == i0 { continue }
      let d = dist(ax: i0x, ay: i0y, bx: coords[2 * i], by: coords[2 * i + 1])
      if (d < minDist && d > Static.Epsilon) {
        i1 = i
        minDist = d
      }
    }
    var i1x = coords[2 * i1]
    var i1y = coords[2 * i1 + 1]
    
    var minRadius = Double.infinity
    
    // find the third point which forms the smallest circumCircle with the first two
    for i in 0..<n {
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
    
    if minRadius == Double.infinity {
      // order collinear points by dx (or dy if all x are identical)
      // and return the list as a hull
      for i in 0..<n {
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
      quicksortRandom(&ids, using:dists, low:0, high:n - 1)
      var j = 0
      var d0 = -Double.infinity
      for i in 0..<n {
        let id = ids[i]
        if (dists[id] > d0) {
          hull[j] = id
          j += 1
          d0 = dists[id]
        }
      }
      hull.removeLast(n - j)
      triangles = []
      halfEdges = []
      return
    }
    
    // swap the order of the seed points for counter-clockwise orientation
    // These needs a lot of checking
    if (orient(rx: i0x, ry:i0y, qx: i1x, qy: i1y, px: i2x, py: i2y)) {
      let i = i1
      let x = i1x
      let y = i1y
      
      i1 = i2
      i1x = i2x
      i1y = i2y
      i2 = i
      i2x = x
      i2y = y
    }
    
    // Get the circumcentre of this triangle
    centre = circumCentre(first:  Point(x:i0x, y:i0y),
                          second: Point(x:i1x, y:i1y),
                          third:  Point(x:i2x, y:i2y))
    
    for i in 0..<n {
      dists[i] = dist(ax: coords[2 * i], ay:coords[2 * i + 1],
                      bx:centre.x, by:centre.y);
    }
    
    // sort the points by distance from the seed triangle circumCentre
    quicksortRandom(&ids, using:dists, low:0, high:n - 1)
    
    // set up the seed triangle as the starting hull
    hullStart = i0
    hullSize = 3
    
    hullNext[i0] = i1; hullPrev[i2] = i1
    hullNext[i1] = i2; hullPrev[i0] = i2
    hullNext[i2] = i0; hullPrev[i1] = i0
    
    hullTri[i0] = 0
    hullTri[i1] = 1
    hullTri[i2] = 2
    for i in 0..<hullHash.count {
      hullHash[i] = -1
    }
    
    hullHash[hashKey(p:Point(x: i0x, y: i0y))] = i0
    hullHash[hashKey(p:Point(x: i1x, y: i1y))] = i1
    hullHash[hashKey(p:Point(x: i2x, y: i2y))] = i2
    
    numberEdges = 0
    let _ = addTriangle(i0, i1, i2, -1, -1, -1)
    
    var xp = 0.0, yp = 0.0
    for k in 0..<ids.count {
      let i = ids[k]
      let x = coords[2 * i]
      let y = coords[2 * i + 1]
      
      // skip near-duplicate points
      if (k > 0 && is_near_zero(x:x-xp) && is_near_zero(x:y-yp)) {continue}
      xp = x
      yp = y
      
      // skip seed triangle points
      // This test exact equality?
      if (i == i0 || i == i1 || i == i2) {continue}
      
      // find a visible edge on the convex hull using edge hash
      var start = 0
      let key = hashKey(p:Point(x: x, y: y))
      for j in 0..<hashSize {
        start = hullHash[(key + j) % hashSize]
        if (start != -1 && start != hullNext[start]) {break}
      }
      
      start = hullPrev[start]
      var e = start
      var q = hullNext[e]
      while (!orient(rx: x, ry: y, qx: coords[2 * e], qy: coords[2 * e + 1], px: coords[2 * q], py: coords[2 * q + 1])) {
        e = q
        if (e == start) {
          e = -1
          break
        }
        q = hullNext[e]
      }
      
      if (e == -1) {continue} // likely a near-duplicate point; skip it
      
      // add the first triangle from the point
      var t = addTriangle(e, i, hullNext[e], -1, -1, hullTri[e])
      
      // recursively flip triangles from the point until they satisfy the Delaunay condition
      hullTri[i] = legalize(count: t + 2)
      hullTri[e] = t  // keep track of boundary triangles on the hull
      hullSize += 1
      
      // walk forward through the hull, adding more triangles and flipping recursively
      var next = hullNext[e]
      q = hullNext[next]
      
      while (orient(rx: x, ry: y, qx: coords[2 * next], qy: coords[2 * next + 1], px: coords[2 * q], py: coords[2 * q + 1])) {
        
        t = addTriangle(next, i, q, hullTri[i], -1, hullTri[next])
        hullTri[i] = legalize(count: t + 2)
        hullNext[next] = next // mark as removed
        hullSize -= 1
        next = q
        q = hullNext[next]
      }
      
      // walk backward from the other side, adding more triangles and flipping
      if (e == start) {
        q = hullPrev[e]
        while (orient(rx: x, ry: y, qx: coords[2 * q], qy: coords[2 * q + 1], px: coords[2 * e], py: coords[2 * e + 1])) {
          
          // Not sure this returns correct number of triangles
          t = addTriangle(q, i, e, -1, hullTri[e], hullTri[q])
          let _ = legalize(count:t + 2)
          hullTri[q] = t
          hullNext[e] = e // mark as removed
          hullSize -= 1
          e = q;
          q = hullPrev[e]
        }
      }
      
      // update the hull indices
      hullStart   = e; hullPrev[i]    = e
      hullNext[e] = i; hullPrev[next] = i
      hullNext[i] = next
      
      // save the two new edges in the hash table
      hullHash[hashKey(p:Point(x: x, y: y))] = i
      hullHash[hashKey(p:Point(x:coords[2 * e], y:coords[2 * e + 1]))] = e
    }
    
    hull = Array(repeating: 0, count: hullSize)
    var e = hullStart
    for i in 0..<hullSize {
      hull[i] = e
      e = hullNext[e]
    }
    
    // trim typed triangle mesh arrays ?????
    triangles.removeLast(triangles.count - numberEdges)
    halfEdges.removeLast(halfEdges.count - numberEdges)
  }
  
  mutating func legalize(count: Int) -> Int {
    var i = 0, ar = 0
    var a = count
    
    // recursion eliminated with a fixed-size stack
    while (true) {
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
       */
      let a0 = a - a % 3
      ar = a0 + (a + 2) % 3
      
      if (b == -1) { // convex hull edge
        if (i == 0) {break}
        i -= 1
        a = EDGE_STACK[i]
        continue
      }
      
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
        if (hbl == -1) {
          var e = hullStart
          repeat {
            if (hullTri[e] == bl) {
              hullTri[e] = a
              break
            }
            e = hullPrev[e]
          } while (e != hullStart);
        }
        
        link(a, hbl)
        link(b, halfEdges[ar])
        link(ar, bl)
        
        let br = b0 + (b + 1) % 3
        
        // don't worry about hitting the cap: it can only happen on extremely degenerate input
        if (i < EDGE_STACK.count) {
          EDGE_STACK[i] = br
          i += 1
        }
      } else {
        if (i == 0) {break}
        i -= 1
        a = EDGE_STACK[i]
      }
    }
    
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
    
    triangles[t] = i0
    triangles[t + 1] = i1
    triangles[t + 2] = i2
    
    link(t, a);
    link(t + 1, b)
    link(t + 2, c)
    
    self.numberEdges += 3
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
    let c = convex(rx:points[d.hull[j]].x, ry:points[d.hull[j]].y,
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

// return 2d orientation sign if we're confident in it through J. Shewchuk's error bound check
func orientIfSure(px: (Double), py: (Double),
                  rx: (Double), ry: (Double),
                  qx: (Double), qy: (Double)) -> Double {
  let l = (ry - py) * (qx - px);
  let r = (rx - px) * (qy - py);
  return abs(l - r) >= Static.Threshold * abs(l + r) ? l - r : 0.0
}

// a more robust orientation test that's stable in a given triangle (to fix robustness issues)
func orient(rx: (Double), ry: (Double),
            qx: (Double), qy: (Double),
            px: (Double), py: (Double)) -> Bool {
  // First test
  var sense = orientIfSure(px: px, py: py, rx: rx, ry: ry, qx: qx, qy: qy)
  if (sense != 0.0) {return sense < -Static.Tolerance}
  
  // Second test
  sense = orientIfSure(px: rx, py: ry, rx: qx, ry: qy, qx: px, qy: py)
  if (sense != 0.0) {return sense < -Static.Tolerance}
  
  // Final test
  sense = orientIfSure(px: qx, py: qy, rx: px, ry: py, qx: rx, qy: ry)
  
  // Final case
  return (sense < -Static.Tolerance)
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

// Is a hull convex?
func convex(rx: (Double), ry: (Double),
            qx: (Double), qy: (Double),
            px: (Double), py: (Double)) -> Bool {
  // First test
  var sense = orientIfSure(px: px, py: py, rx: rx, ry: ry, qx: qx, qy: qy)
  if !is_near_zero(x:sense) {return sense >= Static.Tolerance}
  
  // Second test
  sense = orientIfSure(px: rx, py: ry, rx: qx, ry: qy, qx: px, qy: py)
  if !is_near_zero(x:sense) {return sense >= Static.Tolerance}
  
  // Final test
  sense = orientIfSure(px: qx, py: qy, rx: px, ry: py, qx: rx, qy: ry)
  
  // Final case - different threshold
  return (sense >= -Static.Tolerance)
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
