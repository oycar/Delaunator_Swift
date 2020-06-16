//
//  Zone.swift
//  Delaunator_Swift
//
//  Created by Z Chameleon on 26/5/20.
//  Copyright © 2020 Z Cha. All rights reserved.
//
import Foundation


// Some extra methods for point
extension Point {
  var width: Double {
    get {self.x}
    set {self.x = newValue}
  }
  var height: Double {
    get {self.y}
    set {self.y = newValue}
  }
  init(width w: Double, height h: Double) {
    self.x = w
    self.y = h
  }
}

// This is way to complicated
struct Box: Hashable, Codable {
  var origin = Point(x:0, y:0)
  var size   = Point(width:0, height:0)
  
  // Lots of computed values
  var centre : Point {
    Point(x:origin.x + 0.5 * size.width,
          y:origin.y + 0.5 * size.height)
  }
  var height: Double {
    get {size.height}
    set {size.height = newValue}
  }
  var width:Double {
    get {size.width}
    set {size.width = newValue}
  }
  var scale:Double {max(size.width, size.height)}
  
  var left: Double {
    get {origin.x}
    set {origin.x = newValue}
  }
  var bottom: Double {
    get {origin.y}
    set {origin.y = newValue}
  }
  
  var right:  Double {origin.x + size.width}
  var top:    Double {origin.y + size.height}

}

struct Zone: Hashable, Codable, Identifiable {
  var name: String
  var triangulation: Triangulation? = Triangulation()
  var points: Array<Array<Double>>? = [[Double]]()
    
  // Extras for display code
  var box:Box?
  var id:UUID? = UUID()
  var isSelected: Bool = false
  
  // Variables which can be loaded from a header block
  var action:String? = "testValid"
  var tolerance:Double? = Static.Tolerance
  var comment:String?
  var hull: Array<Int>?
  var scale:Double? = 1
  
  init(name:String, triangulation: Triangulation, bounds box:Box?) {
    self.name = name
    self.triangulation = triangulation
    self.box = box
  }
}

// We want to write out a JSON data file
struct StoredZones : Codable {
  let zones: [Zone]
  
  struct Zone: Codable {
    var name: String
    var points: Array<Array<Double>>?
    
    // A triangulation could be read in too
    var triangulation: Triangulation?
    
    // The header
    var header: Header?
  }
  
  struct Header: Hashable, Codable {
    // Extra stuff
    var action:String?
    var tolerance:Double? = Static.Tolerance
    var comment:String?
    var hull: Array<Int>?
    var scale:Double?
  }
}

// The actual zone structure ...
extension Zone {
  init(from stored: StoredZones.Zone) {
    name = stored.name
    triangulation = stored.triangulation
    
    // Stuff from input
    points = stored.points
    
    action = stored.header?.action
    hull = stored.header?.hull
    scale = stored.header?.scale
    tolerance = stored.header?.tolerance
  }
}


// Process each zone
//
func readZones(using storedData:StoredZones) -> [Zone] {
  var testIt: ((String, Array<Point>, Delaunator_Swift?, Array<Int>?) -> Bool)? = nil
  var outputZones:[Zone] = [Zone]()
  
  // Get each zone - process one at a time
  var inputZones = [Zone]()
  for (_, z) in storedData.zones.enumerated() {
    inputZones.append(Zone(from: z))
  }
  
  // Output JSON
  let encoder = JSONEncoder()
  encoder.outputFormatting = .prettyPrinted
  
  let X = 0
  let Y = 1
  //let C = 2
  do {
    for (_, z) in inputZones.enumerated() {
      var list = [Point]()
      var box: Box?
      var phrase:String = z.name
      let perimeter:Array = z.points ?? [[Double]]()
      let name:String = z.name
      let scale:Double = z.scale ?? 1.0
      let expectedHull:Array = z.hull ?? [Int]()
      
      // Control which test (if any) is run
      // uses the header "action" tag
      switch z.action {
        case "testValid":
          testIt = validateTriangulation
          phrase = name + " Produces correct triangulation"
        case "testEquals":
          testIt = testEqual
          phrase = name + " Produces expected triangulation"
        default:
          testIt = logTriangulation
      }
      
      // Tolerance
      Static.Tolerance = z.tolerance ?? Static.Tolerance
      
      // populate an array of point indices; calculate input data bbox
      var minX =  Double.infinity
      var maxX = -Double.infinity
      var minY =  Double.infinity
      var maxY = -Double.infinity
      
      // Setup the list of points
      for (_, point) in perimeter.enumerated() {
        let p = Point(x:scale * point[X],
                      y:scale * point[Y])
        
        // Condition is not used ATM
        //condition = point[C]
        list.append(p)
        
        // Compute bounds
        if (p.x < minX)  {minX = p.x}
        if (p.y < minY)  {minY = p.y}
        if (p.x > maxX)  {maxX = p.x}
        if (p.y > maxY)  {maxY = p.y}
      }
      
      // Bounding box
      if perimeter.count > 0 {
        box = Box(origin:Point(x:minX, y:minY), size:Point(width: maxX - minX, height: maxY - minY))
      }
      
      // This is the triangulation step
      let delaunay = Delaunator_Swift(from: list)
      
      // Print half edges & triangles for comparison
      if let testFun = testIt {
        let _ = testFun(phrase, list, delaunay, expectedHull)
      }
      
      // Build up output list of zones
      outputZones.append(Zone(name: name, triangulation: Triangulation(using:delaunay, with:list), bounds: box))
    }
    
    let data = try encoder.encode(outputZones)
    print(String(data: data, encoding: .utf8)!)
    
  } catch {
    // Simplest possible error catching
    fatalError("Caught \(error)")
  }
  
  return outputZones
}

