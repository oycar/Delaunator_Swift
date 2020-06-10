//
//  Zone.swift
//  showMe
//
//  Created by Z Chameleon on 26/5/20.
//  Copyright © 2020 Z Cha. All rights reserved.
//

import SwiftUI

// Define a simple point struct
struct Point: Hashable, Codable {
  var x, y: Double
}

// Make point operate as a vector
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
  var id:UUID? = UUID()
  var isSelected: Bool = false
  
  // Variables which can be loaded from a header block
  var action:String? = "testValid"
  var tolerance:Double? = Static.Tolerance
  var comment:String?
  var hull: Array<Int>?
  var scale:Double? = 1
  
  init(name:String, triangulation: Triangulation) {
    self.name = name
    self.triangulation = triangulation
  }
}

// We want to write out a JSON data file
struct StoredZones : Codable {
  let zones: [Zone]
  
  struct Zone: Codable {
    var name: String
    var points: Array<Array<Double>>?
    
    var triangulation: Triangulation?
    
    // The header
    var header: Header?
  }
  
  struct Header: Hashable, Codable {
    // Extra stuff
    var action:String?
    var tolerance:Double = Static.Tolerance
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

  
  var image: Image {
    ImageStore.shared.image(name: name)
  }
}

struct Zone_Previews: PreviewProvider {
  static var previews: some View {
    Text("Hello, Zone World!")
  }
}



