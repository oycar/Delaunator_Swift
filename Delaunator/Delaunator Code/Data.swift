//
//  Data.swift
//  showMe
//
//  Created by Z Chameleon on 26/5/20.
//  Copyright © 2020 Z Cha. All rights reserved.
//


import UIKit
import SwiftUI
import CoreLocation

// This needs cleaning up
// Read input data using JSON
// What happens if it's given triangulated data?
// needs to recognize this... triangulated data... next step

let zoneData: [Zone] = readZones(using: load("zoneData.json"))
let triangulationData: [Triangulation] = zoneData.map({$0.triangulation ?? Triangulation()})

func load<T: Decodable>(_ filename: String) -> T {
  let data: Data
  
  guard let file = Bundle.main.url(forResource: filename, withExtension: nil)
    else {
      fatalError("Couldn't find \(filename) in main bundle.")
  }
  
  do {
    data = try Data(contentsOf: file)
  } catch {
    fatalError("Couldn't load \(filename) from main bundle:\n\(error)")
  }
  
  do {
    let decoder = JSONDecoder()
    return try decoder.decode(T.self, from: data)
  } catch {
    fatalError("Couldn't parse \(filename) as \(T.self):\n\(error)")
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
      
      // Setup the list of points
      for (_, point) in perimeter.enumerated() {
        let p = Point(x:scale * point[X],
                      y:scale * point[Y])
        
        // Condition is not used ATM
        //condition = point[C]
        list.append(p)
      }
      
      // This is the triangulation step
      var delaunay = Delaunator_Swift(from: list)
      delaunay.triangulate()
      
      // Print half edges & triangles for comparison
      if let testFun = testIt {
        let _ = testFun(phrase, list, delaunay, expectedHull)
      }
      
      // Build up output list of zones
      outputZones.append(Zone(name: name, triangulation: Triangulation(using:delaunay, with:list)))
    }
    
    let data = try encoder.encode(outputZones)
    print(String(data: data, encoding: .utf8)!)
    
  } catch {
    // Simplest possible error catching
    fatalError("Caught \(error)")
  }
  
  return outputZones
}


final class ImageStore {
  typealias _ImageDictionary = [String: CGImage]
  fileprivate var images: _ImageDictionary = [:]
  
  fileprivate static var scale = 2
  
  static var shared = ImageStore()
  
  func image(name: String) -> Image {
    let index = _guaranteeImage(name: name)
    
    return Image(images.values[index], scale: CGFloat(ImageStore.scale), label: Text(name))
  }
  
  static func loadImage(name: String) -> CGImage? {
    guard
      let url = Bundle.main.url(forResource: name, withExtension: "jpg"),
      let imageSource = CGImageSourceCreateWithURL(url as NSURL, nil),
      let image = CGImageSourceCreateImageAtIndex(imageSource, 0, nil)
      else {
        return nil
        //fatalError("Couldn't load image \(name).jpg from main bundle.")
    }
    return image
  }
  
  fileprivate func _guaranteeImage(name: String) -> _ImageDictionary.Index {
    if let index = images.index(forKey: name) { return index }
    
    images[name] = ImageStore.loadImage(name: name) ?? ImageStore.loadImage(name: "defaultDelaunay")
    //if let index = images.index(forKey: name) { return index }
    //return images.index(forKey: "ukraine")!
    return images.index(forKey: name)!
  }
}


struct Data_Previews: PreviewProvider {
  static var previews: some View {
    /*@START_MENU_TOKEN@*/Text("Hello, World!")/*@END_MENU_TOKEN@*/
  }
}
